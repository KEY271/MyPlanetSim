# Phase 14 — RHS 評価削減・物理の頻度分離・CPU 性能改善

状態: 実装中（2026-09-09）。P14.01 の計測基盤・集計境界修正・発達後checkpoint比較を実装。
P14.02 のN=6湿潤・固定ICI反復行列を実装・測定済み。独立 ARK2 IMEX 比較器を実装し、
モデル結合比較と P14.03 以後は未着手。
詳細は [検証記録](validation/phase-14.md)。
調査基準: `b6d3bbb`。以下の新しい時間刻み・性能目標・誤差閾値は提案値であり、実測値ではない。

## 1. 到達点と優先順位

Phase 14 は、新しい物理を追加せず、Phase 13 の湿潤モデルを同じ精度で安く実行することに集中する。
順序は **計測 → RHS 評価回数削減 → 物理の頻度分離 → column/RHS 並列化 →
reconstruction–flux 融合 → 配列・modal 演算・tile 最適化 → 統合検証** とする。
実装は一項目ずつ評価し、単独で速くても全体を遅くする変更は標準経路へ採用しない。

標準候補は、力学 `dt=1800 s`、ICI 2反復、BL/surface `600 s`、SBM 診断 `900 s`、
放射診断 `1800 s`。BL `900 s`、SBM `1800 s` は次の比較候補とする。
これらは最終既定値の宣言ではなく、下記の検証で採否を決める実験点である。
飽和調節は温湿度が変化する物理更新の後に実行し、SBM と一緒に低頻度化しない。

MPI、GPU backend、並列 I/O、非静水圧化、水蒸気の放射効果・雲・氷は Phase 15 以後の候補。
GPU 化に利用できる view と batched kernel の境界は本 Phase で用意する。

## 2. 現状から修正すべき前提

knowledge MCP を `MyPlanetSim` で検索したが該当記録はなかった。
判断の根拠は現在のコードと [Phase 13 検証報告](validation/phase-13.md) である。

| 項目 | 現在の実装・実測 | 計画への反映 |
|---|---|---|
| ICI | `dry_hydrostatic_driver.cpp` が反復と最終候補の full RHS を評価。湿潤の後処理で次 step への RHS 再利用が無効になる | 反復数と実際の full RHS 評価数を別々に数える |
| 放射 | `rhs_with_components()` 内。RK stage / ICI candidate ごとに実行され、observer の診断にも別呼出しがある | 300秒の物理ループだけを分けても放射は安くならない。RHS からの分離を独立変更にする |
| 湿潤物理 | `apply_local_physics()` が BL → 乾燥調節 → SBM → 飽和調節を共通 interval で実行。各段階で全球 diagnose/validate、retry 用 state コピーがある | 過程ごとの診断・更新・再検査を分け、後に column 内へ局所化する |
| BL/surface | `implicit_boundary_layer_column()` は既に backward Euler の三重対角 solve。蒸発・潜熱も応答解と bracket 付き二分法で結合 | 主課題は係数固定の時間精度、非線形 solve のコスト、再計算の重複。単なる陽→陰置換ではない |
| OpenMP | reconstruction の level ループのみ、最大4 thread。`MPS_ENABLE_OPENMP` は既定 OFF | edge scatter、共有 workspace、診断加算を整理して範囲を広げる |
| 配列 | scalar は `cell*K+level`、tracer は `(tracer*C+cell)*K+level`。水平演算に level slice のコピーがある | column の連続性も評価しながら blocked/AoSoA を比較する |
| face/flux | 全球 face state を保存し、driver が再走査。tracer flux も左右セルへの加算時に再計算される | 共有 edge で一度だけ flux を計算し、両セルへ同じ値を渡す |
| allocation | moist coupling に cell ごとの `std::vector<Real> vapor(levels)` 等が残る。Phase 13 の300秒 substep 比較では warm allocation/step=1366 | thread ごとの再利用 workspace を先に整える。プロセス全体と kernel の allocation を区別する |

Phase 13 の1200日 pilot は N=6、K=20、上端20,000 Pa、地表細分化 S=5、2 tracer、
力学1800秒・ICI 3反復・物理300秒で、52分13秒、`2.611 s/model-day`、retry=0。
これは過去の実測参照であり、新しい build/CPU の性能比較の分母には同一環境で再測定した baseline を使う。
上端で深い対流が切られる問題と約 `0.3 W m-2` の既存全系エネルギー残差を引き継いで記録する。

## 3. P14.01 — 計測と比較基準の固定

既存の `benchmark_dry_core`、`benchmark_semi_implicit`、`benchmark_moist` を拡張し、
共通 runner `tools/benchmark_phase14.py` と小さい JSON 集計を追加する計画とする。
計測runnerと精度manifestは実装済み。以下の後続項目の新規API・configは、
検証記録に実装済みと明記したものを除き、計画段階である。

### 計測項目

- full RHS を initial / iteration / final / retry / diagnostic に分けて計数。
  総評価数/受理step、総評価数/model-day、破棄試行のコストも出す。
- RHS の排他的時間: diagnose、level pack、gradient/limiter、face reconstruction、flux、
  column coupling、source、diffusion、RHS 内 physics。modal transform、fast operator、
  Helmholtz/GMRES、物理の外側の時間も分ける。
- physics の過程別 column call、重い診断回数、安い increment 適用回数、実更新間隔、
  solve iteration、cache age、retry、state copy の時間と bytes。
- wall/model-day、accepted dt の分布、cell-level updates/s、peak RSS、scratch bytes/thread、
  warm allocation、推定 read/write bytes。使用可能な環境では bandwidth/cache miss と
  compiler の vectorization report を取得し、取得不可は欠測と記す。
- timer の親子関係を定義し、RHS 内 radiation と RHS 全体などを加算して二重計上しない。
  timer 無効の対照実行で計測 overhead も確認する。

### 比較行列

| 軸 | 最小行列 |
|---|---|
| 湿潤・乾燥 | Phase 10 Held–Suarez、Phase 13 aquaplanet、有限陸 fraction＋地形 |
| 水平 | N=6/12（小格子と既存基準）、24/48（並列・cache）、96（短い高解像度試験） |
| 鉛直 | K=20/40、K=80 は modal/layout の短い kernel 試験から開始 |
| tracer | 1/2/4。registry 上の水の位置も入れ替える |
| thread | OpenMP OFF、ON の1/2/4、搭載 CPU に応じて8以上 |
| 大気状態 | 一様初期場と、SBM が活動している発達後の同一 checkpoint |

全組合せを長期積分しない。kernel は代表 snapshot、短い統合比較は N=6/12/48、
30日・1200日は登録済み N=6/K=20 を中心に段階実行する。
高解像度では advective/vertical CFL に従って dt を縮め、1800秒を強制しない。
SBM が不活性な初日だけで湿潤の性能を判定しない。

同じ compiler/Release flags/FP設定、同じモデル期間、同じ初期場、同じ出力頻度で測定する。
warm-up 後の短い run は最低5回、中央値と範囲を報告し、CPU・thread配置・電源条件・
OpenMP runtime・commit・dirty flag・config fingerprint を保存する。
allocation 計測は runtime の初期化を除いた別 run とし、計数自体で速度を変えない。
OpenMP 導入効果は OFF baseline と、同じ新実装の ON/1 thread の両方を分母として出す。

## 4. P14.02 — Semi-implicit の full RHS 評価回数削減

変更箇所: `dry_hydrostatic_driver.cpp`、`dry_hydrostatic_semi_implicit.cpp`、
semi-implicit config・diagnostics、Phase 10/13 の検証 runner。

### 最初の候補: 固定2反復

既存 `semi_implicit.nonlinear_iterations` で 2/3/5 を比較する。最初は数値コードを変えず、
dt=450/900/1800秒、物理300秒を固定した行列で2反復の使用可能範囲を登録する。

retry・observer 呼出しを除く現在の full RHS 回数は次のとおり。

| 経路 | ICI 3反復 | ICI 2反復 | 理想的な RHS 部分の時間削減 |
|---|---:|---:|---:|
| 湿潤＋gray radiation（毎 step initial を再評価） | 5 | 4 | 20% |
| 時刻非依存 RHS、後処理なし、再利用が効く定常 step | 3 | 2 | 33% |
| 同乾燥経路の初回 step | 4 | 3 | 25% |

線形 solve も3回から2回になるが、これをそのまま全体の1.5倍速とは見積もらない。
最終評価は CFL、energy attribution、残差、次stepへの再利用を担うため、先に削除しない。
物理によって state が変わった後に、物理前の `candidate_rhs` を再利用しない。
時刻依存 RHS では、state 配列が同じでも評価時刻が違えば初回候補と initial を同一視しない。

線形残差許容値・発散検査・正値性・CFL は維持する。2反復が不合格ならケースを3反復に戻す。
反復削減を成立させるためだけの diffusion 増強や残差閾値の緩和は採用しない。
新しい predictor を用いる場合は、初期候補での RHS 再利用条件も更新する。

### 第二候補: predictor-corrector / IMEX

2反復の結果を基準に、`N(y,t)=F_dyn(y,t)-L(y-y_ref)` と固定 reference-linear fast operator
を使う二次の predictor-corrector を試作する。独立比較器として既知の二次 additive RK を用意し、
Butcher table、stage 時刻、明示 RHS/implicit solve 回数を ADR に固定する。
方式の比較には [SUNDIALS ARKODE の公式 Butcher table](https://sundials.readthedocs.io/en/v7.7.0/arkode/Butcher_link.html)
を用いる。SUNDIALS 自体の runtime 依存を導入する計画ではない。

目標は **重い nonlinear RHS が通常2回、最終検査を含めても3回以下/受理step**。
`N` の評価には `F_dyn` が必要なので、fast operator の呼出しに名前を替えて評価数を隠さない。
final full RHS を省く候補は、安い state/CFL 検査、stage weight による過程収支、solver 検査を
先に実装する。debug では full RHS と照合し、放射や後処理の変更とは別に比較する。

linear gravity wave の位相・振幅・減衰、fast/slow split の加減算一致、dt 半減時の二次収束、
terrain rest、mountain wave、baroclinic、鋭い tracer と水の非負性を通す。
implicit stage の安定性だけから、全体系の大刻み安定性や tracer の SSP 性を主張しない。
分割物理を含む全体の時間次数は別に測定する。
精度を満たす2反復より end-to-end で10%以上速くならなければ、IMEX は標準採用しない。

## 5. P14.03 — Physics を過程ごとの時間間隔へ分離

変更箇所: driver の `apply_local_physics`、physics coupling、config、checkpoint、診断集計。
`PhysicsSchedule`、`PhysicsWorkspace`、受理済み過程 increment の構造を導入する。
**重い係数・参照場の診断周期と、保存的な state 更新周期を別に表す。**

| 過程 | 最初の比較 | 更新の扱い |
|---|---|---|
| BL＋surface＋蒸発 | 300 → 600 → 900秒 | 大気・地表・bucket を同じ coupled solve で更新 |
| SBM の parcel/CAPE/LCL/参照場 | 300 → 900 → 1800秒 | 全体を間欠更新する比較と、参照場を保持して緩和だけ小刻みに行う比較を分ける |
| 放射 | 現行 RHS 内 → 独立診断900/1800秒 | 診断した interface flux を安い更新へ供給し、大気と地表へ同一 flux を使う |
| 乾燥調節 | 温度更新後の安定性 scan | 必要な柱だけ調節。SBM 診断前にも乾燥不安定を除く |
| 飽和調節 | 力学後・温湿度を変える物理更新後 | 軽い過飽和判定、該当層で solve。SBM off の時も継続 |

### スケジュールと状態更新の契約

1. 最初に legacy schedule を同じ呼出し順で実装し、300秒設定の state/checkpoint 一致を確認する。
   旧 config はこのモードに対応させ、新しい間隔を指定したときだけ multirate を選択する。
2. 初版の間欠更新は、受理候補の力学区間を各過程の最大 interval で分割する。
   例えば力学1800秒で BL600秒は3回、SBM900秒は2回。過程の終端イベントを時刻順に処理し、
   同時刻では放射 increment → BL/surface → 乾燥調節 → SBM → 飽和調節 → 降水移送を基本とする。
   飽和調節は中間の温湿度更新にも必要に応じて挿入し、この分割順を ADR に明示する。
3. 各過程は最後の適用からの実 interval を使い、力学区間末で端数を処理する。
   力学 dt が450秒なら、900秒指定でも間欠更新は450秒以下となる。これは「最大 interval」であり、
   呼出し削減を保証する周期ではない。stage 回数や `step % n` ではスケジュールしない。
4. 力学 step をまたいで重い診断だけを省く cache 経路は、診断時刻・最大 age・参照場・有効性を保持する。
   state の増分適用は毎区間で完結させ、降水や潜熱を将来のイベントまで未計上にしない。
   放射の太陽位置は物理的な評価時刻から求め、夜明け/日没などでは短波を早期更新する。
5. pressure/mass が変わる step 境界をまたいで `K/s` の加熱率をそのまま流用しない。
   放射は `W/m2` の共有 interface flux を持ち、現在の層質量・Exner で保存量増分へ変換する。
   cached flux を使う区間でも温度・地表変化の制約を検査し、必要なら再診断・細分化する。
6. 物理移動後は RHS から該当 tendency と CFL 制約を除き、scheduler 側に制約と収支を移す。
   放射を二重に加えない。RHS 内積分から splitting への変更自体を、周期延長より先に検証する。
7. 失敗時は state、cache、時計、bucket、累積水/熱台帳を一括 rollback する。
   採用 increment だけを収支に反映し、破棄試行の回数・時間は性能診断へ残す。
8. checkpoint に cache と時計を含める場合は schema を version 化する。fingerprint に周期・方式を含め、
   旧 checkpoint は legacy 経路で読む。restart は周期境界以外、dt 短縮、retry の直後でも検証する。
   observer が診断を再計算しても、積分用 cache と過程時計を書き換えない。

最初は現在と同じ dynamics → physics の分割から始める。Strang 化は精度不足が分割順に由来すると
分かった場合の独立候補とし、半step化で診断呼出しが倍増しないことも測る。
新 config の候補名は `physics.schedule`、`boundary_layer.maximum_update_interval_s`、
`convection.diagnostic_interval_s`、`convection.update_mode`、`radiation.diagnostic_interval_s`。
旧 `moisture.maximum_physics_substep_s` との併用・優先順位は曖昧にせず、競合指定を入力時に検出する。

## 6. P14.04 — Simple Betts–Miller の診断頻度を下げる

変更箇所: `simple_betts_miller.cpp/.hpp`、`moist_physics_coupling.cpp/.hpp`。
parcel ascent と参照場診断を `diagnose_sbm_reference()`、有限増分を `apply_sbm_relaxation()`
相当の責務へ分ける。まず両者を毎回呼んで旧結果との一致を確認する。

比較は次の順で行う。

1. 現行 backward Euler 係数を維持し、SBM 全体を300/900/1800秒の間欠更新にする。
   BL・放射の頻度は固定し、時間刻みの効果だけを測る。
2. CAPE/LCL/parcel/支持域を900/1800秒保持し、温湿度への緩和を300秒で適用する。
   毎回、現在の層質量・温湿度で enthalpy 補正と水制約を再確認する。
   古い補正済み `Tref/qref` を無条件に再利用しない。浅い対流の柱水保存が成立しなくなれば再診断する。
3. 固定参照場・支持率 `f` での `a_k=1-exp(-f_k*dt/tau)` による analytic 緩和を別途比較する。
   現行 `a*f` と同じ演算ではない。層ごとに係数が異なる場合には有限増分の保存条件を再導出し、
   温湿度を独立に clip して保存を壊さない。

cache は最大 age に加え、最下層温度・湿度・pressure の変化、参照場の不適合で無効化する。
開始候補の変化閾値は `|Delta T_bottom|>0.5 K`、
`|Delta q_bottom|>max(0.1*q_bottom, 1e-5)`、`|Delta ps|/ps>0.01`。
閾値を1/2倍・2倍して診断回数と誤差を測り、値を登録する。
「CAPE が増えたら再診断」のために CAPE を毎回計算する設計にはしない。
inactive cache にも有効期限を設け、新規対流の立上りを取り逃さないか検査する。

deep/shallow/inactive、top 到達、陸面の乾湿遷移、立上りの遅延、累積対流雨、CAPE 分布、
T/q の鉛直 profile と湿潤 enthalpy を比較する。
特に Phase 13 の発達後は対流降水が約97%を占めるため、全球総降水だけで判定しない。
900秒を先に採用候補へ登録し、1800秒は追加 gate を通した場合だけ使う。

## 7. P14.05 — BL/surface の時間刻み拡大と solve 改善

変更箇所: `boundary_layer.cpp/.hpp`、`dry_mixing_coupling.cpp/.hpp`、surface hydrology。

1. 現行の陰解法のまま dt=75/150/300/600/900秒の独立 column を比較する。
   強安定・強不安定、薄い最下層、弱風、低熱容量の陸、湿潤/枯渇 bucket、蒸発↔露を含める。
   K=20/40/80 と地表細分化を変え、係数固定誤差と solver 失敗を分ける。
2. 三重対角行列の factorization を同じ係数・dt の複数 RHS に再利用する。
   運動量成分、tracer 群、surface response で係数が一致するものだけを共有し、
   熱容量・熱/運動量拡散係数の違う行列は共有しない。
3. 水 flux の二分法を safeguarded Newton＋bisection と比較する。
   bucket の active set が変わる点では bracket を維持し、現行と同じ残差・水制約で受理する。
4. 600/900秒で誤差が支配的なら、Ri/K/bulk 係数の predictor 評価または限定 Picard 更新を試す。
   各反復は substep 始点から解き直し、同じ dt の水・熱を反復回数ぶん加えない。
   係数再計算の増分コストが substep 削減を上回る場合は採用しない。
5. analytic 更新は定係数の地表交換・drag の解析可能な部分から試す。
   大気と地表の相反する交換、運動エネルギー散逸の返還熱、水と潜熱を同一増分で収支化する。

600秒を標準候補、900秒を拡張候補とする。必要な柱だけ300秒以下へ局所細分化できるようにするが、
局所 fallback の割合・最小 dt・wall を報告する。頻繁な fallback を600/900秒の成功と数えない。

## 8. P14.06–P14.07 — Column physics と RHS 全体の OpenMP

### P14.06: column 内で処理を完結する

診断 → BL/surface → 乾燥調節 → SBM → 飽和調節 → 水の移送を、過程イベントに応じて
一つの水平 cell に対して実行できる kernel にする。列内では更新後の T/q/geometry を必要な範囲で再診断し、
古い `Derived` を後続過程へ渡さない。
全球 diagnose・state copy を何度も行う構造を、column workspace と column rollback へ置き換える。

thread ごとに geometry、Thomas solver、SBM profile、vapor、診断 scratch を事前確保する。
読み取り専用の state 入力とセル固有の出力を分け、全球の海洋水量などは column increment として返す。
セル単位の値が確定した後に固定順で集約して台帳へ一度だけ反映する。
未収束 column はその column 内を細分化し、最小 interval でも失敗したら外側の試行を全体棄却する。
OpenMP 領域を跨いで例外を投げず、失敗情報を回収後に通常の retry 経路へ戻す。

まず serial の column 化で旧 schedule の算術順を保ち、次に `schedule(static)` で水平 cell を並列化する。
SBM の負荷差が大きければ tile 単位の dynamic/guided も測る。診断集約は実行順に依存させない。
warm column kernel の heap allocation を0にする。

### P14.07: RHS と solver の残りを並列化する

対象は state diagnose、gradient/limiter、flux、CFL、vertical coupling、pressure gradient、
Coriolis、diffusion、fast operator、modal transform、Helmholtz 内の格子演算、残差計算、残る physics。

- **edge flux → cell gather** を基本形とする。edge に一意の flux を書き、各 cell が隣接 edge を
  固定順で集める。現行の左右 cell への scatter に単純な `omp parallel for` を付けない。
- source、diffusion、fast operator の共有 edge 更新にも同じ所有規約を適用する。
  atomics を hot loop の既定にせず、edge coloring はメモリ/速度比較の代案とする。
- max/min/count と浮動小数点の和を区別する。再現性モードは cell/tile 順の決定的集約、
  高速集約は登録許容差内一致を要求する。thread 間で driver の mutable workspace を共有しない。
- 外側の一つの parallel region と worker workspace を基本とし、reconstruction の既存 level 並列を
  内側で重ねない。4 thread の固定上限は、共通実行方針へ置換してから外す。
- GMRES の反復依存は保ち、独立した配列演算を並列化する。mode 並列 solve は mode 数が少ない場合や
  Krylov workspace が肥大化する場合に不利なので、cell 並列との比較後に選ぶ。

OpenMP のループ割当て・reduction の仕様は
[worksharing-loop](https://www.openmp.org/spec-html/5.2/openmpse66.html) と
[reduction clause](https://www.openmp.org/spec-html/5.2/openmpsu52.html) を参照する。
static schedule だけで浮動小数点 reduction の全 thread 数での byte 一致を保証するとはしない。
OpenMP OFF の build と既存 serial 経路を残し、CMake の option 説明・専用 preset を更新する。
ASan/UBSan、対応する compiler/runtime での race 検査、1/2/4 thread と繰り返し実行を検証する。

## 9. P14.08 — Reconstruction と flux の融合

変更箇所: `dry_hydrostatic_reconstruction.cpp/.hpp`、`dry_hydrostatic_flux.cpp`、driver/workspace。

1. gradient と limiter factor の準備を、face state 生成から分離する。
2. edge kernel 内で左右 state を再構成し、直ちに全 tracer を含む flux と CFL 寄与を計算する。
   左右セルへの加算で tracer flux を二度計算しない。
3. 全球の left/right face state 配列を production 経路から外し、まず小さい edge flux 配列を残す。
   cell gather で race-free に更新する。旧 face state API は fixture/debug 比較器として残す。
4. limiter、接ベクトル変換、共有 flux の反対称性、定数 tracer 保持、seam/corner と地形静止を検証する。
   bit-exact にできる部分と演算順の変わる部分を別コミットにする。

face state の除去だけで全 scratch が不要になるとはしない。gradient/factor、edge flux、gather の
bytes と時間を個別に記録し、N=48/K=40/4 tracer で peak scratch と full RHS wall を比較する。

## 10. P14.09 — 3-D 配列レイアウト

まず bounds/stride を表す `Field3DView`、`ColumnView`、tracer view を導入し、
現在の cell-column 配置のまま旧結果を保つ。checkpoint と可視化の論理順序は view の外側で固定する。

| 比較案 | 形と狙い | 計測する不利な点 |
|---|---|---|
| 現行 | `[cell][level]` | 水平 stride と level コピー |
| level-major | `[level][cell]` | BL/Thomas/SBM の column gather/scatter |
| blocked AoSoA | `[cell_block][level][lane]`、tracer/component は別配列 | padding、stencil gather、column の stride、変換コスト |

AoSoA の lane は4/8/16、block は複数 lane 群で比較する。Vec3 は block 内の x/y/z 別配列も試す。
K は SIMD lane に合わせて一律 padding せず、端数マスクと padding の両方を評価する。
最初は tile/workspace の pack に限定して実測し、利益があれば prognostic/derived の内部配置へ広げる。
kernel ごとに全球 transpose する二重レイアウトは避け、pack は一つの処理区間で使い回す。

水平だけでなく column physics・modal・診断・checkpoint の変換を含む model-day の wall と RSS で採否を決める。
1/2/4 tracer、K の端数、restart round-trip、名前による tracer 対応を検証する。
物理配列を変えてもファイル内の canonical な値の順序は維持し、scheduler 用 schema 変更と混同しない。

## 11. P14.10 — Modal transform / fast operator の SIMD・batch 化

変更箇所: `dry_hydrostatic_semi_implicit.cpp`、`dry_hydrostatic_fast_operator.cpp`、
fast modes/workspace、格子 divergence/gradient kernel。

cell block の鉛直量を `X[level, cell_lane]`、mode を `Z[mode, cell_lane]` とし、
既存の forward/inverse 係数規約を保って `Z=A*X`、`X=B*Z` として処理する。
momentum の3成分は同じ係数行列の複数 RHS として扱う。
行列の転置・pack は reference 更新時に行い、呼出しごとにやり直さない。

最初は mode/level の加算順を維持した cell-lane SIMD と loop blocking、次に小行列の
batched GEMM を比較する。K=20/40/80、selected mode 数1/3/5/全modeを測定し、
forward だけでなく inverse・pack・Helmholtz・残差まで含める。
selected modes だけに縮める場合も、非選択 mode の現行 correction を失わないよう代数的同値性を確認する。
BLAS は任意 backend とし、小行列の起動コストと OpenMP/BLAS の二重 threading を避ける。

固有分解そのものの方式変更は本項目に混ぜない。round-trip、線形性、保存 nullspace、
直接 fast operator との一致、modal equation residual、波の位相・減衰を検証する。
まず compiler report で SIMD を確認し、速度のためだけに全体へ `fast-math` を適用しない。

## 12. P14.11 — 高解像度の tile 処理

P14.08 の共有 edge 規約と P14.09 の view を利用する。各 cubed-sphere panel を
論理8×8/16×16/32×32の水平 tile に分け、実装された全 stencil の依存範囲から halo を決める。
gradient/limiter の入れ子の近傍参照、pressure source、二段 Laplacian を含む diffusion では、
一律「halo 1で十分」としない。

初版は halo/geometry を gather → gradient/limiter → reconstruction/flux → interior gather/source
を tile 内で処理する。必要な水平 flux が揃ってから全鉛直 column の coupling を行う。
tile 境界 edge は一意の owner が計算し、小さい境界 buffer と barrier を介して両側が参照する。
interior edge は tile scratch に置き、全球 face/flux buffer を段階的に縮める。
格子 seam/corner でも同じ edge ID・符号・接ベクトル規約を使う。

tile は in-place に予報 state を書き換えず、同一 stage の immutable input から tendency を作る。
隣接 tile が更新後の state を読んで時間離散が変わることを防ぐ。
全球 Helmholtz solve は局所 tile solve に置き換えず、その matvec の格子演算を tile 化する。

scratch bytes/thread は必要な fields×halo込みcell×K×tracer を含めて算出し、cache サイズに対して
余裕を取る。1/2/4 tracer と K=20/40/80 で tile 幅を選び、overlap 再計算と halo pack のコストも含める。
N=24/48/96 で wall・RSS・bytes/update・thread scaling を比較する。
tile サイズ、処理順、thread 数を変えて同じ保存量と登録許容差内の解を得ることを gate とする。

## 13. 共通の正しさ・性能 gate

### 数値・物理

新しい数値方式と、算術を変えない最適化を分ける。後者は同じ build の1 threadで
state/checkpoint の byte 一致を原則とする。演算順変更には field ごとの scale を決め、
独立 kernel の規格化誤差 `1e-12`、短い結合試験 `1e-10` を開始値として P14.01 で登録する。
単一の相対誤差をゼロ近傍の q/flux に適用せず、絶対 scale と相対 scale の和で判定する。
長期のカオス的瞬時場にはこの閾値を使わない。

時間方式・物理頻度変更の開始 gate は次のとおり。基準は旧3反復/300秒と、
Phase 13 の精細比較（450秒/5反復/物理75秒）。独立 column ではさらに dt を半減して収束を確認する。

| 段階 | 受け入れ案 |
|---|---|
| 独立 column | 水・enthalpy の交換量規格化残差 `<=1e-10`、正値性、飽和・bucket 制約、非降水浅対流の水保存。平滑な解析解では設計次数−0.3以上 |
| 1日、初期場＋活動中対流 | 質量/面積重み付き T RMS `<=0.1 K`、PW差 `<=1%`、時間積分E/P差 `<=5%`、TOA差 `<=0.5 W m-2` |
| 30日 | 全球平均T差 `<=0.2 K`、PW `<=2%`、時間積分E/P `<=5%`、TOA `<=0.5 W m-2`。帯状平均T RMS `<=1 K`、東西風 RMS `<=2 m/s` |
| 湿潤分布 | CAPEの平均/p95、雨量p95/p99の差 `<=20%`、対流降水比とtop到達面積の絶対差 `<=0.05`。雨量分布・鉛直加熱profileも併記 |
| 保存 | 乾燥・受動tracerの相対質量drift `<=1e-11`、1日区間の水・moist enthalpy の規格化残差 `<=1e-10`、q非負、pressure/層厚正、NaN/Infなし |
| retry | 既存 N=6 pilot 条件では原則0。ストレス試験は理由別に記録し、同じ精度の基準より短縮・retryが増えた効果も wall に含める |
| restart | 同じ方式/config/thread条件で連続runと再開runの全state・cache・台帳一致。異なるthread条件は登録数値許容差で比較 |

規格化残差の分母は既存 budget の定義を使い、ゼロ交換では P14.01 に登録する物理単位つき絶対 floor を使う。
E/P や tail 指標の基準値がほぼゼロのケースは相対差ではなく絶対差で判定する。
閾値と floor は候補結果を見る前に検証 manifest に固定し、変更理由を記録する。
全系エネルギー残差は過程の局所保存と別に比較し、長期平均の追加残差 `0.1 W m-2` 以下を開始目標にする。

30日瞬時3-D T RMS は報告するが、気候統計と同一の gate にしない。
最終候補は1200日を実行し、baseline と同じ区間の100日 block を比較する。
既存の `(200日,1200日]` と後半区間を別集計し、spin-up 境界で1日余分に含めない。
既存 pilot 集計器の境界誤差は P14.01 で修正し、raw CSV と突き合わせる。
平均差が block 間変動に比べて判定困難なら、共通の微小摂動3本・300日比較を追加する。
独立性のない日次値を独立サンプルとして信頼区間を過小評価しない。

### 性能目標

以下は同一環境で P14.01 に固定する目標。達成前に速度を確約しない。

| 対象 | 目標 |
|---|---|
| RHS回数 | 湿潤 ICI 2反復で通常5→4。第二候補は検査込み3以下を目指す |
| SBM | 活動中の状態で重い診断の3倍程度以上の削減（900秒候補）。再診断込みの実数で判定 |
| BL | 600秒候補で通常呼出し半減。局所fallbackと非線形反復を含めて評価 |
| serial全体 | N=6/12の湿潤で旧serial比1.5倍以上を第一目標 |
| OpenMP | N=48/K=40、4物理coreで新1thread比2.5倍以上。小格子ではserial閾値を設け5%以上の悪化を避ける |
| layout/fusion/tile | 対象のN=48以上でfull RHS wallを直前実装から20%以上削減、scratch削減量をbytesで報告 |
| 最終全体 | N=24/48の湿潤で旧serial比3倍以上を目標。N=6のspeedupとは別に報告 |
| allocation | warm RHS/column/modal kernelで0。driver・出力・checkpointのallocationは別集計 |

部分改善率を掛け合わせて総合速度を予言しない。`S=1/((1-f)+f/s)` を測定した時間比率に適用し、
改善後に支配箇所を再計測する。性能 gate は専用 runner で正しさの CI と分離する。
full RHS 回数が減っても dt 短縮や retry で同じモデル期間の wall が悪化した候補は採用しない。

## 14. 実装順と成果物

各行は一つのレビュー可能な milestone。必要に応じて算術不変の整理と数値変更を別コミットにする。

| ID | 内容 | 依存 | 完了時の成果物 |
|---|---|---|---|
| P14.01 | counters/timers、baseline、精度manifest、集計境界修正 | なし | benchmark JSON、再現コマンド、精度/性能閾値 |
| P14.02 | ICI2行列、predictor/IMEX比較と採否 | 01 | RHS/solve実回数、波/湿潤精度、採用方式ADR |
| P14.03 | legacy-compatible scheduler、過程間隔・放射分離、rollback/restart | 01 | scheduler API、結合順ADR、config/checkpoint契約 |
| P14.04 | SBM診断分離・900/1800秒・cache比較 | 03 | column/全球比較、診断回数と誤差 |
| P14.05 | BL600/900秒、factorization再利用、非線形/analytic比較 | 03 | 適用可能なdt、fallback分布、収支 |
| P14.06 | column局所化、workspace、column OpenMP | 04,05 | allocation=0、serial/parallel/retry一致 |
| P14.07 | flux gather、source/vertical/fast/solver全域OpenMP | 06 | ON/OFF CI、race検証、thread scaling |
| P14.08 | reconstruction–flux融合 | 07 | face配列削除、保存性、wall/bytes比較 |
| P14.09 | view導入、blocked/AoSoA比較 | 08 | 配列ADR、水平＋column総合評価 |
| P14.10 | modal/fast batch・SIMD | 09 | kernel/solver残差とend-to-end性能 |
| P14.11 | tile・halo・境界edge処理 | 08,09,10 | tileサイズ感度、N=24/48/96の比較 |
| P14.12 | 候補統合、30日・1200日、設定登録 | 採用する02–11 | Phase 14検証報告、results JSON、再現preset |

P14.02 と03は設計上独立して比較できるが、作業は上表の順に進め、方式変更の効果を混ぜない。
SBM cache/analytic、BL再線形化、IMEX、BLAS、全面AoSoA は比較して不採用にできる。
不採用時も根拠を残し、ユーザーが挙げた10項目すべてについて採否と実測を報告する。
RHS削減、物理の頻度分離、OpenMPの基盤、fusion、view/tileの評価は必須成果物とする。

新規の成果物は `docs/validation/phase-14.md`、`phase-14-results.json`、
`configs/phase14_*.cfg`、benchmark/比較 runner、方式・scheduler・配列/並列所有規約の ADR。
巨大な snapshot は output に置き、Git には config・生成手順・小さい集計を保存する。
最終報告には「採用/不採用/適用範囲限定」を明記し、数値正しさ・速度・既知の物理的限界を分けて記録する。

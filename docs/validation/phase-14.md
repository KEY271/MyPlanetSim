# Phase 14 検証記録

状態: P14.01 を実装・検証中。P14.02 の固定反復比較を実装。
predictor/IMEX の独立 ARK2 比較器とモデル結合まで実装。ARK2 は自然な構成の
短時間 screening で不採用。P14.03 以後は未着手。
速度の数値目標を満たすための diffusion、残差許容値、物理パラメータの調整は行わない。

## P14.01: 計測基盤

- RHS の initial / iteration / final を実呼出しで計数。discarded は破棄試行の部分集合であり総数に再加算しない。retry 前に共有する initial は一回だけ数える。線形 solve で先に失敗する試行は RHS を増やさない。
- 任意の RHS profiling を追加。診断/setup、reconstruction、flux/CFL、column coupling、source、diffusion、physics、result copy の排他的8領域。reconstruction 内の pack/limiter、線形 solve 内の modal/Helmholtz はまだ個別に分離していない。
- observer は full RHS を呼ばないので diagnostic full RHS は0。observer の物理診断は積分 RHS の回数に含めない。
- `benchmark_phase14.py` は同じコマンドの warm-up 後に通常5回、profile5回、allocation別1回を実行し、中央値/範囲、環境、binary/config SHA256 を保存する。
- moist benchmark の allocation 計数は `MPS_COUNT_ALLOCATIONS=1` 指定時だけ有効。通常runの0をallocationゼロの根拠にしてはならない。
- [精度 manifest](phase-14-manifest.json) を数値候補の試験前に登録。

```sh
cmake --build --preset release -j 4
ctest --preset release --output-on-failure -j 4
python3 tools/benchmark_phase14.py --output output/phase14/baseline.json -- \
  build/release/tools/benchmark_moist configs/phase13_moist_pilot.cfg \
  1 6 20 1800 3 20000 5 output/phase14/baseline.csv
```

Release の105テストが成功。計測有効/無効を含む継続積分のstate一致を確認。
Phase 13 の既存1200日CSVを修正版で再集計し、(200日,1200日] は
86,400,000秒、100日ブロック10個（最初200–300、最後1100–1200）となった。
積算fluxの区間が境界をまたぐ場合は推定配分せずエラーにする。

P14.01 の残作業: 地形等を含む比較行列、詳細physics/copy/solver計測、
高解像度・thread条件別baseline。初日の測定だけでSBM性能や方式の採否は判定しない。

### 発達後の比較 fixture

`benchmark_moist` の末尾に `INITIAL_CHECKPOINT SOURCE_CONFIG` を指定すると、
元configのfingerprintで検証した湿潤checkpointから指定日数だけ積分する。
物理設定・tracer registry・格子を維持し、時間積分と診断の設定だけ変更できる。
CLIのtop/refinementは元の座標と整合することを確認した上で、丸めを含めた元の
鉛直座標をそのまま用いる。これは比較専用のimportであり、production restartの契約は変えない。
蒸発・降水はimport時点の累積値を差し引き、水台帳も開始時の海洋・外部流出を含む。
共通runnerはcheckpointのSHA256も記録する。

```sh
python3 tools/benchmark_phase14.py --output output/phase14/developed-baseline.json -- \
  build/release/tools/benchmark_moist configs/phase13_moist_pilot.cfg \
  1 6 20 1800 3 20000 5 output/phase14/developed.csv \
  output/phase13_moist_pilot/pilot.chk configs/phase13_moist_pilot.cfg
```

1200日checkpointから1800秒進めるsmokeでRHS=5/step、retry=0、対流降水比0.9845。
`phase14.fixture` は実際に生成したcheckpointを使い、経過秒数・区間蒸発・RHS=4の
2反復import・同条件再実行のsnapshot byte一致・物理/鉛直座標不一致の拒否を検証する。
既存Release 105テストと新しいfixtureテストが成功。

## P14.02: 固定2/3/5反復の最初の比較

[比較契約 ADR 0020](../adr/0020-phase14-fixed-iteration-comparison.md) と
[実測JSON](phase-14-ici-results.json) を登録。数値コード・preset・線形許容値は変更していない。

```sh
python3 tools/compare_phase14_ici.py --output output/phase14/ici-initial
python3 tools/compare_phase14_ici.py --output output/phase14/ici-developed \
  --checkpoint output/phase13_moist_pilot/pilot.chk \
  --source-config configs/phase13_moist_pilot.cfg
```

両ケースでN=6/K=20、450/900/1800秒×2/3/5反復、各1モデル日を測定した。
各点はwarm-up＋通常5回＋profile5回＋allocation別1回。
同一binary、Release、OpenMP OFF。電源・CPU affinityは固定していない。
JSONに中央値と範囲、binary/config/checkpoint SHA256、build設定を記録した。
18条件すべてでretry/physics retry=0、RHSは2/3/5反復に対して4/5/7回/step。
通常/profile/allocationで比較対象の物理診断とRHS回数が完全一致した。
Release buildと107テストが成功。比較器の回帰テストは参照質量による重み付け、
格子・行数不一致、非有限温度の拒否を含む。

以下の温度RMSは450秒・5反復基準。速度は同一dtの3反復に対する実測である。

| 開始状態 | 力学dt [s] | 2反復 [s/day] | 3反復 [s/day] | speedup | 2反復 温度RMS [K] |
|---|---:|---:|---:|---:|---:|
| 初期場 | 450 | 4.5662 | 5.0693 | 1.1102 | 0.000000793 |
| 初期場 | 900 | 3.0011 | 3.2441 | 1.0810 | 0.0002637 |
| 初期場 | 1800 | 2.5323 | 2.6570 | 1.0493 | 0.0008910 |
| 1200日checkpoint | 450 | 4.7425 | 5.2265 | 1.1021 | 0.005092 |
| 1200日checkpoint | 900 | 3.0990 | 3.3697 | 1.0874 | 0.07370 |
| 1200日checkpoint | 1800 | 2.6280 | 2.7343 | 1.0404 | 0.09404 |

現行schedulerは最大300秒を等分するため、450秒の実physics intervalは225秒、
900/1800秒では300秒になる。dtを変えた差にはこの差も含まれる。
同一dtの5反復に対する2反復の温度RMSは、発達後で順に
0.005092 / 0.004543 / 0.04122 K。反復の影響と時間刻みの影響を混同しない。
発達後1800秒・2反復のPW差は0.00923%、降水差0.2265%、TOA差0.1774 W/m2。
測定した1日gateは18条件すべて通過したが、0.09404 Kは登録値0.1 Kに近い。

**採否:** 2反復をこの湿潤1日比較の候補として残す。標準presetは3反復を維持する。
既存線形波テストでは2反復の反復誤差が時間離散誤差を上回るため、全ケースへの一括採用はしない。
1800秒での全体短縮は初期場4.70%、発達後3.89%に留まり、RHSの20%削減を
全体の20%削減として報告しない。目標達成のための追加調整はしていない。

P14.02 の残作業はpredictor/IMEX比較、terrain/baroclinic等の適用範囲、75秒physics基準、
CAPE/降水tail・規格化enthalpy残差・長期統計。P14.01の詳細計測・高解像度baselineと、
P14.03以後の実装も未完了。今回の測定だけでこれらを完了扱いにしない。

### ARK2 IMEX 独立比較器

[ADR 0021](../adr/0021-phase14-ark2-imex-comparator.md) に SUNDIALS ARKODE
v7.7.0 の既定二次 additive pair の完全な表、stage 時刻、評価回数を固定した。
汎用 `Ark2Imex` の split linear ODE で二次収束を確認し、明示 RHS 3回、線形作用3回、
非自明な陰 stage solve 2回を回帰試験にした。これは production driver や preset の採用ではない。
モデル結合は3 stage の全 nonlinear RHS、全鉛直 mode の2回の陰 solve、3回の fast
operator、stage weight による収支を実装した。既存の最終 state/CFL 検査と次 step への
RHS 再利用を維持したため full RHS は ICI2 と同じ4回/step である。乾燥波の保存、stage
回数、solve 回数を回帰試験にした。

正式な12点反復測定の前に、同一 Release binary・初期 N=6/K=20・1日・dt=1800秒で
1回の screening を行った。ICI2 は2.5033秒、ARK2 は2.6243秒で、ARK2 は4.8%遅かった。
ARK2 は linear iteration 126→329、fast operator 0→144回/日となり、full RHS は双方
192回/日だった。平均温度差は `1.13e-5 K`、PW差は `4.37e-6 kg m-2`、retryは双方0。
これは分散を評価した性能測定ではなく、長い比較行列は以下のコマンドで利用者が実行する。

```console
python3 tools/compare_phase14_ici.py --output output/phase14/ici-with-ark2
```

**採否:** ARK2 は production preset に採用しない。full RHS を安全に3回以下へ減らすには
最終検査・収支・再利用契約をさらに変更する必要があり、現状は ICI2 より仕事量が多い。
10%以上高速という事前 gate に対して逆方向であり、係数、許容値、mode数の fine tuning は
行わない。比較器と runner は将来の方式検証用に残す。

## P14.03: 過程別 scheduler

[ADR 0022](../adr/0022-phase14-physics-scheduler.md) に legacy と過程別 event union、
力学step末での端数処理、同時刻の順序、旧共通 substep との競合エラーを固定した。
`make_physics_events` の 1800秒（BL 600 / SBM 900 / radiation 1800）と450秒の
短縮を unit test で確認した。driver の legacy 経路は既存 substep 数と呼出し順を維持する。
過程別経路では BL 200秒、SBM 300秒、力学600秒の event union が4イベント、BL 3回に
なること、retryなし、水 inventory 保存を explicit driver で確認した。失敗イベントは
state、bucket、累積台帳、step診断をまとめて戻してから due process を二分する。

過程別経路では gray radiation を full RHS から除き、event の物理時刻で1回診断して
その実 interval の保存量・地表温度更新へ適用する。600秒の explicit/semi-implicit
結合試験で radiation column call が従来の stage 回数ではなく全球1回になること、放射収支が
有限、水保存、retry 0を確認した。各力学step末で全 increment を適用し終えるため cache/時計を
checkpointへ持ち越さず、旧 moist layout のまま legacy/過程別 × explicit/semi-implicit の
連続1200秒と600秒再開が byte 一致する。observer の瞬時放射診断は積分状態を変更しない。

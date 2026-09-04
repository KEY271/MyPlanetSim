# Phase 9 実装計画 — 正しさの回復、性能、平衡ベンチマーク

**状態:** 計画。Phase 8 完了時点（`c3246a9`）で実施したコード全体レビューが見つけた欠陥と、
Phase 5--7 検証報告および ADR 0009 に記録済みの production gap を、新しい物理を追加せずに閉じる。
[docs/plan.md](plan.md) の当初 Phase 9 候補（灰色放射、湿潤過程、非静水圧、OpenMP/MPI）は
Phase 10 以後へ送る。

## 1. 到達点と設計原則

Phase 8 までで、cubed-sphere 幾何、球面輸送、shallow-water、hybrid 鉛直座標、乾燥 3D 静水圧コア、
固定地形、理想化物理、惑星・軌道・陸海 surface の契約と CI gate は揃った。しかし
**(a) 誤った演算子が一つ混じっており、(b) 設定が黙って無視される経路があり、(c) 時間刻み制限が
shallow-water と乾燥コアで異なる意味を持ち、(d) 「定常」ベンチマークが定常でなく、(e) 性能が
production matrix の実行を不可能にしている**。この五つは全て「基盤」の性質であって物理の不足では
ないため、次の物理を積む前に閉じる。

```text
Phase 8 完了 + 全体レビュー(2026-09-05)
        |
        v
Milestone 0: 契約の確定
  ADR 0011 / 0012 / 0013 と ADR 0009 amendment
        |
        v
Milestone 1: 出力を変えない正しさの修正
  vector Laplacian -> dry horizontal diffusion の配線
        |
        v
Milestone 2: 出力を変える正しさの修正
  cell 単位 CFL -> 壊れた DCMIP coordinate -> 平衡初期場と惑星依存 terrain
        |                                  (ここで baseline を一度だけ再登録)
        v
Milestone 3: 性能（結果を変えないことを gate にする）
  静的 geometry cache -> LS stencil -> workspace -> observer 間隔化 -> cost envelope
        |
        v
Milestone 4: 検証強化と衛生 -> Phase 5/6/7 production gate の unblock
```

設計原則を四つ置く。

1. **出力を変える修正を先に完了させ、baseline の再登録は一度だけ行う。** 性能作業を先に行うと、
   性能変更と物理変更の両方で fingerprint が動き、どちらが結果を変えたか切り分けられなくなる。
2. **性能作業は「同じ式を一度だけ評価する」memoization と「式を変形する」refactor を分離する。**
   前者は bit-exact を gate にできる。後者は登録許容差を明示する。両者を同じコミットに混ぜない。
3. **各修正は、その修正がなければ落ちる test を同じコミットで持つ。** 今回の欠陥はいずれも
   「形だけを検査する test」を素通りしており、後追いの test 強化では同じ欠陥を再び通す。
4. **性能ゲートは正しさゲートと分離する**（[docs/plan.md](plan.md) 4.1）。cell-update/s と
   メモリ/cell を独立に登録し、正しさ test の合否に性能を混ぜない。

## 2. 解決する問題

### 2.1 欠陥（本レビューで新たに検出）

| ID | 箇所 | 内容 | 実測による確認 |
|---|---|---|---|
| D1 | `src/numerics/spherical_operators.cpp:270-284` | 球面ベクトル Laplacian の回転成分の符号が逆。曲率項も欠落 | 単位球剛体回転 `u=ẑ×r` で `<∇²u·u>/<u·u> = +1.996`（厳密解 `−1`） |
| D2 | `src/dynamics/dry_hydrostatic_driver.cpp:196-199` | 水平時間刻み制限が cell ではなく edge 単位。四辺形セルで約4倍緩い | `phase7_held_suarez.cfg` の登録 candidate は表示 CFL 0.45 に対し実効 cell CFL ≈ 1.02 |
| D3 | `src/dynamics/dry_hydrostatic_benchmarks.cpp:173-176` | UMJS14 の平衡温度場が未実装。定数等温場 + 剛体回転風 + 一様 `ps` | `Linf(du/dt)` が `N=8/16/32` で `2.662/2.656/2.653e-3 m/s²` と収束せず（230 m/s/day） |
| D4 | 同上 `:134-150` | `umjs14_steady` と `umjs14_baroclinic` の初期場が bit-identical。摂動が無い | `max|ΔM| = max|Δθ·M| = max|Δps| = 0` |
| D5 | `src/dynamics/surface_orography.cpp:240-241` | JW06 地形が `config.planet` を無視し `R`, `Ω` を地球値で固定 | 非地球惑星で `orography.kind=jw06` を選ぶと風と釣り合わない地形が黙って作られる |
| D6 | `src/vertical/hydrostatic_column.cpp:40-43`, `src/physics/surface_energy_balance.cpp:70-75` | 「residual」診断が代数的恒等式で常に 0。欠陥を検出できない | 定義式から自明 |
| D7 | `src/diagnostics/dry_hydrostatic_diagnostics.cpp:35-43` | 乾燥コアの全球保存量が naive summation。shallow-water 側は `compensated_sum` | 41 万項で相対誤差 ~1e-13、測りたい drift と同オーダー |
| D8 | `src/physics/surface_energy_balance.cpp` | 地表 reservoir の陽的積分に安定性制約が無い。正値性しか検査していない | 小さい `C_surface` で `T_s` が振動しても gate を通る |
| D9 | `web/apps/gateway/src/main.js:172,200-201` | 完了 run の event/stderr と `.runs/` が解放されない。`error`→`close` の順で `activeRun` を二重解除し得る | repository に 63 個の run ディレクトリが残存 |

### 2.2 記録済みだが未解決の gap（Phase 9 が閉じる）

| ID | 出典 | 内容 |
|---|---|---|
| G1 | ADR 0009 L117-121 | 乾燥コアが `dry_hydrostatic.diffusion_*` を parse するが適用しない。このため Phase 7 の 6-run production matrix が blocked |
| G2 | [validation/phase-6.md](validation/phase-6.md) L94-99 | `configs/phase6_dcmip_2_0_0.cfg` が `a_half[0]=20544.8` と一様 `b_half` の組み合わせで起動しない |
| G3 | [validation/phase-5.md](validation/phase-5.md), [validation/phase-6.md](validation/phase-6.md) | Phase 5 の 3D 定量収束、UMJS14 envelope、Phase 6 の DCMIP 2-0-0 / Williamson 5 / JW06 production gate |
| G4 | [validation/phase-7.md](validation/phase-7.md), [validation/phase-8.md](validation/phase-8.md) | 1200 日 climate matrix、slow/rapid/tidally locked の長時間積分 |

G3/G4 の実行そのものは Phase 5--7 の成果物である。Phase 9 の責務は**実行を阻む二つの blocker
（G1 の拡散欠落と 2.3 の実行コスト）を除去し、実行可能であることを 1 本の pilot で示すこと**に限る。

### 2.3 性能

`configs/phase7_held_suarez.cfg`（`N=12`、`K=20`、864 cell × 20 level = 17,280 cell-level）、
clang `-O2`、単一スレッド、Apple Silicon での実測。

```text
100 steps = 22.0 s  ->  0.220 s/step  ->  7.85e4 cell-level-update/s
```

`sample` による self time 内訳（6410 サンプル）。

| 割合 | 内訳 |
|---:|---|
| 43% | `atan2` 1113、`acos` 605、`normalize` 324、`tangent_basis` 255、`inverse_map` 158、`tan` 100、`sin` 64、`map_to_unit_sphere` 52 |
| 12% | `neighbor_across` 340、`cell_index` 200、`cell_id` 93、`cell_edges` 58、`cell` 49 |
| 9% | `reconstruct_dry_hydrostatic_face_states` 本体 |
| 9% | `least_squares_gradient` + `least_squares_vector_gradient` |
| 4% | `pow` |

**55% 超が時間不変な格子ジオメトリとトポロジの再導出**である。原因は構造的で、次の三つに集約される。

- `reconstruct_tangent_vector`（`spherical_operators.cpp:218-242`）が呼ばれるたびに
  `inverse_map`（`atan2`×2）、`tangent_basis`（`tan`×2 + `normalize`×3）、
  `logarithmic_displacement`（`acos`+`sin`）を再評価する。これらは `(cell, edge)` ペアごとに固定量。
- RHS 1回あたりの最小二乗勾配の全格子パスは、reconstruction が `4K` スカラー + `K` ベクトル、
  source が `2K` スカラーで、`K=20` なら **140 パス / RHS、420 パス / step**。各パスが全 cell で
  上記の超越関数を再評価する。
- `couple_dry_hydrostatic_columns`（`dry_hydrostatic_coupling.cpp:21-62`）が cell あたり
  約 30 個の `std::vector` を確保する。`DryHydrostaticReconstruction` は `N=48`/`K=30` で
  **RHS ごとに約 93 MB** を確保・解放する。

観測される帰結。D2 の修正で登録 candidate の `dt` は 600 s から約 264 s へ下がるため、
1200 日 1 run は現状 **約 24 時間**、ADR 0009 の 6-run matrix は **150--250 core-hour** になる。
これが G4 が実行されていない実質的な理由である。

### 2.4 欠陥ではないもの（意図的な契約として維持する）

次はレビューで「地球固定」に見えるが、ADR で明示的に決定済みであり Phase 9 では変更しない。
変更したい場合は新しい ADR が必要であり、既存 benchmark kernel の exactness を壊してはならない。

- Held--Suarez / planetary Newtonian forcing の `day = 86400 s`、`315/60/10/200 K`、`σ_b=0.7`、
  `40/4/1 day`（ADR 0009 「No forcing constant is configurable」、ADR 0010 §5.1）。これは
  published benchmark を再現するための固定 kernel である。惑星ごとに時定数を変えたい場合は、
  既存 kernel を触らずに**別の設定可能な forcing 経路**を後から追加する。
- `ClimateStatisticsAccumulator` の 200 日 spin-up 境界（ADR 0009 「The standard window is
  day 200 through day 1200」）。
- Rusanov flux の散逸速度に音速を使うこと（ADR 0003 で選んだ基準スキーム）。ただし
  2.5 のとおり Phase 9 は**その散逸量を測って報告する**。

### 2.5 測定して報告するが Phase 9 では変更しないもの

- Rusanov の散逸係数が `|u_n| + √(γRT) ≈ 355 m/s` であり移流固有値 `u_n ≈ 30 m/s` の約 10 倍である
  こと。G1 を閉じると陽的拡散が初めて使えるようになるため、**暗黙散逸と陽的拡散を同じ単位で比較した
  測定**を Phase 9 の検証報告に載せる。flux 関数自体の変更（HLL/HLLC、advective-only dissipation）は
  この測定を根拠とする後続 ADR の対象とする。
- 地形追随座標の pressure-gradient error。予備測定（`K=8` の代替 A/B、DCMIP 2-0-0 の静止大気、
  厳密解は恒久的に `u≡0`）では

  ```text
   N     Linf|du/dt| [m/s2]     L2 [m/s2]
   8       7.39e-4              1.12e-4     (地形が未解像: 1 cell / 半波長)
  16       1.67e-3              1.89e-4
  32       1.64e-3              3.18e-4
  ```

  となり、地形が解像された `N=16 -> 32` でも `Linf` は減らず `L2` は増える。これは
  [docs/plan.md](plan.md) 4.3 の「解像度を上げても誤差が減らない」に該当し、Phase 6 完了ゲートの
  「静止大気の偽風と pressure-gradient error が解像度とともに減少」を満たさない。Phase 9 は
  G2 を直して**正規の座標で本測定を行い gate として登録する**ところまでを担い、離散化の変更
  （鉛直座標の再構成、Simmons--Burridge 型 PGF、参照状態差分）は測定結果を根拠とする後続 ADR とする。

## 3. スコープ境界

### 3.1 含めるもの

- D1--D9 の修正と、それぞれを落とす falsifying test
- G1（乾燥コア水平拡散の配線）と G2（DCMIP 座標の修復）
- 静的格子ジオメトリの前計算、最小二乗ステンシルの前計算、RHS の workspace 化、observer の間隔化
- cell-update/s とメモリ/cell の測定 harness、および登録コスト envelope
- 平衡ベンチマーク初期場の契約（JW06 / UMJS14 / DCMIP 2-0-0）と定常性・収束 gate
- gateway の retention と lifecycle 修正、Earth surface ingest の線形時間化
- 既存 CI gate の byte-exact 維持と、影響を受ける preset の一度きりの再登録

### 3.2 含めないもの

- 新しい物理過程。灰色放射、境界層混合、乾燥対流調節、湿潤過程、雲、詳細放射は Phase 10 以後
- 非静水圧方程式、semi-implicit / IMEX / split-explicit 時間積分
- Riemann solver の変更、鉛直座標の再構成、PGF 離散化の変更（2.5 のとおり測定と ADR まで）
- OpenMP / MPI / accelerator。**理由**: 現状は実行時間の 55% 超が冗長計算であり、それを並列化しても
  無駄を並列化するだけである。Milestone 3 のデータレイアウト整理が並列化の前提になる。ただし
  Phase 9 の refactor は RHS を隠れたグローバル状態から自由にし、cell 単位の作業を独立に保つことで
  後続の並列化を妨げないこと**は**要件に含める
- benchmark forcing 定数の設定化（2.4）と、惑星ジェネリック forcing の新経路
- binary checkpoint。テキスト checkpoint は `N=48`/`K=30` で約 50 MB/回になるが、production matrix の
  checkpoint 頻度が実測で問題になった時点で別途扱う
- Web UI の機能追加、FrameV3、production ETOPO/GSHHG 派生物の作成
- Phase 5--7 の production matrix そのものの完走（Phase 9 は unblock と 1 本の pilot まで）

## 4. 正しさの契約

### 4.1 球面ベクトル Laplacian（D1）

接ベクトル場 `v` に対し、半径 `a` の球面上で

```text
D    = div(v)
zeta = k_hat . curl(v)

lap(v) = grad(D) + k_hat x grad(zeta) + v / a^2
```

とする。第2項の符号が現在と逆であり、第3項は現在存在しない。第3項は曲率補正であり、
回転成分・発散成分のどちらに対しても同じ `+v/a^2` になる。次数 `n` の球面調和から作った場では

```text
rotational  v = k_hat x grad(psi_n) :  lap(v) = (-n(n+1)+1)/a^2 * v
divergent   v = grad(chi_n)         :  lap(v) = (-n(n+1)+1)/a^2 * v
```

となり、`n=1` の剛体回転では `lap(v) = -v/a^2` である。

gate は次を要求する。

- 剛体回転 `u = Omega x r` に対して `lap(u)` が `-u/a^2` に登録許容差内で一致する
- `n=1,2` の rotational/divergent 解析場で観測収束次数が 1 を下回らない
- 定数ベクトル場では第1・2項が exact zero、残るのは曲率項だけである
- `kLaplacian` の運動量拡散が滑らかな回転場の運動エネルギーを**減らす**
- `kBiharmonic` は符号反転が二回打ち消されて現在も正しい結果を出しているため、修正後も同じ符号で
  減衰することを確認する（この打ち消しが欠陥を隠していた）

### 4.2 水平時間刻み制限の正規化（D2）

有限体積の Courant 条件は面ごとではなく cell ごとに

```text
dt * sum_over_faces( lambda_f * L_f ) / A_cell <= cfl
```

である。`spherical_transport.cpp` と `shallow_water_driver.cpp` は既にこの形であり、乾燥コアだけが
`min_f( cfl * A / (L_f * lambda_f) )` を使っている。乾燥コアを同じ形へ統一し、
`dry_hydrostatic.cfl` と `shallow_water.cfl` と `transport.cfl` が**同一の意味**を持つことを
ADR 0011 で固定する。

影響範囲は測定済みで、次のとおり限定される（`lambda = |u| + sqrt(gamma*Rd*T)`、`gamma = cp/(cp-Rd)`）。

| preset | `N`/`K` | 現行 `run.time_step_s` | 修正後の cell CFL 上限（概算） | 影響 |
|---|---|---:|---:|---|
| `phase5_*`（8 件） | 4 / 8 | 10 s | ~780 s | なし |
| `phase6_dcmip_2_0_0` | 8 / 15 | 60 s | — | G2 修復後に再評価 |
| `phase6_jw06_steady`, `phase6_jw06_baroclinic` | 16 / 8 | 300 s | ~188 s | **あり** |
| `phase6_linear_mountain_wave` | 16 / 8 | 30 s | ~204 s | なし |
| `phase7_held_suarez_short` | 2 / 2 | 10 s | ~1700 s | なし |
| `phase7_held_suarez`, `phase8_earth_like`, `phase8_rapid_rotator`, `phase8_slow_rotator`, `phase8_tidally_locked` | 12 / 20 | 600 s | ~264 s | **あり** |
| `phase8_earth_geography` | 4 / 4 | 60 s | ~840 s | なし |
| shallow-water / transport 全件 | — | — | — | なし（既に正しい） |

影響する 7 preset のうち production run が既に登録されているものは無い（ADR 0009 の candidate も
Phase 6 の JW06 production も未実行）。したがって**無効化される既存 baseline は存在せず、全ての CI
gate は byte-exact のまま通る**。修正コストが最小である今が実施の適期であり、production matrix を
走らせた後では 6 run 分の再実行を伴う。

同じコミットで、地表 reservoir の陽的積分にも時間刻み制約を課す（D8）。

```text
tau_surface(c) = C_surface[c] / (4 * epsilon * sigma_SB * T_s[c]^3)
dt <= surface.cfl * min_c tau_surface(c)
```

`surface.cfl` は既定 0.5 とし、正値性 retry ではなく事前制約で振動を防ぐ。

### 4.3 乾燥コアの水平拡散（G1）

ADR 0006 は「Existing horizontal Laplacian/biharmonic machinery may be applied level by level by a
benchmark preset」と定めているが実装されていない。level ごとに既存演算子を適用し、次を守る。

- `diffusion_kind = none` は**現在の RHS と byte-exact に一致**する。既存 preset は全て `none` である
- 拡散は `air_mass`、`surface_pressure`、地形を変えない。`M*u`、`M*theta`、`M*q` にのみ作用する
- `stable_diffusion_time_step` を乾燥コアの `dt` 制約に加える
- diffusion tendency を advection / pressure gradient / Coriolis / vertical transport とは別項として
  budget に出す（ADR 0006「attribution of dissipation rather than hiding it in a residual」）
- 4.1 の修正後の演算子を使う。**壊れた演算子を配線してはならない**ため、D1 の修正が先行条件である

これで ADR 0009 の 6-run matrix の blocker が外れる。

### 4.4 平衡ベンチマーク初期場（D3, D4, D5, G2）

「steady」を名乗る preset は、**初期 RHS の残差が解像度とともに減少すること**を満たさなければ
ならない。現状の測定は次のとおり。

```text
Linf |du/dt| [m/s2]         N=8        N=16       N=32     判定
isothermal rest (flat)      0.0e+00    0.0e+00    0.0e+00  exact
JW06 steady (jw06 terrain)  7.19e-04   4.14e-04   3.03e-04 収束が 1 次未満で頭打ち
UMJS14 steady (flat)        2.66e-03   2.66e-03   2.65e-03 収束しない = 初期場が非平衡
```

- **UMJS14**（D3, D4）: 現在の実装は `u = 35 cos(phi)` の剛体回転風、定数等温場、一様 `ps`、
  平坦地形であり、UMJS14 でも steady state でもない。`steady` と `baroclinic` は bit-identical で
  摂動が無い。ADR 0013 で二択のいずれかを決める。
  1. Ullrich--Melvin--Jablonowski--Staniforth 2014 の解析的平衡温度場・地表気圧場と、規定の
     baroclinic 摂動を実装する。Phase 5 の未完了 production gate（G3）を閉じる正規の道
  2. UMJS14 を名乗るのをやめ、`solid_body_unbalanced` として非平衡 smoke test に降格する
  既定は 1 とする。2 を選ぶ場合は Phase 5 検証報告と plan.md の該当ゲートを同時に取り下げる。
- **JW06**（D5）: 地形の解析式は正しいが `R = 6.371229e6`、`Omega = 7.29212e-5`、`u0 = 35`、
  `eta0 = 0.252` がハードコードされ、benchmark 側と二重定義されている。定数を一箇所に集約し、
  `R` と `Omega` は `config.planet` から取る。JW06 は地球定数に対して定義された benchmark なので、
  惑星定数が JW06 の登録値と異なる config で `orography.kind=jw06` を選んだ場合は**黙って
  非平衡な地形を作らず明示的に拒否する**。
- **DCMIP 2-0-0**（G2）: `a_half = [20544.8, 0, ..., 0]` と一様 `b_half` の組み合わせは
  `p_half[1] = ps/15 = 6667 Pa < p_half[0] = 20544.8 Pa` となり validation で落ちる。
  `uniform_sigma_coefficients(20544.8, 15)` と同じ減少ランプに直し、preset が実際に起動することを
  smoke gate に入れる。修復後、2.5 の pressure-gradient error 収束測定を正規の座標で行う。

平衡性の gate は「残差が小さい」ではなく「残差が解像度とともに減る」で定義する。絶対値のしきい値は
Phase 9 の pilot 測定後に ADR 0013 で登録し、後から緩めない。

### 4.5 診断の健全性（D6, D7）

- 恒等式になっている residual を撤去し、独立に計算した量との差に置き換える。`HydrostaticColumn`
  では half-level の Φ を full-level 側から再構成した値と比較する。surface budget では、時間積分した
  reservoir エネルギー変化と、区間積算した flux の差を出す
- 乾燥コアの全球保存量（`dry_mass`、`potential_temperature_mass`、`tracer_mass`、`total_energy`、
  `axial_angular_momentum`）を `compensated_sum` に統一する。丸め由来の見かけ drift が
  測定対象の drift を汚さないこと自体を test で示す

## 5. 性能の契約

### 5.1 前計算する静的量

`CubedSphereGrid` の構築時（または grid が所有する派生 cache）に次を保持する。全て時間不変である。

```text
per cell            TangentBasis          (inverse_map + tangent_basis を1回)
per cell            1 / area_m2
per (cell, edge)    logarithmic_displacement  Vec3
per (cell, edge)    (x, y) = displacement を cell 接基底へ射影した成分
per (cell, edge)    neighbour flat index, edge sign
per (cell, edge)    outward unit normal, edge tangent basis
per cell            sum_f( sign * L_f * n_f )   shallow-water の pressure geometry correction
per cell            least-squares 逆行列から導いた隣接重み (wx, wy) x 4
```

`neighbor_across` / `edge_sign_for_cell` / `cell_index` / `cell_id` / `cell_edges` の境界検査つき
経路は hot loop から外し、フラット配列で置き換える。既存の検査つき API は残す。

### 5.2 bit-exact 段階と tolerance 段階の分離

| コミット | 種別 | gate |
|---|---|---|
| P9.07 静的 geometry cache | memoization。同じ式を同じ順序で一度だけ評価して保存する | **全 preset の出力が bit-exact** |
| P9.08 最小二乗ステンシル | refactor。`a00/a01/a11` の解を隣接重みへ畳み込む | 登録許容差（相対 `1e-12` を初期値とし pilot で確定） |
| P9.09 workspace 化 | 割り当て除去。算術順序は不変 | **bit-exact** |
| P9.10 observer 間隔化 | 診断の呼び出し頻度のみ変更 | **state と checkpoint が bit-exact**、診断行は間隔分だけ減る |

P9.08 だけが数値を動かす。ここを独立したコミットにするのは、fingerprint が動く原因を一つに限定
するためである。

### 5.3 割り当て契約

RHS 評価は steady state で追加のヒープ割り当てを行わない。

- `DryHydrostaticDriver` が RHS 用の workspace（level slice、column slice、reconstruction buffer、
  tendency buffer）を所有し、RHS は span で受け取る
- `couple_dry_hydrostatic_columns` と `vertical_scalar_rhs` は呼び出し側 buffer に書く。
  `interface_coordinate` / `center_coordinate` は column あたり一度だけ計算し 5 個のスカラー間で共有する
- `AtmosphericHybridCoordinate::geometry` の結果を diagnose 時に一度だけ作り、physics kernel は
  それを受け取る（現在は `diagnose` と `held_suarez_tendency` / `planetary_newtonian_tendency` で
  cell あたり二重に計算している）
- 割り当て回数そのものを test で固定する（グローバル `operator new` を数える test fixture）

同じコミットで hot path のスカラー処理も整理する。`speed_factor` の 54 回二分探索を二次方程式の
閉形式へ、`Vec3::operator[]` を throw しない `noexcept` アクセサへ、`pow(x,4)` を積へ、
`asin` の直後の `sin`/`cos` を `z` と `sqrt(1-z^2)` へ置き換える。いずれも値は同じか許容差内である。

### 5.4 observer の間隔化

`DryHydrostaticDriver::advance` は `obs` が非 null なら受理ステップごとに `diagnose(state)` と
physics kernel をもう一度フル評価する。一方 `PhysicsDiagnosticsAccumulator` は
`step % diagnostics.interval_steps == 0` の行しか保存しない（既定 144）。

- observer の signature に「この呼び出しで derived と physics rate が必要か」を渡す
- 毎ステップ必要なのは stage 重み付き累積量（`thermal_energy_contribution_j`、
  `rayleigh_drag_energy_contribution_j`）だけであり、これは RHS 側で既に得られている
- `ClimateStatisticsAccumulator` は時間重み付き積分なので、間隔サンプリングでも定義は保たれる。
  ただし**サンプリング間隔が統計値を変えないこと**を、二つの間隔で同じ run を比較して確認する

### 5.5 測定 harness と登録目標

正しさゲートから分離した性能ゲートを追加する。

```text
tools/benchmark_dry_core        固定 preset を固定ステップ数だけ回す
出力: cell-level-update/s, s/step, peak RSS/cell, RHS あたりの割り当て回数
```

registered target（`N=12`/`K=20`、単一スレッド、clang `-O2`、`phase7_held_suarez` 相当）:

| 指標 | 現状 | Phase 9 目標 | stretch |
|---|---:|---:|---:|
| cell-level-update/s | 7.85e4 | **>= 8.0e5**（10x） | 1.6e6（20x） |
| s/step | 0.220 | <= 0.022 | 0.011 |
| RHS あたりの heap 割り当て | 数十万 | **0**（steady state） | 0 |

根拠は 2.3 の内訳である。43% の超越関数と 12% のトポロジ経路を除去し、割り当てと observer の
重複を外すと 6 割強が消える。残りの算術に対して 2 倍程度の改善を見込んで 10x を必達、20x を stretch とする。

コスト帰結（1200 日、D2 修正後の `dt ≈ 264 s`、392,700 step）:

```text
現状          1 run 約 24 h      ADR 0009 の 6-run matrix 150--250 core-hour
10x 達成後    1 run 約 2.4 h     6-run matrix 15--25 core-hour
```

6 run は互いに独立なので、プロセス並列だけで実時間 3 時間程度になる。**MPI/OpenMP なしで
production matrix が実行可能になる**ことが Phase 9 の性能目標の意味である。

回帰防止として、CI では小さい preset の cell-update/s を測り、登録値から 30% 以上劣化したら失敗させる。
絶対値は機械依存なので、しきい値は同一 job 内で測った基準ケースとの比で表す。

## 6. 検証の強化

### 6.1 各修正に付随する falsifying test

今回の欠陥は全て「形を見る test」を通っている。修正コミットは次を同時に持つ。

| 対象 | 現状の test | 追加する test |
|---|---|---|
| D1 | ゼロ場のみ（`test_shallow_water_operators.cpp:59`） | 剛体回転の解析解、調和次数依存、KE 減衰の符号 |
| D2 | `stable_dt < time_step_s` の定性判定のみ | 解析セルでの `sum_f(lambda*L)/A` 一致、三つのソルバの CFL 定義が同値 |
| D3/D4 | 「有限で温度が 180--330 K」（`test_jw06.cpp`） | steady 残差の解像度収束、steady と perturbed が異なること |
| D5 | なし | 非登録惑星定数で JW06 地形が拒否されること |
| D6 | 恒等式を検査 | 意図的に壊した Φ 積分・flux で residual が非ゼロになること |
| G1 | なし | `none` の byte-exact、拡散有効時の分散減衰と budget 項の分離 |
| G2 | なし | 全 preset の起動 smoke gate |

### 6.2 既存の弱い gate の補強

- `test_spherical_operators.cpp:35-53` は誤差の単調減少しか見ていない。**観測収束次数**を
  3 解像度から算出して登録値と比較する（[docs/plan.md](plan.md) 4.1 Convergence の要求）
- `phase6.dcmip_2_0_0` を新設し、静止大気の偽風 `Linf`/`L2` が解像度とともに減少することを gate にする
- `phase5.umjs14` を新設し、4.4 の決定に応じて平衡性か非平衡の明示のいずれかを固定する
- 割り当て回数と cell-update/s の性能 gate（5.5）

### 6.3 CI 予算

`ctest --preset dev` は現在 69 test / 77 s である。Phase 9 の追加分を **+40 s 以内**に収める。
収束 gate は `N=8/16/32` の初期 RHS 残差だけを見れば足り、長時間積分は不要である。
DCMIP の 6 日積分と 1200 日 matrix は CI に入れず、検証報告の再現コマンドとして記録する。

## 7. 周辺の衛生

- **gateway**（D9）: 完了 run を LRU で上限個数まで保持し、超過分の event/stderr を解放する。
  `.runs/` 配下の run ディレクトリを retention と同じ規則で削除する。`child.on("close")` が
  自分の run だけ `activeRun` を解除するよう `activeRunId` で判定する。既存の loopback / Bearer /
  path traversal 対策は変更しない
- **Earth surface ingest**: `parse_earth_source`（`surface_boundary.cpp:72-75`）が行ごとに
  `std::ranges::find` で線形探索するため `O(rows * unique_lon)` になる。1 度格子（64,800 行）では
  問題ないが、production ETOPO（0.25 度 = 約 100 万行）では 1e9 回の比較になる。ハッシュ集合に
  置き換えて線形時間にする。派生物の内容と fingerprint は変えない
- **checkpoint**: テキスト形式は維持する。`N=48`/`K=30` で約 50 MB/回になることと、binary 化の
  判断基準（production matrix で checkpoint I/O が壁時計の 5% を超えたら）を検証報告に記録する

## 8. 変更境界

新規に想定するファイル。最終名は ADR 0011--0013 で固定する。

```text
docs/adr/0011-vector-operators-and-time-step-normalization.md
docs/adr/0012-static-grid-cache-and-performance-gates.md
docs/adr/0013-balanced-benchmark-initial-states.md
docs/validation/phase-9.md
include/myplanetsim/grid/cubed_sphere_cache.hpp
src/grid/cubed_sphere_cache.cpp
include/myplanetsim/dynamics/dry_hydrostatic_workspace.hpp
src/dynamics/dry_hydrostatic_workspace.cpp
tools/benchmark_dry_core.cpp
tests/test_vector_operators.cpp
tests/test_time_step_normalization.cpp
tests/test_dry_hydrostatic_diffusion.cpp
tests/test_balanced_benchmarks.cpp
tests/test_rhs_allocation.cpp
tests/test_preset_smoke.cpp
```

既存の変更は次に留める。`spherical_operators`（4.1）、`dry_hydrostatic_driver`（CFL・observer・
workspace）、`dry_hydrostatic_coupling` / `vertical_transport`（buffer 受け取り）、
`dry_hydrostatic_reconstruction` / `shallow_water_reconstruction`（cache 利用と閉形式 limiter）、
`dry_hydrostatic_benchmarks` / `surface_orography`（平衡初期場と惑星定数）、
`surface_energy_balance`（時間刻み制約）、`dry_hydrostatic_diagnostics` / `hydrostatic_column`
（residual と compensated sum）、`configs/phase6_dcmip_2_0_0.cfg` と影響 preset の `run.time_step_s`、
`web/apps/gateway/src/main.js`、`surface_boundary.cpp` の parser、CMake / test 登録。

ADR 0009 は amendment を必要とする（production candidate の `dt` と、blocked だった diffusion 設定）。
ADR 0006 は文言どおりに実装が追いつくだけなので変更しない。README と
[docs/plan.md](plan.md) には Phase 9 へのリンクと、Phase 9 候補が Phase 10 へ移動したことを記す。

## 9. コミット列

```text
1.  P9.01 docs: fix the correction, performance, and benchmark contracts
      ADR 0011/0012/0013、ADR 0009 amendment、影響 preset 表、登録許容差と性能目標

2.  P9.02 numerics: correct the spherical vector Laplacian
      符号 + 曲率項、解析解 gate。既存 preset は全て diffusion=none なので出力不変
3.  P9.03 dynamics: apply the configured dry horizontal diffusion
      level ごとの適用、dt 制約、budget 項の分離、none の byte-exact。ADR 0009 の matrix を unblock

4.  P9.04 dynamics: normalize the horizontal and surface time-step limits per cell
      4 面合算 CFL、三ソルバの CFL 定義統一、surface reservoir の tau 制約、影響 7 preset の再登録
5.  P9.05 config: repair the DCMIP 2-0-0 hybrid coordinate
      減少 A ランプ、全 preset 起動 smoke gate
6.  P9.06 dynamics: bind analytic terrain and steady states to the configured planet
      JW06 定数の一元化と惑星整合検査、UMJS14 の平衡初期場（または明示的降格）と摂動

7.  P9.07 grid: precompute static cell, edge, and neighbour geometry        [bit-exact]
8.  P9.08 numerics: precompute the least-squares gradient stencils          [許容差]
9.  P9.09 dynamics: give the RHS a preallocated workspace                   [bit-exact]
10. P9.10 driver: sample observer diagnostics at the configured interval    [state bit-exact]
11. P9.11 perf: register the cell-update benchmark and cost envelope

12. P9.12 test: gate balanced steady states and terrain pressure-gradient convergence
      JW06/UMJS14 残差収束、DCMIP 静止大気の偽風収束、演算子の観測収束次数
13. P9.13 diagnostics: replace tautological residuals and compensate global sums
14. P9.14 housekeeping: bound gateway retention, run directories, and surface ingest cost
15. P9.15 validation: record the Phase 9 correctness, performance, and unblocked-gate matrix
      暗黙散逸と陽的拡散の比較測定、pilot 1 run のコスト、Phase 5/6/7 gate の実行可否判定
```

```text
entry -> P9.01 -> P9.02 -> P9.03 ---------------------+
                    |                                 |
                    +-> P9.04 -> P9.05 -> P9.06 ------+-> P9.07 -> P9.08 -> P9.09 -> P9.10 -> P9.11
                                            |                                                  |
                                            +----------------> P9.12 <-------------------------+
                                                                 |
                                        P9.13, P9.14 (独立) ----+-> P9.15
```

依存の理由。

- P9.02 は P9.03 の前提。壊れた演算子を配線してはならない
- P9.04 は登録 candidate の `dt` を変えるので ADR 0009 amendment（P9.01）が先行する
- P9.05 は P9.12 の DCMIP 収束測定の前提。preset が起動しなければ測れない
- P9.06 は P9.12 の JW06/UMJS14 定常性 gate の前提
- P9.07--P9.10 は P9.02--P9.06 で確定した結果に対して bit-exact / 許容差を gate にする。
  順序を逆にすると baseline の再登録が二度必要になる
- P9.11 のコスト測定が P9.15 の「production matrix を実行してよいか」の判断材料になる

## 10. 完了チェックリスト

### 正しさ

- [ ] 球面ベクトル Laplacian が剛体回転と球面調和の解析解に一致し、Laplacian 運動量拡散が KE を減らす。
- [ ] 乾燥・shallow-water・輸送の CFL が同一定義であり、`cfl` の値が三者で同じ意味を持つ。
- [ ] `dry_hydrostatic.diffusion_*` が実際に適用され、`none` は byte-exact に従来と一致する。
- [ ] 全 preset が起動し、`phase6_dcmip_2_0_0` の座標が単調である。
- [ ] `steady` を名乗る preset の初期 RHS 残差が解像度とともに減少するか、非平衡であることを名前と
      文書で明示している。
- [ ] `baroclinic` preset が対応する `steady` preset と異なる初期場を持つ。
- [ ] 解析地形が `config.planet` と整合し、非整合な組み合わせを黙って受理しない。
- [ ] residual 診断が恒等式でなく、全球保存量が compensated summation である。

### 性能

- [ ] cell-level-update/s が登録目標（10x）を満たし、benchmark が再現手順とともに記録されている。
- [ ] steady state の RHS 評価が追加のヒープ割り当てを行わない。
- [ ] 静的 geometry cache と workspace の導入が全 preset で bit-exact である。
- [ ] 数値を動かす唯一のコミット（P9.08）の許容差が登録され、超過理由が説明されている。
- [ ] observer のサンプリング間隔が state と気候統計を変えないことが二間隔比較で示されている。
- [ ] メモリ/cell が記録され、`N=48`/`K=30` の見積もりが現実的である。

### 検証と引き渡し

- [ ] 各修正が、その修正がなければ落ちる test を同じコミットに持つ。
- [ ] 演算子の**観測収束次数**が登録され、単調減少だけの gate が残っていない。
- [ ] 暗黙 Rusanov 散逸と陽的拡散が同じ単位で比較・報告されている。
- [ ] 地形追随 PGF error の収束測定が正規の座標で行われ、結果が gate または後続 ADR の根拠として
      記録されている。
- [ ] ADR 0009 の 6-run matrix が blocker なしで実行可能であり、1 run の pilot コストが実測されている。
- [ ] CI 追加時間が +40 s 以内、GCC/Clang、ASan/UBSan、`format-check`、全 Phase 0--8 regression が通る。
- [ ] gateway の retention と lifecycle が test で固定されている。

## 11. Phase 9 が既存の判断に与える影響

| 対象 | 影響 |
|---|---|
| ADR 0006 | 変更なし。「diffusion may be applied level by level」の文言に実装が追いつく |
| ADR 0009 | **amendment 必要**。production candidate の `dt` が 600 s から cell CFL 由来の値へ変わる。blocked だった 6-run matrix の前提が満たされる。forcing 定数と 200 日窓は変更しない |
| ADR 0003 | 変更なし。Rusanov 基準スキームは維持。暗黙散逸の測定値を後続 ADR の入力として記録する |
| ADR 0008 | 変更なし。ただし PGF 収束測定が accepted 化の判断材料になる |
| ADR 0010 | 変更なし。benchmark forcing 定数は 2.4 のとおり固定のまま |
| validation/phase-5,6,7 | Known gaps のうち G1/G2 と実行コストの blocker が解消される。matrix の実行結果は各 phase の報告に追記する |
| docs/plan.md | 当初 Phase 9 候補（放射、湿潤、非静水圧、並列化）を Phase 10 へ移す。5 節「直近の実装順」を更新する |

## 12. 参照

- I. M. Held and M. J. Suarez (1994),
  [A Proposal for the Intercomparison of the Dynamical Cores of Atmospheric General Circulation Models](https://www.gfdl.noaa.gov/bibliography/related_files/ih9401.pdf):
  ADR 0009 の 1200 日 matrix の出典。Phase 9 はその実行 blocker を外す。
- P. A. Ullrich, T. Melvin, C. Jablonowski, A. Staniforth (2014),
  [A proposed baroclinic wave test case for deep- and shallow-atmosphere dynamical cores](https://doi.org/10.1002/qj.2241):
  UMJS14 の平衡初期場と摂動の定義。4.4 の選択肢 1 の出典。
- C. Jablonowski and D. L. Williamson (2006),
  [A baroclinic instability test case for atmospheric model dynamical cores](https://doi.org/10.1256/qj.06.12):
  JW06 の平衡初期場、地表 geopotential、定常性の受け入れ条件。
- P. A. Ullrich et al. (2012),
  [DCMIP 2012 Test Case Document v1.7](https://public.websites.umich.edu/~cjablono/DCMIP-2012_TestCaseDocument_v1.7.pdf):
  test 2-0-0（地形上の静止大気）と、pressure-gradient error を解像度で評価する手順。
- A. J. Simmons and D. M. Burridge (1981),
  [An Energy and Angular-Momentum Conserving Vertical Finite-Difference Scheme and Hybrid Vertical Coordinates](https://doi.org/10.1175/1520-0493(1981)109%3C0758:AEAAMC%3E2.0.CO;2):
  4.5 の hydrostatic residual を恒等式でない形にするための参照定式化。
- C. Schär, D. Leuenberger, O. Fuhrer, D. Lüthi, C. Girard (2002),
  [A New Terrain-Following Vertical Coordinate Formulation for Atmospheric Prediction Models](https://doi.org/10.1175/1520-0493(2002)130%3C2459:ANTFVC%3E2.0.CO;2):
  2.5 で観測した地形追随座標の pressure-gradient error と、その後続 ADR の候補となる定式化。
- W. C. Skamarock (2004),
  [Evaluating Mesoscale NWP Models Using Kinetic Energy Spectra](https://doi.org/10.1175/MWR2830.1):
  暗黙数値散逸と陽的拡散を同じ単位で比較する方法。2.5 の測定の設計指針。

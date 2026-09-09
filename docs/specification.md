# MyPlanetSim 実装仕様書

この文書は、現在のリポジトリに実装されている振る舞いを説明する。将来計画や開発 Phase の履歴は扱わない。仕様とコードが食い違う場合は、`include/`、`src/`、`apps/`、`web/` および自動テストを正とする。

## 1. 目的と適用範囲

MyPlanetSim は、理想化した惑星大気を対象とする研究・検証用シミュレーターである。実装範囲は次のとおり。

- cubed-sphere 上の球面 tracer 輸送
- hybrid sigma-pressure 座標の単一鉛直 column
- 乾燥および希薄水蒸気を含む 3 次元静水圧大気
- 固定地形、惑星・軌道条件、陸海混合地表
- 理想化 Newtonian forcing、灰色放射、乾燥・湿潤対流、bulk 境界層
- checkpoint/restart、診断 CSV、期間平均 dataset
- 保存済み dataset を読むローカル Web viewer

これは研究開発段階の実験的ソフトウェアであり、予報業務や production climate model としての妥当性を保証しない。全球 shallow-water model、実行中の対話制御、ブラウザーからの計算開始・初期値編集・瞬時 frame 再生は実装対象外である。

## 2. 構成

| 場所 | 責務 |
|---|---|
| `apps/my_planet_sim.cpp` | CLI、実験モード選択、出力、停止・再開 |
| `include/myplanetsim/`, `src/` | C++ model core、数値計算、物理、診断、I/O |
| `configs/` | 実行可能な設定 fixture |
| `data/earth/` | 地形・陸海面入力の checked-in fixture |
| `tests/` | unit、回帰、CLI、checkpoint、dataset のテスト |
| `tools/` | benchmark、比較、集計、可視化補助 |
| `web/packages/protocol/` | VisualDatasetV1 の読込・検証 |
| `web/apps/ui/` | React/Vite 製の読み取り専用 viewer |
| `output/` | 実行時生成物。Git 管理対象外 |

`configs/`、`tests/`、`tools/` 内の `phase*` 名は既存 fixture の安定した識別子として残っているだけで、仕様の段階分けを意味しない。

## 3. 共通規約

内部計算は SI 単位を用い、`mps::Real` は `double` である。設定名と公開 C++ member は、次元がある場合に原則として `_m`、`_s`、`_pa`、`_k`、`_m_s` などの単位 suffix を持つ。

惑星固定 Cartesian 座標は右手系で、`+z` が北極、`+x` が経度 0 の赤道、`+y` が東経 90 度の赤道を表す。経度は東向き、緯度は北向きを正とする。正の自転角速度に対する Coriolis 加速度は `-2 Omega x u` である。

鉛直 index `k` は上端から地表へ向かって増える。`K` 層に対し interface は `0..K`、full level は `0..K-1` である。interface 圧力は次式で定義する。

```text
p_half[k] = A_half_pa[k] + B_half[k] * surface_pressure_pa
```

`A` と `B` はそれぞれ `K + 1` 要素を持ち、圧力は下向きに単調増加し、各層の `Delta p / g` は正でなければならない。最終 step は `run.end_time_s` を越えないよう短縮される。

## 4. 実行モード

`experiment.kind` が実行経路を決める。

| 値 | 内容 | 主な出力 |
|---|---|---|
| `ode` | 指数減衰の製造 ODE。Euler または SSP-RK3 | 標準出力、任意 checkpoint |
| `sphere_transport` | cubed-sphere 上の prescribed-flow tracer 輸送 | `tracer.csv`、標準出力、任意 checkpoint |
| `vertical_column` | hybrid 座標の独立鉛直 column | `column_profile.csv`、`column_diagnostics.csv`、任意 checkpoint |
| `dry_hydrostatic` | 3D 静水圧大気と任意の物理・地表・水循環 | 診断 CSV、viewer dataset、任意 checkpoint |

`experiment.kind` を省略した場合の構造体既定値は `ode` だが、設定 parser は各モードに必要な key 一式を検証する。

## 5. 水平格子と輸送

水平格子は球面を 6 panel に分けた cell-centered cubed-sphere である。各 panel は `N x N` cell を持ち、総 cell 数は `6 N^2`。flat 配列は panel-major で、panel 内は `i` が最速で変化する。panel seam をまたぐ vector は Cartesian 空間を介して接平面へ変換する。

scalar の水平輸送は有限体積 flux で計算する。再構成は piecewise constant または linear、limiter は `none` または Barth-Jespersen を選択できる。独立輸送モードは `upwind` / `linear` scheme と、constant、Gaussian hill、cosine bell、slotted cylinder の初期場を持つ。流れ場は solid-body、deformational、divergent の test case から選択する。

乾燥コアでは unique edge に対して flux を計算し、tracer、質量、運動量、温位の conservative tendency を組み立てる。任意の Laplacian / biharmonic 水平拡散、格子量 cache、再利用 workspace、tile 分割、OpenMP kernel を備える。

## 6. 鉛直 column

hybrid 座標から interface/full-level 圧力、Exner 関数、層質量を診断し、温位から静水圧 geopotential と温度を積分する。予報量は地表気圧、層温位質量、tracer 質量である。鉛直輸送は donor-cell または linear reconstruction、limiter は `none` または `minmod` を選べる。

単独 column の test case は `isothermal`、`dry_adiabatic`、`moving_surface_pressure`、`manufactured_transport`。上端・地表 mass flux、質量・温位・tracer budget、CFL、静水圧残差、非有限値を診断する。

## 7. 3D 静水圧コア

### 7.1 状態

cell ごとの予報状態は以下からなる。

- 地表気圧
- 各 full level の水平 momentum mass（3 成分 Cartesian vector）
- 各 full level の potential-temperature mass
- tracer ごと・cell ごと・level ごとの tracer mass
- 地表物理を有効にした場合の地表温度
- bucket hydrology を有効にした場合の陸水量
- 湿潤実行時の累積対流性降水、格子スケール降水、蒸発、runoff、海洋水変化、外部流出

診断状態は layer 圧力・質量、接線風、温位、tracer mixing ratio、温度、geopotential を含む。各 step で地表気圧範囲、最小層厚、温度 floor、接線運動量、有限値、必要に応じ tracer 非負性を検証する。

### 7.2 初期・検証ケース

`dry_hydrostatic.test_case` は次を受理する。

```text
isothermal_rest
solid_body_transport
dcmip_deformational
dcmip_hadley
linear_wave
dcmip_2_0_0_rest
linear_mountain_wave
jw06_steady
jw06_baroclinic
umjs14_steady
umjs14_baroclinic
held_suarez
```

一部の case は対応する地形、Earth 定数、または物理設定を必須とする。整合しない組み合わせは config 読込時に拒否される。

### 7.3 時間積分と安定性

既定の乾燥コアは explicit SSP-RK3。`dry_hydrostatic.time_integrator` を指定すると、`semi_implicit` または比較用 `ark2_imex_comparison` を使える。semi-implicit 経路は鉛直 mode 分解で高速 Lamb・重力波成分を implicit に扱い、移流成分は explicit のまま残す。reference state は `fixed` または `per_step`、線形系は GMRES/Helmholtz solver を用いる。

要求 step は水平 fast-wave、移流、拡散、鉛直輸送、地表 reservoir、放射、物理過程の制約で短縮される。CFL・状態 invariant・線形 solver に失敗した step は、許容される最小 step まで縮小して retry する。semi-implicit の nonlinear iteration 数は固定回数で、残差は診断値である。

## 8. 地形・惑星・地表

固定地形は既定の flat surface のほか、`dcmip_2_0_0`、`linear_bell`、`jw06`、`latlon_csv` を実装する。flat は `orography.kind` を省略して選ぶ。CSV 入力は config からの相対 path、FNV-1a 64-bit fingerprint、smoothing 回数で同一性を管理する。

惑星定数は半径、自転角速度、重力加速度、乾燥気体定数、定圧比熱、基準圧力を設定する。軌道を使う場合は周期、離心率、赤道傾斜、近点経度、初期平均近点角、初期 substellar 経度、軌道長半径位置での恒星 flux を一組で指定する。

地表 geography は一様な陸率または Earth CSV fixture。Earth 入力は config 相対 path、fingerprint、cell 内 quadrature order、smoothing 回数を検証する。地表は陸・海の熱容量、初期温度、albedo、emissivity、内部熱 flux、粗度を持つ。checked-in Earth data は結合経路を検証する解析的 fixture であり、実観測から作った production dataset ではない。

## 9. 物理過程

### 9.1 大規模 forcing と放射

`physics.kind` は `none`、`held_suarez`、`planetary_newtonian`、`surface_energy_balance`、`gray_radiation`。

- Held-Suarez は Newtonian cooling と下層 Rayleigh drag を与える。
- planetary Newtonian forcing は `axisymmetric` または、軌道を使う `substellar` geometry を持つ。
- surface energy balance は陸海混合熱容量を持つ zero-depth surface と大気の熱交換を解く。
- gray radiation は短波吸収・散乱近似、地表反射、上下長波を同じ column flux から計算し、大気・地表 energy budget を記録する。

### 9.2 鉛直 mixing と対流

乾燥対流は `convection.kind = dry_adjustment`。bulk 境界層は `boundary_layer.kind = bulk_k_profile` と `backward_euler` integrator を用い、Richardson 数、turbulent Prandtl 数、gustiness、陸海別の momentum/heat roughness を設定する。どちらも gray radiation を前提とする。境界層を有効にすると legacy の地表顕熱交換を重複させないため `surface.air_exchange_coefficient_w_m2_k` は 0 でなければならない。

物理 scheduler は互換用 `legacy` と `process_intervals` を持つ。後者は受理された dynamics step 境界に合わせ、境界層・対流・放射の最大更新／診断間隔を短縮する。Simple Betts-Miller では `intermittent` または `cached_relaxation` update を選べる。

### 9.3 希薄水蒸気

複数 tracer は `tracers.names` で順序を定義し、各 tracer に `passive` / `water_vapor` role、初期 mixing ratio または初期相対湿度、非負制約、水平拡散の有無を指定する。湿潤物理は `moisture.kind = dilute_water` で有効になり、water-vapor tracer はちょうど 1 個必要である。

実装済み過程は saturation adjustment、Simple Betts-Miller 型対流、bulk 地表水蒸気交換、bucket 陸面水文である。bulk 水蒸気交換は bulk 境界層・地表・bucket を必要とする。潜熱、相変化、蒸発、降水、runoff、海洋 reservoir と外部流出を水・enthalpy budget に計上する。水蒸気の水平拡散を有効にする場合、乾燥コアの拡散は Laplacian でなければならない。

## 10. 設定ファイル

設定は UTF-8 の行指向 `key = value` 形式で、空行と `#` comment を許容する。未知 key、重複 key、不正な型、非有限値、必須 key の欠落、相互依存に反する組み合わせは error になる。list は comma 区切り、boolean は `true` / `false`。入力ファイル path は config ファイルの directory を基準に解決する。

主要 group は次のとおり。

| group | 内容 |
|---|---|
| `experiment` | 実行モード |
| `planet`, `orbit`, `star` | 惑星定数、軌道、恒星 flux |
| `run` | 開始・終了時刻、要求 step、乱数 seed |
| `grid` | panel あたり cell 数 |
| `transport`, `vertical`, `dry_hydrostatic`, `semi_implicit` | 数値 model と solver |
| `tracers`, `moisture` | tracer registry と水蒸気過程 |
| `orography`, `surface` | 地形、陸海面、熱・水 reservoir |
| `physics`, `radiation`, `convection`, `boundary_layer` | 物理過程と scheduler |
| `diagnostics`, `statistics` | sample 間隔と期間平均 |
| `output` | 生成先 directory |

全 key と既定値の型定義は `include/myplanetsim/config/experiment_config.hpp`、受理値と依存関係は `src/config/experiment_config.cpp`、実行例は `configs/` を参照する。config の canonical serialization から fingerprint を生成し、run metadata と restart 互換性検証に使う。

## 11. CLI と終了契約

```text
my_planet_sim --config PATH [options]

--integrator euler|ssprk3
--checkpoint PATH
--restart PATH
--stop-after-step N
--progress-interval-s SEC
--help
```

`--config` は必須。`--integrator` の Euler は ODE 専用で、他の mode は SSP-RK3 のみ受理する。乾燥コアの explicit/semi-implicit 選択は CLI ではなく config で行う。

乾燥コアは既定で 10 秒ごとに step、model day、進捗率、elapsed、ETA を標準 error へ出す。`--progress-interval-s 0` で無効化できる。`--stop-after-step` は absolute step 番号で停止する。SIGINT/SIGTERM は受理 step 境界で graceful cancellation し、指定済み checkpoint と partial diagnostics を書き、status `cancelled`、exit code 130 を返す。通常成功は 0、設定・実行 error は 1。

標準出力には model/build/config metadata、最終 status、時刻、step、主要診断を key-value 形式で出す。

## 12. 出力と restart

通常の生成先は `output.directory`。既存の同名 CSV と visual manifest は上書きされるため、run ごとに専用 directory を指定する。

乾燥コアの追加ファイルは有効機能に応じて生成される。

| 条件 | ファイル |
|---|---|
| Held-Suarez / planetary Newtonian / gray radiation | `physics_diagnostics.csv`, `climate_statistics.csv` |
| surface energy balance | `surface_state.csv`, `surface_diagnostics.csv` |
| gray radiation | `surface_state.csv`, `radiation_diagnostics.csv`, `radiation_column.csv` |
| semi-implicit | `semi_implicit_diagnostics.csv` |
| dry convection | `convection_diagnostics.csv`, `mixing_column.csv` |
| boundary layer | `boundary_layer_diagnostics.csv`, `mixing_column.csv` |
| multi-tracer / moisture | `tracer_diagnostics.csv`, `moisture_diagnostics.csv`, `moist_convection_diagnostics.csv`, `moist_column.csv`, `moist_surface_state.csv` |

一部の列指向 snapshot は該当過程が有効な場合だけ作られる。CSV schema は writer と CLI integration test が契約を保持する。

checkpoint は時刻、absolute step、flat state、config fingerprint、layout ID を保存する。restart 時は fingerprint、layout、state size を照合し、不一致を拒否する。地表・複数 tracer・湿潤 state は異なる layout ID を持つ。期間平均を有効にして checkpoint を書く場合は、同じ path に `.statistics` suffix を付けた sidecar を併記する。sidecar は checkpoint hash と時刻・step を照合してから accumulator と公開済み期間を復元する。

## 13. VisualDatasetV1 と viewer

乾燥コアは `output.directory/viewer/` に次を作る。

```text
viewer/
├── manifest.json
├── terrain.bin
└── means/
    └── period_000000.bin ...
```

`terrain.bin` と期間平均 binary は little-endian float64 payload と schema version 1 header を持つ。manifest は model kind、格子・鉛直座標、planet、config/terrain fingerprint、field offset、期間、byte length、coverage、完了状態を記録する。配列順は field-cell-level。

期間平均は `statistics.enabled = true` のときだけ生成される。受理 step の終端値を実時間幅で重み付けする rectangle rule を使い、`statistics.start_time_s` 以降を `statistics.period_s` の半開区間に分ける。途中終了区間は coverage と `complete = false` を持つ。

選択可能な標準 field は `surface_pressure`、`temperature`、`pressure`、`potential_temperature`、`eastward_wind`、`northward_wind`、`wind_speed`、利用可能なら `surface_temperature`、および `tracer.<name>`。`statistics.fields` が空なら利用可能な全 field、指定時は列挙した field だけを保存する。

viewer はユーザーが選んだ `viewer` directory を File System Access API で読み、期間・field・model level・cell の選択を 2D map、3D globe、鉛直 profile で共有する。viewer は simulator process を起動せず、dataset を変更しない。Canvas/WebGL の実描画は手動 smoke 対象で、protocol、数値変換、UI state は Vitest で検証する。

## 14. ビルド、並列化、検証

必要条件は CMake 3.25 以上、Ninja、C++20 compiler。preset は `dev`、`release`、`release-openmp`、`asan-ubsan` を提供する。OpenMP は `MPS_ENABLE_OPENMP=ON` の build で dry-core と column kernel に使われる。MPI と分散 I/O は実装していない。

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
cmake --build build/dev --target format-check
```

Web 側は npm workspace で protocol と UI を管理する。

```sh
cd web
npm ci
npm run lint
npm run typecheck
npm run test
npm run build
```

`just setup` は C++ build と npm install、`just dataset` は小格子 dataset の生成、`just dev` は localhost の viewer、`just check` は C++/Web のブラウザー不要 gate をまとめて実行する。

自動テストは geometry、輸送、鉛直座標、乾燥・湿潤力学、各物理過程、保存量、CFL、checkpoint/restart、config validation、CLI CSV、visual dataset、TypeScript protocol/UI を対象とする。benchmark executable と Python 比較 script は性能・候補比較用であり、通常の pass/fail test とは分離される。

## 15. 現在の制約

- hydrostatic primitive-equation model であり、non-hydrostatic 現象は扱わない。
- 水は dilute water-vapor 近似であり、雲微物理、雲 condensate の予報、氷相、詳細放射はない。
- 灰色放射、bulk 境界層、Simple Betts-Miller、bucket 地表は理想化 parameterization である。
- checked-in Earth surface は解析的 fixture で、production geography ではない。
- OpenMP shared-memory 並列のみで、MPI、GPU、並列 file format はない。
- viewer は保存済み期間平均専用で、live control や瞬時 state animation はない。
- 数値的な unit/regression gate を通していることと、現実の気候再現性は別である。

## 16. 変更時の仕様維持

公開挙動を変更する場合は、この文書、該当 config fixture、C++/Web test を同じ変更で更新する。特に config key、checkpoint layout、CSV header、VisualDatasetV1 schema を変える場合は、互換性を暗黙に破らず、version または明示的な移行処理を設ける。

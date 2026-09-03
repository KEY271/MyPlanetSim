# Phase 0 実装計画 — 科学・ソフトウェア基盤

**状態:** 2026-09-04 にローカル実装・検証完了。結果は
[Phase 0 validation report](validation/phase-0.md)を参照。

## 1. Phase 0 の到達点

Phase 0 では流体方程式をまだ解かない。Phase 1 以降で数値誤差と基盤バグを切り分けられるよう、次の土台を完成させる。

- C++20 の library / application / test ターゲットと、Debug・Release・sanitizer 構成
- CI 上の GCC / Clang ビルド、unit test、sanitizer test
- SI 単位、座標、添字、符号、浮動小数点比較に関する科学規約
- 検証付き惑星パラメータと、再現可能な実験設定
- halo を持つ最小 `Field2D` と periodic boundary
- Forward Euler と SSP-RK3 の明示時間積分
- 再現可能な reduction、誤差 norm、run metadata、乱数 seed
- version 付き checkpoint/restart
- 製造解 ODE による時間収束試験と Phase 0 検証報告

Phase 0 完了時の実行プログラムは、設定を読み、製造解 `dy/dt = -lambda*y` を積分し、解析解との誤差と metadata を出力できるものとする。2D 移流、shallow-water、cubed-sphere は Phase 1 以降で実装する。

## 2. 現在地

2026-09-04 時点の基準は次の2コミットである。

- `315794e docs: add staged development plan`
- `33fb523 build: add minimal C++20 executable`

後者をこの文書の P0.01 とみなす。CMake configure、Apple Clang 21 による C++20 build、空プログラム実行、CTest の smoke test までは完了済みである。

## 3. Phase 0 内で固定する判断

### 3.1 依存関係

Phase 0 の production code と test code は C++ 標準ライブラリだけで実装する。ネットワークからの build-time download は行わない。理由は、数値基盤の最初の参照実装を小さく保ち、offline build と再現性を確保するためである。

- test は CTest と小さな in-tree test support を使う。
- 設定は UTF-8 の strict `key = value` 形式とし、既知の flat key のみ受理する。
- checkpoint は `max_digits10` で round-trip 可能な version 付きテキスト形式から始める。
- NetCDF/HDF5、TOML/YAML/JSON library、MPI、OpenMP、logging library は Phase 0 に入れない。

構造化配列が必要になる hybrid level 設定や大規模 field I/O の前に、改めて外部ライブラリを ADR で選定する。

### 3.2 C++ 規約

- namespace は `mps`、公開 header は `include/myplanetsim/` 以下に置く。
- 浮動小数点型は `using Real = double` に集約する。
- 公開パラメータ名には `_m`, `_s`, `_kg`, `_pa`, `_rad_s` など SI 単位を含める。
- 所有配列には `std::vector`、非所有 view には `std::span` を使う。
- recoverable な入力エラーは説明付き例外とし、application 境界で捕捉して非0終了する。
- library から `std::cout` / `std::cerr` へ直接書かず、結果値または stream を呼出側から渡す。
- Debug では index と有限値を検査し、Release の検査省略は profile 後に判断する。

### 3.3 構成ファイル

Phase 0 の設定 key は最低限次を持つ。

```text
planet.radius_m
planet.rotation_rate_rad_s
planet.gravity_m_s2
planet.gas_constant_j_kg_k
planet.heat_capacity_cp_j_kg_k
planet.reference_pressure_pa
run.start_time_s
run.end_time_s
run.time_step_s
run.random_seed
ode.initial_value
ode.decay_rate_s_1
output.directory
```

空 key、未知 key、重複 key、非有限値、単位違いを暗黙に無視しない。設定 object から canonical な key 順で同じ形式へ書き戻せるようにする。

### 3.4 `Field2D` の最小責務

`Field2D<Real>` は Cartesian logical field の保存領域だけを担当する。

- interior size `(nx, ny)` と全周同幅 `halo` を constructor で確定
- row-major の contiguous storage
- logical index は interior を `0 <= i < nx`, `0 <= j < ny` とし、halo は負添字を許す
- `at(i,j)` は常に範囲検査、`operator()(i,j)` は Debug 検査
- interior view、全保存領域 view、`fill`、shape query を提供
- copy は deep copy、move は所有権移動

座標値、セル面積、物理境界条件、flux は持たせない。これにより Phase 1 の grid geometry と storage を分離する。

## 4. コミット運用規約

各コミットは次の条件を満たす。

1. 一つのレビュー可能な関心事だけを扱う。
2. 新しい挙動とその unit test を同じコミットへ入れる。
3. `cmake --build` と `ctest` が成功する状態だけをコミットする。
4. API rename と機能追加を混ぜない。必要なら rename を先行コミットにする。
5. generated file、`build/`、一時出力をコミットしない。
6. commit message は Conventional Commits 型の `type: imperative summary` を使う。
7. 後続コミットがなくても、そのコミット単体の目的と検証方法を履歴から理解できるようにする。

以下の message は予定であり、実装中に名前が変わった場合は意味を保って調整する。複数項目を一つに squash するのは、片方だけでは build 不能な場合に限る。

## 5. コミット列

### P0.01 — 最小 C++20 実行環境（完了）

予定/実績 message:

```text
build: add minimal C++20 executable
```

変更:

- `CMakeLists.txt`, `src/main.cpp`, `.gitignore`
- C++20 executable `my_planet_sim`
- `smoke.empty_program` CTest
- README の configure/build/run/test 手順

検証:

- clean configure と build
- executable が出力なし、終了コード0
- CTest 1件合格

### P0.02 — library / app / test のターゲット分離

message:

```text
build: establish project target structure
```

変更:

- `myplanetsim_core` library と `my_planet_sim` application を分ける。
- `include/myplanetsim/`, `src/core/`, `apps/`, `tests/` を実際に使う時点で作る。
- executable の入口を `apps/my_planet_sim.cpp` へ移す。
- install はまだ行わないが、公開/非公開 include path を区別する。

テストと受け入れ条件:

- core library を link した application が終了コード0。
- tests 有効/無効の両 configure が成功。
- source tree 内へ build artifact が生成されない。

### P0.03 — 依存なし unit test support

message:

```text
test: add dependency-free unit test support
```

変更:

- `tests/support/test.hpp` に test registration、`CHECK`、`CHECK_EQ`、近似比較、例外確認を実装する。
- 失敗時に file、line、式、actual/expected を表示する。
- test binary は失敗数を終了コードへ反映する。

テストと受け入れ条件:

- 成功する self-test を CTest に登録。
- 意図的失敗用 binary を `WILL_FAIL` で登録し、非0終了を確認。
- 1 test 内の複数失敗をまとめて報告できる。

### P0.04 — warning と format 規約

message:

```text
build: centralize compiler warnings and formatting
```

変更:

- target 単位の warning helper を `cmake/` に分離する。
- GCC/Clang は `-Wall -Wextra -Wpedantic`、MSVC は `/W4 /permissive-` を維持する。
- CI だけ warning-as-error にできる CMake option を追加する。
- `.clang-format` と `format-check` target を追加する。利用 tool の version は文書化する。

テストと受け入れ条件:

- library、application、test のすべてへ warning が適用される。
- warning-as-error 構成が成功。
- repository 対象ファイルの format check が成功。

### P0.05 — configure preset と sanitizer

message:

```text
build: add development and sanitizer presets
```

変更:

- `CMakePresets.json` に `dev`, `release`, `asan-ubsan` を定義する。
- sanitizer flag は専用 interface target 経由で core/app/test に適用する。
- compiler 非対応時には configure 時に明示エラーにする。

テストと受け入れ条件:

- `dev` と `release` preset が clean build できる。
- `asan-ubsan` で全 CTest が leak/UB report なしに合格。
- Release build に sanitizer flag が誤って入らない。

### P0.06 — GCC / Clang CI

message:

```text
ci: verify GCC Clang and sanitizer builds
```

変更:

- CI workflow に GCC Debug、Clang Debug、Clang ASan/UBSan、Release を登録する。
- configure/build/test/format-check を独立した step にする。
- compiler と formatter の version を log に残す。

テストと受け入れ条件:

- 全 matrix が clean checkout から成功。
- test failure と sanitizer failure が job failure になる。
- CI が `build/` や生成物の commit を必要としない。

### P0.07 — 科学・数値規約の文書化

message:

```text
docs: define scientific and numerical conventions
```

変更:

- `docs/conventions.md` に SI 単位、右手系、緯度経度、自転軸、Coriolis 符号を記載する。
- full/half level と添字順は将来仕様として図示する。
- missing value、NaN、許容誤差、保存量の relative drift の定義を固定する。
- 主要方程式の記号と C++ field 名の対応表を作る。

受け入れ条件:

- 同じ物理量に複数の単位・符号規約が残っていない。
- Phase 1 の変数名、座標方向、誤差 norm をこの文書だけで決められる。
- 未決事項は本文へ混ぜず、明示した open decision にする。

### P0.08 — 共通型と validation utility

message:

```text
core: add scalar types and validation utilities
```

変更:

- `Real`, `Index`, `Seed` と、有限値・正値・範囲を検査する utility を追加する。
- error message に parameter 名、actual value、必要条件を含める。
- `pi` など数学定数は `<numbers>` から一元参照する。

テストと受け入れ条件:

- NaN、Inf、ゼロ、負値、境界値の受理/拒否を unit test する。
- validation failure の message が対象 parameter を含む。
- `Real` が IEEE 754 binary64 という build-time 前提を検査する。

### P0.09 — 惑星パラメータ

message:

```text
core: add validated planet parameters
```

変更:

- `PlanetParameters` と Earth-like factory を追加する。
- 半径、自転角速度、重力、`Rd`, `cp`, 基準気圧を SI で保持する。
- `cv = cp - Rd`, `kappa = Rd/cp`, rotation period などの派生値を提供する。

テストと受け入れ条件:

- Earth-like 値と派生量が所定許容誤差内。
- `radius > 0`, `gravity > 0`, `cp > Rd > 0`, `reference_pressure > 0` を強制。
- 自転0は有効値として受理し、負の自転は逆回転として明示的に受理する。

### P0.10 — 決定的な実験設定

message:

```text
config: parse deterministic experiment settings
```

変更:

- strict `key = value` parser、`ExperimentConfig`、canonical writer を追加する。
- comment、空行、UTF-8、locale 非依存の数値 parsing を定義する。
- `configs/phase0_ode.cfg` を追加する。
- application に `--config PATH` と `--help` を追加する。

テストと受け入れ条件:

- sample config の全値が一致し、write/read round-trip が一致。
- 未知・重複・欠落 key、不正型、NaN/Inf、不可能な時間範囲を拒否。
- `LC_ALL` に依存せず小数点 `.` を解釈。
- `--help` は0、設定エラーは非0で終了し原因を標準エラーへ出す。

### P0.11 — checked `Field2D` storage

message:

```text
grid: add checked two-dimensional field storage
```

変更:

- contiguous row-major `Field2D<T>` と shape/halo query を追加する。
- interior/halo を含む logical index mapping と span access を追加する。
- overflow を考慮して allocation size を検証する。

テストと受け入れ条件:

- halo 0/1/2、非正方形 field の index mapping を全要素で確認。
- `fill`、deep copy、move、contiguous storage を確認。
- 0 extent、負 extent、範囲外、size overflow を拒否。
- sanitizer 下で境界アクセス test が clean。

### P0.12 — periodic halo boundary

message:

```text
grid: add periodic field halo filling
```

変更:

- scalar `Field2D` の x/y periodic halo fill を追加する。
- corner halo の fill 順序に依存しない写像を定義する。
- boundary API は後の cubed-sphere exchange と差し替え可能にする。

テストと受け入れ条件:

- 各 interior cell に固有番号を入れ、辺と四隅の期待値を全比較。
- halo width 1/2 と `nx != ny` を検証。
- halo fill が interior 値を変更しない。
- 定数 field は全保存領域で定数のまま。

### P0.13 — 明示時間積分器

message:

```text
numerics: add Forward Euler and SSPRK3 steppers
```

変更:

- `std::span<Real>` state と RHS callable を受ける Forward Euler / SSP-RK3 を追加する。
- stepper の scratch allocation は構築時または resize 時に行い、各 step では行わない。
- `t`, `dt`, state size の契約と alias rule を明記する。

テストと受け入れ条件:

- `dy/dt = 0` で bitwise に定数保持。
- 定数 RHS で1 step の期待値と一致。
- 複数成分が独立に更新される。
- `dt <= 0`、size mismatch、非有限 RHS を所定方針どおり拒否。
- test 用 allocation counter で steady-state step 中に heap allocation しない。

### P0.14 — reduction と誤差 norm

message:

```text
diagnostics: add reproducible reductions and error norms
```

変更:

- 固定順序 compensated sum、min/max、有限値検査を追加する。
- weighted `L1`, `L2`, `Linf` と relative drift を追加する。
- 空 field とゼロ reference norm の意味を明記する。

テストと受け入れ条件:

- 手計算できる scalar/vector 例と一致。
- 桁落ちしやすい列で naive sum より誤差が悪化しない。
- 同じ入力を反復して bitwise に同じ診断値。
- weight の負値、size mismatch、NaN/Inf を検出。

### P0.15 — run metadata と乱数 seed

message:

```text
io: record reproducible run metadata
```

変更:

- canonical config、Git revision/dirty flag、compiler、build type、開始条件を metadata に記録する。
- Git が利用不能な source archive build では `unknown` と明示する。
- `std::mt19937_64` は config の明示 seed だけで初期化し、wall clock seed を禁止する。

テストと受け入れ条件:

- 同じ config/seed の乱数列と metadata の決定的部分が一致。
- seed 変更で乱数列が変わる。
- metadata に全必須 key が一度だけ canonical order で出る。
- timestamp など非決定値は比較対象から分離される。

### P0.16 — version 付き checkpoint/restart

message:

```text
io: add versioned text checkpoints
```

変更:

- magic、format version、simulation time、step、state size、state、config identity を保存する。
- `max_digits10` と classic locale で `Real` を round-trip する。
- temporary file へ書いてから rename し、不完全 file を正式 checkpoint にしない。

テストと受け入れ条件:

- 複数の有限 `Real` が bitwise round-trip。
- truncated file、未知 version、不正 size、config mismatch を説明付きで拒否。
- temporary directory 内だけで I/O test を実行し、終了時に片付ける。
- 読込失敗時に既存 state を部分更新しない。

### P0.17 — 製造解 ODE application

message:

```text
app: run the Phase 0 manufactured ODE
```

変更:

- config から `y0`, `lambda`, `dt`, `t_end` を読み `dy/dt = -lambda*y` を積分する。
- `--integrator euler|ssprk3`, `--checkpoint`, `--restart` を接続する。
- 解析解、絶対/相対誤差、step 数、metadata を machine-readable な `key = value` で出力する。

テストと受け入れ条件:

- sample config の実行が0終了し、全出力値が有限。
- 最終 step は `t_end` を超えず、割り切れない場合の短い最終 `dt` を正しく扱う。
- clean run と中間 restart run の最終 `(time, step, state)` が一致。
- 不正 CLI と書込不能 output で非0終了。

### P0.18 — 時間収束と restart の自動ゲート

message:

```text
test: enforce Phase 0 convergence and restart gates
```

変更:

- `dt`, `dt/2`, `dt/4`, `dt/8` を実行する convergence test を追加する。
- 観測次数と誤差表を test failure 時に表示する。
- checkpoint なし/あり、同一 seed 反復の end-to-end test を追加する。

テストと受け入れ条件:

- Forward Euler の観測次数が `>= 0.7`。
- SSP-RK3 の観測次数が `>= 2.7`。
- `dt` 細分化ごとに誤差が単調減少。
- restart 前後の state は bitwise 一致。metadata の決定的部分も一致。
- ASan/UBSan を含む CI matrix で全 gate が合格。

### P0.19 — Phase 0 検証報告

message:

```text
docs: record Phase 0 verification results
```

変更:

- `docs/validation/phase-0.md` に compiler、preset、test 件数、収束表、sanitizer 結果を記録する。
- 実行コマンド、使用 config、commit hash を記録する。
- Phase 1 へ持ち越す制約と open decision を列挙する。
- README の現在状態を Phase 0 完了へ更新する。

受け入れ条件:

- clean checkout から報告値を再生成できる。
- P0.01〜P0.18 の全テストと CI が green。
- [全体計画](plan.md) の Phase 0 完了ゲートを満たす証拠が報告内に対応付けられている。

## 6. 依存関係と並行可能性

```text
P0.01
  └─ P0.02 ─ P0.03 ─ P0.04 ─ P0.05 ─ P0.06
       ├─ P0.07 ─ P0.08 ─ P0.09 ─ P0.10
       ├─ P0.11 ─ P0.12
       └─ P0.13 ─ P0.14

P0.10 + P0.13 + P0.14 ─ P0.15 ─ P0.16
P0.16 + P0.17 ─ P0.18 ─ P0.19
```

番号順の実装を基本とする。P0.07 文書は P0.08 以降の命名を拘束するため先行させる。`Field2D` 系と時間積分系は API 境界が固まれば並行に検討できるが、共有ファイルの競合とレビュー負荷を避けるため、履歴上は番号順に統合する。

## 7. 各コミット共通の検証コマンド

通常コミット前:

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev --output-on-failure
cmake --build --preset dev --target format-check
git diff --check
```

数値計算、メモリ、I/O に触れるコミットでは追加で:

```sh
cmake --preset asan-ubsan
cmake --build --preset asan-ubsan
ctest --preset asan-ubsan --output-on-failure
```

P0.05 より前は preset が存在しないため、README にある `cmake -S . -B build` 形式を使う。テスト結果を通すために許容誤差を緩める場合は、同じコミットに理由と数値根拠を記録する。

## 8. Phase 0 完了チェックリスト

- [x] clean checkout から GCC と Clang の Debug/Release build が成功する。
- [x] ASan/UBSan で全テストが成功する。
- [x] format-check と warning-as-error が成功する。
- [x] 単位・座標・符号・添字・誤差 norm の規約が文書化されている。
- [x] `PlanetParameters` と config が対象とする不正値を明示的に拒否する。
- [x] `Field2D` と periodic halo の全辺・全 corner がテストされている。
- [x] Forward Euler と SSP-RK3 が製造解で期待収束次数を示す。
- [x] 診断値、乱数列、canonical metadata が同じ入力で再現する。
- [x] checkpoint/restart 後の最終 state が連続実行と一致する。
- [x] Phase 0 検証報告を記載し、再生成コマンドを示している。
- [x] Phase 1 の 2D conservative advection に必要な基盤 API が揃っている。

## 9. Phase 0 から除外するもの

次は Phase 0 完了のために実装しない。

- Cartesian grid の座標・セル面積・flux divergence
- 2D advection と CFL controller
- shallow-water 方程式
- cubed-sphere geometry と panel exchange
- hybrid sigma の `A/B` 配列
- NetCDF などの大規模 field 出力
- MPI/OpenMP/GPU 最適化
- 実時間性能の合否判定

これらを先に追加すると Phase 0 の完了条件が動くため、必要性が見つかった場合も Phase 1 以降の issue または ADR として記録する。

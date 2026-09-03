# Phase 2 実装計画 — cubed-sphere 全球 shallow-water

**状態:** 実装済み。Phase 1 完了コミット `c3d349d` を開始点とし、結果は
[Phase 2 検証報告](validation/phase-2.md) と
[ADR 0003](adr/0003-shallow-water-horizontal-discretization.md) に記録する。
未達項目は検証報告の «Known gaps» に列挙する。

## 1. 到達点

Phase 2 では Phase 1 の全球 geometry、unique edge、Cartesian 接ベクトル、
共有 flux、診断、checkpoint を再利用し、回転球面上の1層 shallow-water
方程式を実装する。完了時には、基準となる cell-centered finite-volume 法と
compatible/mimetic 候補を同じ実験 API で比較し、標準ケース、保存量、
seam/corner 誤差、計算量に基づいて Phase 4 以降で採用する水平離散化を
ADR で固定できることを目標とする。

成果物は次を含む。

- 正の layer depth と Cartesian 接運動量を持つ version 付き全球 state
- 質量・運動量の unique-edge numerical flux と Coriolis/球面拘束 source
- piecewise-constant および limited linear reconstruction
- Rusanov 基準 flux と、選定した compatible/mimetic flux 候補
- gravity-wave CFL、SSP-RK3、depth positivity、failure diagnostics
- Laplacian/biharmonic diffusion と mass/energy/enstrophy budget
- Williamson test 2、Williamson test 6、Galewsky–Scott–Polvani jet
- mass、energy、potential enstrophy、axial angular momentum、PV、
  balance、seam/corner、grid-imprinting 診断
- shallow-water 用 CLI、設定、checkpoint/restart、CSV、Phase 2 検証報告

Phase 1 の transport 実装と設定は回帰対象として維持する。Phase 2 の作業中も
Phase 0/1 gate を無効化しない。

## 2. Phase 2 に含めないもの

- orography と topographic source term。Williamson test 5 は Phase 6 へ回す。
- multilayer、baroclinic、hydrostatic primitive equations（Phase 5）
- hybrid sigma-pressure、鉛直移流、surface pressure（Phase 4/5）
- moisture、放射、境界層、外力による気候積分
- MPI/OpenMP/GPU と大規模 parallel I/O
- semi-implicit、split-explicit、exponential time integration
- adaptive mesh refinement、高次 DG/SE の追加比較

明示 SSP-RK3 で gravity-wave CFL が律速する範囲を Phase 2 の性能基準とする。
時間積分方式の高度化は、水平離散化の誤差を分離した後に判断する。

## 3. 連続方程式と変数

### 3.1 基本方程式

地形なし、一定重力の回転 shallow-water 系を用いる。

```text
partial_t h + div(h u) = 0

partial_t u + (zeta + f) k cross u
  + grad(g h + 0.5 |u|^2) = D_u
```

ここで `h` は layer depth、`u` は球面接速度、`k` は外向き単位法線、
`zeta = k dot curl(u)`、`f = 2 Omega dot k` である。基準 finite-volume 法は
同値な質量・接運動量 flux form を使い、各 stage 後に Cartesian 運動量を
cell center の接平面へ射影する。球面曲率項を暗黙に落とさないよう、Cartesian
tensor flux の surface divergence と接平面射影の関係を P2.01 の ADR に式で
固定する。

予報量は cell average の

```text
h_c                 [m]
m_c = h_c u_c       [m^2 s^-1]
```

とする。`m_c` は3成分 Cartesian vector だが、独立自由度は接平面内の2成分で
あり、`dot(m_c, cell_center)=0` を invariant とする。速度は `u=m/h` として
診断する。`h <= depth_floor_m`、NaN/Inf、非接運動量を stage 境界で検出する。

### 3.2 保存量

密度を一定かつ共通因子として省略し、次を診断する。

```text
mass = integral h dA
energy = integral (0.5 h |u|^2 + 0.5 g h^2) dA
PV = (zeta + f) / h
potential_enstrophy = 0.5 integral h PV^2 dA
axial_angular_momentum = integral h [ (r cross u)_z
                              + Omega R^2 cos(latitude)^2 ] dA
```

Rusanov 法では mass のみを構造的保存量とし、energy/enstrophy の単調性や drift
を報告する。compatible 候補では mass に加え、離散 energy と potential
enstrophy のどちらを exact/near-exact に保存する設計かを API と ADR で明示する。
「保存的」という語だけで両方を主張しない。

## 4. 固定する数値設計

### 4.1 基準 finite-volume 法

- 配置は collocated cell-centered `h` と Cartesian tangent momentum `m`。
- edge center で left/right primitive state を共通局所基底
  `(edge_normal, edge_tangent)` へ変換する。
- 1D normal shallow-water 系へ local Lax–Friedrichs/Rusanov flux を適用し、
  spectral radius は `max(|u_n| + sqrt(g h))` とする。
- mass flux と Cartesian momentum flux は unique edge ごとに一度だけ計算し、
  left/right cell へ反対符号で scatter する。
- pressure flux は `0.5 g h^2 edge_normal`。tangential momentum は同じ mass
  flux で移流し、mass–momentum consistency を保つ。
- 一次方式を最初の reference とし、Phase 1 の tangent-plane least-squares を
  scalar/vector primitive reconstruction へ拡張する。
- depth limiter は face state が `depth_floor_m` 以上になるよう slope を縮小する。
  velocity slope は局所 neighbor bounds と kinetic-energy spike の両方を検査する。

Rusanov dissipation と明示 diffusion を分けて budget 化する。数値 flux の散逸を
明示 diffusion と誤認しない。

### 4.2 Coriolis と球面拘束

planet-fixed Cartesian frame の回転 vector を
`Omega=(0,0,planet.rotation_rate_rad_s)` とし、Coriolis acceleration は

```text
-2 project_tangent(Omega cross u, k)
```

とする。flux tendency、Coriolis、diffusion を合成した後、各 SSP-RK3 stage の
運動量を接平面へ射影する。射影前後の radial momentum を診断し、射影が大きな
誤差を隠していないことを gate にする。

静止一様 depth は bitwise または丸め誤差 scale で静止を保つ。Williamson test 2
では pressure gradient と Coriolis の個別 tendency、およびその残差を出力する。

### 4.3 compatible/mimetic 候補

基準法とは別に、primal cell の `h`、primal edge の法線 velocity/mass flux、
primal vertex 周りの dual cell で `zeta/PV` を持つ C-grid 候補を実装する。
Phase 1 の頂点・辺 incidence から次を構成する。

- signed cell-edge incidence `D2`
- signed edge-vertex incidence `D1`
- dual cell area、cell-center–edge-center segment、dual edge length
- `D2 D1 = 0` の exact topology identity
- edge-normal velocityから cell kinetic energy と vertex circulation への写像
- non-orthogonal cubed-sphere edge の metric/Hodge correction

PV flux operator は antisymmetry と constant-PV consistency を unit test で固定する。
energy/enstrophy conserving formを候補とするが、cube corner の3-valent vertex と
非直交 dual metric で条件を満たせない場合は、条件を黙って緩めず P2.18 の比較
結果に記録する。

P2.01 では基準 A-grid を採用する判断だけを固定する。compatible 候補の最終採否は
Williamson/Galewsky 比較後の P2.19 ADR で行い、結果を先取りしない。

### 4.4 diffusion

設定で `none`, `laplacian`, `biharmonic` を切り替える。

- Laplacian coefficient は `m^2 s^-1`、biharmonic coefficient は `m^4 s^-1`。
- vector diffusion は surface `grad(div(u)) - k cross grad(zeta)` を基本とし、
  Cartesian component ごとの平面 Laplacianだけで済ませない。
- depth diffusion を有効にする場合は全球積分が厳密に0となる edge flux formを使う。
- diffusion tendencyによる kinetic/potential energy、enstrophy の変化を別 budget
  として記録する。
- coefficient は解像度と e-folding time から導出できる helper を用意し、
  config へ根拠のない resolution-dependent magic number を置かない。

### 4.5 time step と positivity

各 cell について incident edge の

```text
edge_length * (|u_normal| + sqrt(g h))
```

を合計し、cell area と `shallow_water.cfl` から stage 共通の time step を求める。
設定上の `run.time_step_s` を上限とし、最終 step は終了時刻へ揃える。
diffusion 有効時は Laplacian/biharmonic operator の安定上限も同じ helper で評価し、
gravity-wave CFL と diffusion CFL の小さい方を採る。

Rusanov 一次方式が正値となる CFL を unit test で確認する。線形再構築は face depth
positivity limiter を必須とし、更新後の小さな負値を clip して質量を失う実装は
禁止する。`h < depth_floor_m` が発生した場合は cell、stage、CFL、最小値を含む
error で停止する。

## 5. 設定 schema と state layout

`experiment.kind = shallow_water` を追加する。最低限の schema は次とする。

```text
experiment.kind = shallow_water
grid.cells_per_panel
shallow_water.test_case
shallow_water.scheme
shallow_water.reconstruction
shallow_water.limiter
shallow_water.cfl
shallow_water.mean_depth_m
shallow_water.depth_floor_m
shallow_water.diffusion_kind
shallow_water.diffusion_coefficient
diagnostics.interval_steps
```

test case 固有値は preset 内で canonical value を明記し、論文既定値を code 内の
暗黙定数にしない。ただし初期条件を生成する名前付き factory は標準値を返せる。
未知 key、別 experiment 用 key、非有限値、`h<=floor`、`0<CFL<=1` 違反、係数の
単位/order 不整合を strict に拒否する。

checkpoint の state layout は version を持ち、panel-major cell order で

```text
depth[all cells], momentum_x[all cells],
momentum_y[all cells], momentum_z[all cells]
```

とする。compatible scheme が edge state を持つ場合は別 layout ID を使う。
既存 version 1 checkpoint の ODE/transport reader を壊さず、layout と expected
state size の不一致を restart 前に拒否する。

## 6. 検証方針

### 6.1 state・operator gate

- flatten/unflatten、deep copy、checkpoint text round-trip が bitwise 一致。
- 全 momentum の接線性 `|dot(m,k)|` が `64 epsilon max(|m|,m_scale)` 以下。
- signed incidence が `D2 D1 = 0` を全解像度・全cornerで exact に満たす。
- 定数 depth/zero velocity の RHS が丸め誤差 scale。
- arbitrary mass flux の全球 mass tendency が flux scaleの`5e-13`以下。
- manufactured `h,u` に対する pressure gradient、divergence、vorticity、vector
  Laplacian の誤差が `N=8,16,32,64` で低下する。
- seam/corner RMS が global RMS の3倍を継続的に超えない。

### 6.2 linear-wave・balance gate

- `Omega=0` の小振幅 inertia-gravity wave で位相速度が `sqrt(g H)` へ収束。
- spatial resolution と time step を別々に細分化し、空間/時間次数を分離する。
- 一様静止 state を1000 step保持し、mass drift `<=5e-14`、velocity は丸め誤差
  scale、depth の新極値なし。
- geostrophic adjustment で重力波と balanced mode を区別し、NaN/負 depth なし。
- Coriolisだけの inertial oscillation が SSP-RK3 の設計時間次数を示す。

### 6.3 Williamson test 2

test 2 の解析的な定常 zonal geostrophic flowを、回転軸と flow axis が一致する
場合および oblique に回転した場合で実行する。

- `N=12,24,48,96`、少なくとも5日積分。
- depth と velocity の `L1/L2/Linf`、balance residual、mass/energy/AAM drift。
- depth/velocity 誤差の観測次数は linear reconstruction で原則 `>=1.7`。
- mass drift `<=5e-13`、depth positive、全 state finite。
- oblique case の global error は axis-aligned case の2倍以内。
- seam/corner RMS は global RMS の3倍以内で解像度とともに低下。

### 6.4 Williamson test 6

Rossby–Haurwitz wave は厳密な時間依存 shallow-water 解として扱わず、公開初期条件
から生成した高解像度 reference と14日までの診断を比較する。

- day 1, 7, 14 の height/vorticity/PV norm と少数 field checksum。
- wave-4 の位相、振幅、zonal spectrum、grid-relative wave-4 contamination。
- mass、energy、potential enstrophy、AAM の時系列。
- reference は生成 compiler、scheme、resolution、time step、checksum を記録し、
  binary field を根拠なく Git へ固定しない。

### 6.5 Galewsky–Scott–Polvani jet

論文の balanced jet と局所 height perturbation を再現し、6日積分を主 gate とする。

- 初期 gradient-wind balance の数値積分誤差を、時間積分開始前に検査する。
- day 1–6 の relative vorticity/PV contour statistics、zonal spectrum、jet maximum、
  minimum depth、保存量を記録する。
- 少なくとも `N=48,96,192` で自己収束を確認し、最高解像度または公開 spectral
  referenceとの差を報告する。
- panel orientation を回転したcaseで instability onset、dominant wavenumber、
  conservation driftが同等であることを確認する。
- cube cornerだけのvorticity spike、seam反射、checkerboard/4-grid-waveを専用
  region/spectrum診断で検出する。

### 6.6 scheme比較 gate

基準 Rusanov と compatible 候補を同じ grid、initial state、CFL target、diffusion
e-folding timeで比較する。

| category | measure |
|---|---|
| accuracy | height/velocity/PV norm、phase、self-convergence |
| conservation | mass、energy、enstrophy、AAM driftとbudget residual |
| grid artifacts | seam/corner ratio、wave-4、checkerboard amplitude |
| robustness | minimum depth、non-finite、limiter/rejection count |
| cost | wall time、cell-step/s、edge flux/s、peak bytes/cell |

異なる diffusion を使って一方だけを有利にしない。安定化に必要な最小 diffusion
も比較結果の一部とする。P2.19 ADR は採用/不採用だけでなく、棄却した候補、既知の
欠点、Phase 5 で再評価する条件を記録する。

## 7. コミット列

### P2.01 — shallow-water 方程式・配置 ADR

```text
docs: define shallow-water state and discretization conventions
```

変更:

- `docs/adr/0002-shallow-water-state-and-staggering.md` を追加。
- 連続方程式、符号、`h/m/u/PV`、Cartesian tangent momentum、flux/source 分割を固定。
- 基準 collocated FV と比較対象 C-grid の責務境界を固定。

受け入れ条件: mass/momentum flux、Coriolis、pressure、sphere projection の式と
配置が一意で、test 2 の符号を手計算できる。

### P2.02 — shallow-water 設定 schema

```text
config: add shallow-water experiment schema
```

変更:

- `experiment.kind=shallow_water` と variant parameters を追加。
- scheme、reconstruction、limiter、CFL、depth、diffusion、diagnostic interval を
  strict parse/write。
- `configs/phase2_williamson2.cfg` を追加。

テスト: canonical round-trip、別 schema key、不正 depth/floor/CFL/係数を拒否し、
Phase 0/1 config fingerprint の回帰を維持する。

### P2.03 — shallow-water state と checkpoint layout

```text
dynamics: add versioned shallow-water state storage
```

変更:

- `ShallowWaterState`、depth、Cartesian momentum、primitive view を追加。
- flatten/unflatten と named checkpoint layout ID を追加。
- stage境界の有限性、positivity、接線性 validation を追加。

テスト: cell固有pattern、bitwise round-trip、shape/layout mismatch、radial momentum、
非有限値、負 depthを検出する。

### P2.04 — dual topology と incidence

```text
grid: build cubed-sphere dual topology and incidence
```

変更:

- cell-edge、edge-vertex incidence、vertex incident ring、dual area/lengthを構築。
- cube corner の3-valent dual cellを明示処理。

テスト: Euler count、ring orientation、positive dual metrics、`D2 D1=0`、全球 dual
area、`N=1..64` property gate。

### P2.05 — vector reconstruction と shallow-water operators

```text
numerics: add tangent-vector shallow-water operators
```

変更:

- Cartesian tangent vectorのleast-squares gradient/reconstructionを追加。
- edge normal/tangent projection、vorticity、PV、vector Laplacianを追加。

テスト: solid rotation、spherical harmonics/manufactured fields、constant/zero fields、
seam/corner convergence。

### P2.06 — shallow-water diagnostics

```text
diagnostics: add shallow-water invariants and balance budgets
```

変更:

- mass、energy、potential enstrophy、AAM、PV extrema、radial momentumを実装。
- flux、Coriolis、pressure、diffusion別 tendency/budget と deterministic reduction。

テスト: 解析的定数/solid rotation、回転軸、単位scale、budget component sum。

### P2.07 — Rusanov edge Riemann flux

```text
dynamics: add shallow-water Rusanov edge flux
```

変更:

- left/right stateをedge basisへ変換し、mass/normal/tangential momentum fluxを計算。
- wave speed、共通 Cartesian flux、left/right scatter APIを追加。

テスト: identical state、lake at rest、反転対称性、単一edge手計算、supercritical limit、
seamとpanel内edgeの同一路径。

### P2.08 — first-order global RHS

```text
dynamics: assemble conservative shallow-water tendencies
```

変更:

- first-order flux divergence、pressure、Coriolis、球面射影を合成。
- source/flux tendencyを個別取得できるdiagnostic APIを用意。

テスト: zero velocity/uniform depth、global mass cancellation、接線性、回転なし1D
dam-break symmetry、NaN propagation。

### P2.09 — gravity-wave CFL と driver

```text
dynamics: integrate shallow-water states with gravity-wave CFL
```

変更:

- edge wave speedからCFL stepを算出し、SSP-RK3へ接続。
- final-step調整、stop/restart、metadata、CLI dispatchを追加。

テスト: CFL上下、終了時刻、stage validation、連続/restart bitwise一致、sample smoke。

### P2.10 — depth-positive linear reconstruction

```text
dynamics: reconstruct bounded shallow-water face states
```

変更:

- tangent-plane linear reconstructionをdepth/velocityへ拡張。
- Barth–Jespersen bounds と face-depth positivity limiterを実装。

テスト: constant/linear field、smooth extrema、near-dry state、seam、mass不変性、
limiter activity count。

### P2.11 — linear inertia-gravity waves

```text
test: validate shallow-water wave propagation
```

変更:

- non-rotating linear-wave initial conditionとgeostrophic adjustment caseを追加。
- spatial/time refinementを分離したphase/amplitude診断。

受け入れ条件: section 6.2 のwave、rest、inertial oscillation gateを通す。

### P2.12 — Williamson test 2 initial state

```text
dynamics: add Williamson test 2 balanced flow
```

変更:

- standard constants、任意flow-axis回転、解析depth/velocityを追加。
- pressure/Coriolis balance residualを初期化時に検査。

テスト: pole/equator、axis rotation、解析mass/energy、初期PV、invalid constants。

### P2.13 — Williamson test 2 convergence

```text
test: validate steady geostrophic shallow-water flow
```

変更:

- `N=12,24,48,96` の5日 regression/convergence gate。
- height/velocity、balance、保存量、axis、seam/corner tableを生成。

受け入れ条件: section 6.3 を満たす。

### P2.14 — explicit diffusion と budget

```text
numerics: add budgeted shallow-water diffusion
```

変更:

- Laplacian/biharmonic scalar/vector diffusionとe-folding coefficient helper。
- diffusion energy/enstrophy tendencyを診断へ接続。

テスト: global mass-neutrality、constant/solid rotation、mode damping rate、符号、
resolution scaling、disabled時bitwise同一。

### P2.15 — Williamson test 6

```text
test: add Rossby-Haurwitz wave benchmark
```

変更:

- wave-4 initial condition、day 1/7/14 snapshot、reference generatorを追加。
- phase/spectrum/checksum と conservation regressionを追加。

受け入れ条件: section 6.4 のtableを再生成でき、解像度増加でreference差が低下する。

### P2.16 — Galewsky jet

```text
test: add unstable shallow-water jet benchmark
```

変更:

- balanced jet height integral、perturbation、6日driver、reference comparisonを追加。
- vorticity/PV spectrum、grid imprinting、minimum depth診断を追加。

受け入れ条件: section 6.5 を満たし、初期balance誤差とreference provenanceを記録。

### P2.17 — compatible C-grid diagnostic core

```text
numerics: add compatible shallow-water diagnostic operators
```

変更:

- edge-normal velocity、dual circulation/vorticity、cell kinetic energy、PV interpolation。
- non-orthogonal metric/Hodge operatorを実装。

テスト: exact incidence identities、constant PV、solid rotation、adjoint/skew symmetry、
cube corner、manufactured convergence。

### P2.18 — compatible PV/mass flux scheme

```text
dynamics: add compatible shallow-water flux scheme
```

変更:

- shared mass flux、Bernoulli gradient、PV flux/Coriolis operatorをdriverへ接続。
- positivity、diffusion、checkpoint layout、diagnosticsを基準法と共通化。

テスト: mass、constant PV、energy/enstrophy budget、Williamson 2 short gate、restart。

### P2.19 — scheme comparison ADR

```text
docs: select the horizontal shallow-water discretization
```

変更:

- `docs/adr/0003-shallow-water-horizontal-discretization.md` を追加。
- W2/W6/Galewsky の精度、保存、artifact、robustness、cost tableを比較。
- Phase 4/5で採用するscheme、fallback、既知の制約を決定。

受け入れ条件: section 6.6 の条件が同一で、結論を再生成できるraw diagnosticsと
commandが記録されている。

### P2.20 — CLI、出力、Phase 2 報告

```text
docs: record Phase 2 shallow-water validation
```

変更:

- machine-readable summary、interval diagnostics CSV、小規模field snapshotを完成。
- `docs/validation/phase-2.md`、README、全体計画の現在地を更新。
- clean-checkout validation commandとreference生成手順を固定。

受け入れ条件: 全CTest、Phase 0/1回帰、GCC/Clang、Debug/Release、ASan/UBSan、
formatがgreenで、公開した表をclean checkoutから再生成できる。

## 8. 依存関係

```text
P2.01 ─ P2.02 ─ P2.03 ─ P2.06
   │                 │
   ├──── P2.04 ─ P2.05 ─ P2.07 ─ P2.08 ─ P2.09 ─ P2.10
   │                                           │
   │                              P2.11 ─ P2.12 ─ P2.13
   │                                           │
   │                              P2.14 ─ P2.15 ─ P2.16
   │                                           │
   └────────────── P2.17 ─ P2.18 ─────────── P2.19 ─ P2.20
```

番号順を基本とする。dual topologyは基準法のcritical pathから外せるが、共有 header、
CMake、state layout の競合を避けるため P2.04 までに構造を固定する。高解像度
reference生成は P2.15/P2.16 の code review後に行い、未検証binaryを先に固定しない。

## 9. 各コミット共通の検証

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
cmake --build build/dev --target format-check
git diff --check
```

state、geometry、flux、time integration、checkpointへ触れるコミットは:

```sh
cmake --preset asan-ubsan
cmake --build --preset asan-ubsan
ctest --preset asan-ubsan
```

短い unit/property test と高解像度 validation を分け、後者は `phase2_gate` labelを
付ける。reference更新だけ、許容値緩和だけ、diffusion増加だけのコミットは禁止し、
数値根拠とresolution/budget tableを同時に記録する。

## 10. 完了チェックリスト

- [x] shallow-water state、配置、flux/source、Coriolis符号がADRで固定されている。
- [x] depth/momentum stateと全checkpoint layoutがbitwise round-tripする。
- [x] dual topologyとincidence identityが全panel/cornerで検証されている。
- [x] Rusanov unique-edge fluxが全球massを構造的に保存する。
- [x] uniform rest、linear wave、inertial oscillation、geostrophic adjustmentが通る。
- [x] face/update depth positivityをclipなしで維持する。
- [x] pressure、Coriolis、diffusionと保存量budget residualを分離できる。
- [~] Williamson test 2がaxis、収束、mass gateを通る。seam/corner gateは未実装。
- [~] Williamson test 6の14日runとwave-4診断を再生成できる。独立reference比較は未実施。
- [~] Galewsky jetの6日runと初期balance検査を再生成できる。自己収束/reference比較は未実施。
- [x] diffusionのmode damping、mass、energy/enstrophy budgetが検証されている。
- [x] compatible候補のtopology、PV flux、保存特性がunit testを通る。
- [x] 基準/compatible比較が同じCFL・diffusion条件で行われている。
- [x] Phase 4/5採用schemeとfallbackが結果ベースのADRで決定されている。
- [x] restartが連続実行とbitwise一致する。
- [x] Phase 0/1の全回帰、GCC/Clang、Release、ASan/UBSan、formatがgreen。
- [x] Phase 2検証報告と全tableをclean checkoutから再生成できる。

## 11. 一次資料

- D. L. Williamson et al., “A standard test set for numerical approximations to
  the shallow water equations in spherical geometry,” JCP 102, 211–224 (1992),
  <https://doi.org/10.1016/S0021-9991(05)80016-6>.
- J. Galewsky, R. K. Scott, and L. M. Polvani, “An initial-value problem for
  testing numerical models of the global shallow-water equations,” Tellus A 56,
  429–440 (2004), <https://doi.org/10.3402/tellusa.v56i5.14436>.
- T. Ringler et al., “A unified approach to energy conservation and potential
  vorticity dynamics for arbitrarily-structured C-grids,” JCP 229, 3065–3090
  (2010), <https://doi.org/10.1016/j.jcp.2009.12.014>.

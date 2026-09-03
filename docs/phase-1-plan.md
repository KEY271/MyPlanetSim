# Phase 1 実装計画 — cubed-sphere 幾何と球面輸送

## 1. 到達点

Phase 1 では全球 cubed-sphere の幾何と球面輸送を実装する。完了時には、任意軸の prescribed wind による保存型 scalar transport を6面全体で実行し、標準ケースの収束・質量保存・極値・seam/corner誤差を自動評価できることを目標とする。

成果物は次を含む。

- equiangular gnomonic cubed-sphere の座標写像、逆写像、panel topology
- セル中心・頂点・面積・辺長・辺接法線・局所基底を持つ全球 geometry
- 全球で一意な共有辺と、隣接セルへ反対符号で適用される flux
- panel-local scalar field と edge halo exchange
- gradient、divergence、curl、Laplacian の基礎離散演算子
- prescribed wind と conservative tracer equation
- 一次風上参照方式、線形再構築、limiter、SSP-RK3、CFL制御
- Williamson test 1、solid-body rotation、Lauritzen et al.変形流
- mass、`L1/L2/Linf`、min/max、filament、mixing、seam/corner診断
- transport用CLI、設定、checkpoint/restart、Phase 1検証報告

3次元Cartesianベクトルは球面を埋め込む共通表現として使用する。

## 2. Phase 1に含めないもの

- shallow-waterの運動量、自由表面、Coriolis力（Phase 2）
- 地形、鉛直座標、3D大気
- MPI/OpenMP/GPU並列化
- NetCDFなどの大規模field形式
- 高次DG/SE、semi-Lagrangian、semi-implicit方式の比較

既存の`Field2D` periodic haloはstorage/boundaryの回帰テスト用utilityとして残すが、Phase 1の科学ケースには使わない。

## 3. 固定する数値設計

### 3.1 panel座標と写像

各panelの論理座標を

```text
-pi/4 <= alpha <= pi/4
-pi/4 <= beta  <= pi/4
```

とし、各方向を`N`セルへ等分する。panelごとの中心法線`n`と接方向`e_alpha`, `e_beta`から単位球面位置を

```text
x = normalize(n + tan(alpha) * e_alpha + tan(beta) * e_beta)
```

で定義する。常に`e_alpha cross e_beta = n`となる右手系を使う。P1.01のADRでは以下を基準案として固定する。

| Panel | `n` | `e_alpha` | `e_beta` |
|---|---|---|---|
| `PX` | `(1,0,0)` | `(0,1,0)` | `(0,0,1)` |
| `PY` | `(0,1,0)` | `(-1,0,0)` | `(0,0,1)` |
| `NX` | `(-1,0,0)` | `(0,-1,0)` | `(0,0,1)` |
| `NY` | `(0,-1,0)` | `(1,0,0)` | `(0,0,1)` |
| `PZ` | `(0,0,1)` | `(1,0,0)` | `(0,1,0)` |
| `NZ` | `(0,0,-1)` | `(1,0,0)` | `(0,-1,0)` |

逆写像は最大絶対Cartesian成分からpanelを選ぶ。共有辺・cube corner上のtie-breakは診断用にだけcanonical順を定義し、topologyやflux接続を逆写像のtie-breakへ依存させない。

### 3.2 geometry表現

`CubedSphereGrid(N, radius_m)`は不変なgeometryを所有する。

- cellは`(panel, i, j)`で識別し、flatten順はpanel-major、次に`j`, `i`とする。
- cell cornerは解析写像から求める。
- 球面四角形面積は符号付き球面三角形2枚のsolid angleから求める。
- 辺長はgreat-circle arc lengthとする。
- edge centerは端点の正規化和、edge tangent normalは球面接平面内でleft cellからright cellへ向ける。
- 全球edgeは一度だけ保存し、closed sphereなので各edgeは必ず2セルを持つ。

`N`に対してcell数は`6N^2`、unique edge数は`12N^2`、unique vertex数は`6N^2+2`になる。これをtopologyの構造テストに使う。

### 3.3 scalar fieldとhalo

tracerはcell average `q`として6個の`Field2D<Real>`に保持する。Phase 1の再構築はface-neighborだけを使うため、panel edge haloだけを交換し、物理的に一意でないcorner ghostを計算へ使わない。corner haloはDebug時にpoison値のままとし、誤参照を検出する。

面をまたぐscalarは値をコピーできるが、`i/j`の入替えと反転は`PanelTopology`が明示する。将来のvector haloで同じ写像を再利用できるAPIにする。

### 3.4 共有辺flux

予報量はcell tracer mass

```text
M_c = area_c * q_c
```

とし、各unique edgeで面積流量`F_e`とupwind face valueを一度だけ計算する。同じ`F_e*q_face`をleft cellへ負、right cellへ正で加える。これによりpanel seamを含む全球質量保存を構造的に保証する。

prescribed wind APIの主要量はpoint velocityではなく、向き付きedge-integrated flux `F_e(t)`とする。非発散風ではvertex streamfunction差からfluxを構成できる経路を用意し、cellごとの離散divergenceを抑える。一般の風は同一edge上の共通quadratureで積分する。

### 3.5 輸送方式

1. 参照方式: piecewise constant upwind flux。低精度だが保存性、向き、CFL、面接続の比較基準とする。
2. 目標方式: cell-centered least-squares gradientによる線形再構築。球面接平面の物理変位を使い、panel-local座標差を直接使わない。
3. limiter: Barth–Jespersen型を初期候補とし、face-neighborの局所min/max内へ再構築値を制限する。滑らかな収束試験ではunlimitedとlimitedを分けて報告する。
4. time integration: Phase 0のSSP-RK3を使用する。Forward Eulerはflux/RHSの局所テストだけに使う。
5. time step: 各cellのoutgoing edge fluxから求める有限体積CFLと、設定上の最大`dt`の小さい方を使う。最終stepは終了時刻へ揃える。

### 3.6 設定schema

Phase 1で`experiment.kind`を導入し、Phase 0 ODEとsphere transportのschemaを区別する。transport設定は最低限次を含む。

```text
experiment.kind = sphere_transport
grid.cells_per_panel
grid.halo_width
transport.test_case
transport.initial_condition
transport.scheme
transport.limiter
transport.cfl
transport.rotation_axis_x
transport.rotation_axis_y
transport.rotation_axis_z
transport.angular_speed_rad_s
```

schemaごとに必須・許可keyを変え、未知keyや別experiment用keyを拒否する。既存`phase0_ode.cfg`とcheckpoint回帰は維持する。

## 4. 検証方針

### 4.1 geometry gate

- 全位置vectorの長さ誤差が`32 * epsilon`以下。
- cell interiorの写像→逆写像誤差が`128 * epsilon rad`以下。
- 共有vertex位置と共有edge端点が向きによらず同一許容誤差内。
- 全area/edge lengthが正で有限。
- 全球面積和の`4*pi*R^2`に対する相対誤差が`5e-13`以下。
- topology数が`cells=6N^2`, `edges=12N^2`, `vertices=6N^2+2`と一致。
- 各edgeがちょうど2セルを持ち、left/right法線規約が反対称。

`N=1, 2, 4, 16, 64`をunit/property testで使い、cube cornerを必ず含める。

### 4.2 operator gate

- 定数scalarのgradientが丸め誤差scale。
- 任意の共有edge fluxについて全球`sum(area * divergence)`がflux scaleの`5e-13`以下。
- smooth manufactured fieldのgradient/divergence/Laplacian誤差が解像度とともに減少。
- `curl(gradient)`などの恒等式は、採用離散化でexactなものと収束で0になるものを区別して報告。
- panel seam/corner近傍のRMS誤差が全球RMSの3倍を継続的に超えず、解像度とともに減少。

### 4.3 transport gate

- 全標準ケースでNaN/Infを出さない。
- unique edge flux方式の全球tracer mass driftが標準ケースで`5e-13`以下。
- 非発散solid-body windで定数tracerを許容誤差内に保持。
- smooth Gaussian tracerの一周後誤差が、一次方式で観測次数`>=0.8`、unlimited線形再構築で`>=1.7`。
- limiter有効時に初期/global boundaryが許す範囲を`1e-12`超える新極値や負値を作らない。
- 回転軸を座標軸、斜め軸へ変えても収束次数と質量保存が同じgateを通る。
- restart後のflattened tracer stateが連続実行とbitwise一致。

Williamson cosine bellは微分可能性が低いため二次収束gateには使わず、公開診断との回帰比較に使う。Lauritzen suiteでは通常normに加え、filament preservationとunmixing/overshooting/real mixing診断を記録する。

## 5. コミット規約

Phase 0と同じく、各コミットは一つの関心事と対応テストを含み、`dev` build、全CTest、format-checkを通す。geometry、index、flux、time integration、I/Oへ触れるコミットはASan/UBSanも通す。許容値変更だけのコミットは禁止し、数値根拠とresolution tableを同時に記録する。

## 6. コミット列

### P1.01 — panel規約ADR

```text
docs: fix cubed-sphere panel and orientation conventions
```

変更:

- `docs/adr/0001-cubed-sphere-panel-conventions.md`を追加。
- panel enum、basis table、alpha/beta方向、edge名、left/right、tie-breakを固定。
- 数式と小さな手計算例で隣接面の向きを確認。

受け入れ条件: 全12 cube edgeの隣接panel、index入替え、反転有無が表で一意に決まる。

### P1.02 — experiment別設定schema

```text
config: add cubed-sphere transport experiment schema
```

変更:

- `experiment.kind`とODE/transportのvariant設定を追加。
- grid、scheme、limiter、CFL、test case、回転軸をstrict parse/write。
- `configs/phase1_solid_body.cfg`を追加。

テスト:

- ODE設定の後方回帰。
- transport設定のcanonical round-trip。
- schema外key、`N<=0`、不正軸、`0<CFL<=1`違反、不明schemeを拒否。

### P1.03 — 3D vectorと球面基本演算

```text
geometry: add Cartesian vector and spherical primitives
```

変更:

- `Vec3`、dot、cross、norm、normalize、接平面射影、safe angleを追加。
- componentアクセスと演算はconstexpr/noexcept可能な範囲を明示。

テスト: 直交基底、右手系、正規化、反平行vector、zero vector拒否、NaN拒否をGCC/Clangで確認。

### P1.04 — gnomonic panel写像

```text
geometry: map equiangular panel coordinates to the sphere
```

変更:

- `Panel`, `PanelBasis`, forward/inverse mapping、局所接基底を追加。
- panel interiorとcanonical boundary tie-breakを分離。

テスト: 中心・四隅・edge、forward/inverse round-trip、単位長、basis接線性とhandednessを全panelで確認。

### P1.05 — panel topology

```text
grid: encode cubed-sphere panel topology
```

変更:

- 4 edge × 6 panelのadjacency、neighbor edge、index swap/reversalを追加。
- edge indexの往復写像APIを追加。

テスト: 全接続の往復が恒等、共有点が同じCartesian位置、各cube edgeが2回だけ参照される。

### P1.06 — cell/edge geometry

```text
geometry: build cubed-sphere cell and edge metrics
```

変更:

- cell center/corner/area、edge center/length/tangent normalを生成。
- `N`, `radius_m`、overflow、有限性を検証。

テスト: 正値性、半径scale、panel対称性、共有edge一致、`N=1` corner caseを確認。

### P1.07 — 全球topologyとgeometry gate

```text
test: enforce cubed-sphere geometry invariants
```

変更:

- unique vertex/edge構築とEuler count property test。
- 全球面積、写像誤差、共有位置、法線向きのmulti-resolution gate。
- failure時にpanel/index/誤差表を表示。

受け入れ条件: 4.1の全geometry gateを`N=1..64`で通す。

### P1.08 — cubed-sphere scalar field

```text
grid: add panel-based cubed-sphere scalar fields
```

変更:

- 6 panelの`Field2D<Real>`を所有する`CubedSphereField`を追加。
- panel-major flatten/unflatten、interior iteration、fill/copyを追加。
- checkpoint順をversion付きで固定。

テスト: 全cell固有値、flatten round-trip、deep copy、shape mismatch拒否、checkpoint state sizeを確認。

### P1.09 — scalar edge halo exchange

```text
grid: exchange scalar halos across panel edges
```

変更:

- `PanelTopology`に従うedge halo exchangeを追加。
- corner ghostは使用禁止状態を維持。

テスト: 全panel/edgeへ固有patternを置き、halo幅1/2、swap/reversal、interior不変性を全要素比較。

### P1.10 — vector basis変換

```text
geometry: transform tangent vectors across panel bases
```

変更:

- Cartesian tangent vectorとpanel contravariant/covariantまたは物理局所成分の変換を追加。
- 採用成分をAPI名で区別し、暗黙変換を禁止。

テスト: round-trip、接線性、norm、共有edge上の法線/接線成分連続性、極・corner近傍を確認。

### P1.11 — unique edge graphと共有flux accumulator

```text
numerics: assemble unique cubed-sphere edge fluxes
```

変更:

- edge ID、left/right cell、orientation、flux accumulatorを追加。
- 同一fluxを隣接2セルへ反対符号でscatterする。

テスト: edge count、全cell connectivity、任意fluxの全球cancellation、seamとpanel内edgeの同一コード経路を確認。

### P1.12 — 基礎離散演算子

```text
numerics: add cubed-sphere differential operators
```

変更:

- finite-volume divergence、least-squares scalar gradient、curl、Laplacianを追加。
- operatorごとに必要なstencil/haloを明示。

テスト: 定数、解析的球面調和関数またはmanufactured Cartesian polynomial、全球積分、multi-resolution収束を確認。

### P1.13 — prescribed windと初期tracer

```text
transport: add spherical wind and tracer test cases
```

変更:

- solid-body rotationのCartesian velocity/streamfunction/edge fluxを追加。
- constant、Gaussian hill、cosine bell、slotted cylinder初期条件を追加。
- 軸vectorを正規化し、0 vectorを拒否。

テスト: 接線風、角速度、周期、非発散flux、回転軸変更、初期mass/min/maxを確認。

### P1.14 — 一次保存型輸送RHS

```text
transport: add first-order conservative upwind fluxes
```

変更:

- cell average→mass、upwind edge value、unique edge flux、mass tendencyを追加。
- 定数保持とmass budgetを診断へ接続。

テスト: zero wind、符号反転、単一edge手計算、panel seam横断、global tendency sumを確認。

### P1.15 — CFLと全球transport driver

```text
transport: integrate cubed-sphere tracer states
```

変更:

- outgoing fluxに基づくcell CFL、最大`dt`、最終step調整を追加。
- Phase 0 SSP-RK3、metadata、checkpoint/restartへflattened stateを接続。

テスト: CFL上下、停止時刻、NaN/負面積拒否、連続/restart bitwise一致、mass driftを確認。

### P1.16 — solid-body rotation回帰

```text
test: validate first-order spherical solid-body transport
```

変更:

- smooth GaussianとWilliamson test 1の一周実験を追加。
- `N=12,24,48,96`の誤差・mass・axis別診断を自動化。

受け入れ条件: 一次方式の観測次数`>=0.8`、mass gate、全軸で有限・安定、seam誤差が解像度とともに減少。

### P1.17 — 線形再構築

```text
transport: add tangent-plane linear reconstruction
```

変更:

- face-neighborを用いるleast-squares gradientと左右face値を追加。
- rank不足、condition不良を明示検出。

テスト: 定数・接平面線形manufactured field、panel seam、回転対称ケース、unlimited smooth transport収束を確認。

### P1.18 — bounded limiter

```text
transport: bound reconstructed tracer states
```

変更:

- Barth–Jespersen型limiterとpositivity optionを追加。
- limiter activity countを診断へ出す。

テスト: smooth extrema、step/slotted cylinder、負値防止、定数保持、mass不変性を確認。

### P1.19 — Lauritzen変形流suite

```text
test: add deformational spherical transport benchmarks
```

変更:

- 非発散/発散変形流、smooth/rough tracerを追加。
- filament preservationとmixing分類診断を実装。
- resolution tableと回帰許容値を設定。

受け入れ条件: mass/bounds gateを通し、smooth caseの設計収束を示し、filament/mixing曲線を再生成できる。

### P1.20 — CLI、出力、Phase 1報告

```text
docs: record Phase 1 spherical transport validation
```

変更:

- applicationを`experiment.kind`でODE/transportへdispatch。
- machine-readable diagnostic summaryと小規模CSV snapshotを追加。
- `docs/validation/phase-1.md`へtoolchain、geometry、operator、収束、保存、axis、seam、Lauritzen結果を記録。
- READMEと全体計画の現在状態を更新。

受け入れ条件: clean checkoutからsample実行と全表を再生成でき、全CTest・ASan/UBSan・format・GCC/Clang matrixがgreen。

## 7. 依存関係

```text
P1.01 ─ P1.03 ─ P1.04 ─ P1.05 ─ P1.06 ─ P1.07
                  │                 │
P1.02 ────────────┼──── P1.08 ─ P1.09 ─ P1.10
                                    │
                         P1.11 ─ P1.12
                                    │
                         P1.13 ─ P1.14 ─ P1.15 ─ P1.16
                                             │
                                   P1.17 ─ P1.18 ─ P1.19 ─ P1.20
```

原則は番号順に統合する。geometryと設定schemaは独立に実装可能だが、共有headerやCMake競合を避けるため履歴上は直列化する。

## 8. 各コミット共通の検証

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev
cmake --build build/dev --target format-check
git diff --check
```

geometry、field、flux、time integration、checkpointへ触れる場合:

```sh
cmake --preset asan-ubsan
cmake --build --preset asan-ubsan
ctest --preset asan-ubsan
```

multi-resolution testは通常unit testと分けて`phase1_gate` labelを付ける。短いCI gateと高解像度validationを分ける場合も、両者の設定と期待値を同じ文書に記録する。

## 9. 完了チェックリスト

- [ ] panel規約と全edge接続がADRで固定されている。
- [ ] geometry gateが`N=1..64`で合格する。
- [ ] scalar/vector変換とedge haloが全panel、seam、corner近傍で検証されている。
- [ ] unique edge fluxの全球cancellationが構造テストを通る。
- [ ] 基礎離散演算子が解析場へ収束する。
- [ ] 一次upwindがsmooth tracerで観測次数`>=0.8`を示す。
- [ ] unlimited線形再構築がsmooth tracerで観測次数`>=1.7`を示す。
- [ ] limiterがmassを変えず、新極値と負値を許容幅内に抑える。
- [ ] Williamson test 1とLauritzen suiteが再現可能な設定で実行できる。
- [ ] 全標準ケースのmass driftが`5e-13`以下である。
- [ ] 回転軸変更時も収束・保存・seam gateを満たす。
- [ ] checkpoint/restartが連続実行とbitwise一致する。
- [ ] GCC/Clang、Debug/Release、ASan/UBSan、formatがgreen。
- [ ] Phase 1検証報告をclean checkoutから再生成できる。
- [ ] Phase 2 shallow-waterが共有geometry/edge fluxを再利用できるAPIになっている。

# Phase 6 実装計画 — 固定地形と下部境界

**状態:** 設計済み・未実装。実装開始前に [Phase 5 検証報告](validation/phase-5.md) の
定量的 convergence と UMJS14 reference-envelope gap を閉じ、実装が
[ADR 0006](adr/0006-dry-hydrostatic-state-and-coupling.md) の SSP-RK3、再構築、CFL 契約と
一致することを再確認する。Phase 6 の判断は
[ADR 0008](adr/0008-fixed-orography-and-lower-boundary.md) に提案状態で記録する。

## 1. 到達点と開始条件

Phase 6 は Phase 5 の乾燥・静水圧 3D core に、時間不変の cell-centred surface
geopotential と不透過な地形追随下端を加える。既存 generalized pressure-coordinate
pressure gradient を傾いた hybrid 面で評価し、静止大気、標準山岳付き shallow-water、
乾燥 baroclinic case で誤差と budget を測る。

```text
Phase 5 entry gate
  quantitative 3D convergence + UMJS14 reference envelope
  + ADR 0006 implementation audit
                    |
                    v
Milestone A: fixed-terrain C++ contract
  orography field -> dry diagnosis/source -> flat equivalence
                    |
                    v
Milestone B: analytic benchmark gate
  DCMIP 2-0-0 -> Williamson 5 -> linear hydrostatic wave -> JW06
                    |
                    v
Milestone C: one imported-data path
  strict lat-lon CSV -> interpolation -> bounded seam-safe smoothing
  -> full C++ validation report
```

Phase 5 entry gate は Phase 6 の成果に数えない。未完了の既存 solver を terrain 固有の補正で
覆わず、前段の責務として直す。特に現在の validation report が bounded
initialization/regression としている DCMIP/UMJS case を、Phase 6 完了の根拠へ読み替えない。

Phase 6 の成果物は次に限定する。

- `6*N*N` の不変な surface-geopotential field と出典 fingerprint
- 既存 hydrostatic column への cell 別 lower-boundary geopotential の供給
- 既存 pressure-gradient identity の sloping-hybrid-surface 評価
- Rusanov shallow-water path の地形 source と Williamson test 5
- DCMIP 2-0-0、短時間の線形 hydrostatic mountain wave、JW06 preset
- pressure-gradient error、spurious wind、mountain torque、terrain pressure work の診断
- strict rectilinear lat-lon CSV 一形式の補間と、一種類の bounded smoothing
- C++ unit/property/convergence/regression と `docs/validation/phase-6.md`

## 2. スコープ境界

### 2.1 含めるもの

- 時間不変、cell-centred、scalar の `surface_geopotential_m2_s2[C]`
- flat、DCMIP 2-0-0、Williamson 5、linear bell、JW06 の named analytic field
- `Phi_s=g*z_s` の単位変換と surface-pressure/layer-positivity validation
- Phase 4 hydrostatic integral、Phase 5 state/continuity/vertical transport の再利用
- Rusanov shallow-water momentum への `-h*grad(Phi_s)` source
- config-relative local CSV、content fingerprint、periodic-longitude bilinear interpolation
- unique-edge を使う conservative/convex smoothing と設定に残る pass count
- CLI、CSV diagnostics、metadata、checkpoint/restart の terrain-aware regression

### 2.2 含めないもの

- Phase 5 の未完了 convergence を terrain test で代替すること
- `DryHydrostaticStateV2`、terrain prognostic、checkpoint layout 変更
- moving/deforming terrain、land/sea mask、surface pressure の prescribed forcing
- surface drag、heat/moisture flux、boundary layer（Phase 7）
- vertical momentum、nonhydrostatic pressure、acoustic mode
- DCMIP 2-0-1、2-1、2-2 の conformance、Rayleigh sponge
- exact-rest 専用 perturbation solver、benchmark 固有 source cancellation
- compatible shallow-water candidate の terrain 対応
- NetCDF、GeoTIFF、GDAL、online data download、欠損値補完、複数 interpolation 法
- adaptive/spectral filter、subgrid orography、gravity-wave drag parameterization
- `FrameV3`、terrain mesh、Web/gateway/control protocol の変更

DCMIP 2-1/2 は reduced-size planet で非静水圧応答と sponge layer を要求するため、乾燥・静水圧で
surface physics を持たない Phase 6 の検証には使わない。線形山岳波は full-size planet の
hydrostatic long-wave regime に限定し、上端反射が診断領域へ戻る前に終了する。

## 3. state と固定地形 contract

### 3.1 run-owned field

地形は atmospheric state から分離し、driver が次を一つ所有する。

```text
SurfaceOrography
  surface_geopotential_m2_s2[C]
  source_fingerprint

surface_height_m[c] = surface_geopotential_m2_s2[c] / gravity_m_s2
```

順序は既存の `panel-major, then j, then i` とする。shape、finite、source fingerprint を state
生成前に検査する。benchmark/imported field は非負とするが、共通 container は reference datum の
違いを妨げないよう finite のみを要求する。

地形は固定 forcing/configuration なので `DryHydrostaticState`、`ShallowWaterState`、checkpoint
payload へ複製しない。restart は同じ canonical config と terrain input fingerprint から field を
再生成する。Phase 4 の scalar `vertical.surface_geopotential_m2_s2` は column fixture のままにし、
全球地形 amplitude として再利用しない。

### 3.2 config

既定値は `flat` とし、既存 Phase 0--5 config に新しい必須 key を加えない。flat default の
canonical writer 出力も変えず、既存 fingerprint を保つ。非 flat のときだけ次を canonical text と
metadata に書く。

```text
orography.kind = dcmip_2_0_0 | williamson5 | linear_bell | jw06 | latlon_csv
orography.input_file = data/orography/example.csv
orography.input_fingerprint_fnv1a64 = 0123456789abcdef
orography.smoothing_passes = 0
```

`input_file`、fingerprint、smoothing は `latlon_csv` のみで必須/有効とする。named analytic field は
論文値を initializer に固定し、smoothing を拒否する。山高、幅、中心、taper を独立 knob にせず、
benchmark を変更したい場合は別 test case として明示する。

relative path は config file の親 directory から解決する。input bytes の FNV-1a 64-bit 値を load
前に照合し、期待値を canonical config に含める。これは再現性の accidental-change detector であり、
暗号学的 integrity は主張しない。

## 4. 乾燥 hydrostatic core への接続

### 4.1 hydrostatic diagnosis

`diagnose_dry_hydrostatic_state` は `surface_geopotential_m2_s2[C]` を受け、cell `c` の既存 column
積分を次の lower boundary から開始する。

```text
Phi_half[c,K] = Phi_s[c]
Phi_half[c,k], Phi_full[c,k] = existing Phase 4 hydrostatic integral
```

`p_half=A+B*ps`、`M=Delta p/g`、Exner、temperature の定義は変えない。pressure/geopotential を
state に保存せず、各 RHS stage の現 state と固定 field から再診断する。

### 4.2 pressure gradient

運動量 source は ADR 0006 のままとする。

```text
G[c,k] = grad_eta(Phi_full[:,k])[c]
       + (Rd*T[c,k]/p_full[c,k]) * grad_eta(p_full[:,k])[c]
S_pg[c,k] = -M[c,k] * project_tangent(G[c,k])
```

二つの gradient は同じ `least_squares_gradient` を使う。`grad(Phi_s)` を別 source としてもう一度
加えない。既存 `Phi_full` が lower-boundary shift を含むためである。

DCMIP 2-0-0 で cancellation error が設計次数どおり減らない場合は、case 固有補正や速度 clip を
入れない。誤差 component と conditioning を記録し、perturbation form などの代替は別 ADR で比較する。

### 4.3 lower boundary と continuity

地形は固定であり、top/bottom coordinate-relative mass flux は引き続きゼロとする。horizontal layer
convergence から `ps_dot` と `F` を作る Phase 4 recurrence は変更しない。地表に沿う水平流は
geometric vertical velocity を持ち得るが、terrain-following coordinate に対する flux はゼロである。

absolute interface pressure velocity が山岳波診断に必要な場合だけ、既存 RHS quantity から

```text
omega_abs = B*ps_dot + u_half dot grad_eta(p_half) + g*F_relative
```

を diagnostic として作る。`w` を予報せず、`omega_relative=g*F` の既存意味も変えない。

## 5. shallow-water 回帰

Williamson test 5 の initial state は test 2 と同じ zonal flow/depth 規約を再利用し、published conical
mountain を `SurfaceOrography` として構築する。Rusanov RHS に独立 component を一つ加える。

```text
S_orography[c] = -depth[c] * grad(Phi_s)[c]
total.momentum += S_orography
```

mass flux、wave-speed、depth checkpoint は変えない。flat field では source を exact zero にする。
Phase 2 で不採用となった compatible path は変更せず、non-flat + compatible を config validation で
明示的に拒否する。

Williamson 5 は解析解がないため、画像の見た目だけを gate にしない。day 10/15 の free-surface
geopotential、zonal/meridional wind について出典付き spectral reference/envelope を用い、area-weighted
norm、extrema、主要波の位置を記録する。reference data が用意できるまで「実行できる」ことだけで
完了扱いしない。

## 6. 地形生成、補間、平滑化

### 6.1 analytic fields

次の profile を named initializer として、文献の座標・定数・great-circle distance 式から直接作る。

- `dcmip_2_0_0`: moderately steep Schar-like circular mountain
- `williamson5`: isolated conical mountain
- `linear_bell`: small-amplitude smooth bell for the hydrostatic linear-response gate
- `jw06`: Jablonowski--Williamson smooth surface geopotential

中心を回転できる内部 test helper は持つが、production config の自由 knob にはしない。rotation test は
同じ scalar profile を cell centre で再評価し、panel interior、seam、corner 近傍に中心を置く。

### 6.2 one external format

外部入力は次の exact header を持つ rectilinear lat-lon CSV 一形式だけとする。

```text
longitude_deg,latitude_deg,height_m
```

longitude は一周期を cover し、latitude と longitude の全 tensor-product pair が一度ずつ存在し、軸は
strict monotone、値は finite とする。duplicate、missing point、non-finite、単位違いを拒否し、暗黙の
穴埋めをしない。longitude は periodic bilinear interpolation、latitude は bounded interpolation とし、
cell centre で `height_m` を標本化して `g` を一度だけ掛ける。

### 6.3 one filter

optional smoothing は cubed-sphere の unique edge を一度ずつ走査する synchronous diffusion pass
一種類だけとする。edge `i-j` で equal-and-opposite な area-integrated increment を交換し、global
coefficient は全 cell で convex combination になる上限の固定 fraction とする。

各 pass は次を property test で保証する。

- constant field を exact に維持する。
- area-weighted mean を rounding 内で維持する。
- input min/max を越えない。
- panel seam/corner を interior edge と同じ経路で処理する。

configurable parameter は integer `smoothing_passes` だけとし、kernel、coefficient、cutoff の組合せを
増やさない。analytic benchmark は filter を通さない。

## 7. benchmark

### 7.1 flat/uniform regression

最初の gate は新しい物理 case ではなく Phase 5 identity とする。

- `Phi_s=0` で derived state、RHS component、stable dt、step 後 state、CSV、checkpoint が既存 fixture
  と exact に一致する。
- uniform nonzero `Phi_s` で全 geopotential が同じ定数だけ shift し、pressure/Coriolis/transport
  tendency は rounding 内で不変である。
- potential-energy difference は `Delta Phi * total dry mass` と一致する。

この gate が通るまで analytic mountain initializer を driver に接続しない。

### 7.2 DCMIP 2-0-0

published lapse-rate atmosphere、terrain-dependent surface pressure、zero wind、zero rotation を用いる。
6 日積分で少なくとも次を出す。

- level/column ごとの pressure-gradient cancellation `L1/L2/Linf`
- maximum horizontal wind と global mean kinetic energy の時系列
- surface pressure、temperature、hydrostatic residual
- panel/seam/corner 別 error

水平 `N` を 3 点、鉛直は少なくとも `K=15,30` で独立に refinement する。開始前に baseline data から
observed-order と absolute ceiling を validation report に固定し、途中で threshold を緩和しない。
DCMIP 文書が述べる「特殊な perturbation form なら静止を構造的に保つ」方式は採用しないため、
machine-zero ではなく収束する truncation error を評価する。

### 7.3 linear hydrostatic mountain response

non-rotating full-size planet、低い smooth bell、長い水平波長、弱い一様 background flow を用いる。
山高を 1/2 にしたときの pressure/wind response、dominant wavelength/phase、pressure-coordinate momentum
flux の線形性を調べる。run は上端反射が mountain 周辺の診断域へ戻る前に終了する。

この case は hydrostatic long-wave regression であり DCMIP 2-1/2 conformance ではない。Rayleigh
sponge、vertical mixing、nonhydrostatic reference を追加しない。

### 7.4 Williamson 5 と rotated terrain

標準 conical mountain を Rusanov shallow-water で day 15 まで進め、day 10/15 の reference envelope、
mass/energy/AAM drift、orographic source budget と比較する。続いて中心を panel interior、seam、corner
近傍へ剛体回転した小規模 case を実行し、回転後に戻した field の norms、extrema、budget が同程度で
あることを確認する。

### 7.5 JW06

Phase 5 の UMJS14 flat case が quantitative gate を通過した後にだけ JW06 steady/perturbed state を
追加する。published smooth surface geopotential、surface pressure、temperature、wind を一体の named
initializer とし、個別 knob で組み替えない。

- unperturbed state: zonal symmetry、`ps/u/T` error、energy/AAM/torque、grid imprint の refinement
- perturbed state: published diagnostic days の `ps`、low-level temperature、vorticity の
  extrema/location/norm を reference envelope と比較

JW06 が失敗したとき UMJS14 も同時に失敗するなら Phase 6 terrain defect と分類せず、共有 Phase 5
dynamics の regression として切り分ける。

## 8. diagnostics と I/O

既存 dry diagnostics に次を additive に出す。

- min/max `surface_height_m`、max resolved `|grad(Phi_s)|/g`
- pressure-gradient cancellation norm、max spurious wind、kinetic energy
- model pressure-gradient axial torque、boundary mountain torque、その residual
- terrain pressure-work conversion と total pressure-work component
- mountain-wave `omega_abs` と pressure-coordinate momentum flux（該当 case のみ）

mountain torque と terrain work は force を二重に加える source ではない。固定・不透過地形から外部へ
moving-boundary work は入らず、terrain work は kinetic と potential/internal energy 間の conversion と
明記する。

state CSV は既存列を維持し、non-flat run に限って各行に derived
`surface_geopotential_m2_s2` と `surface_height_m` を additive に出す。diagnostics CSV は unit
suffix 付き列を追加する。地形専用 binary、volume format、frame schema は作らない。

non-flat run の metadata は orography kind、input relative path、input fingerprint、smoothing pass、
min/max/slope を記録する。flat run の metadata schema と canonical config は変えない。checkpoint
payload/layout は不変で、restart 時に再生成した field と fingerprint を検査する。

## 9. 検証 gate

### 9.1 unit/property

- field shape/order/finiteness、`Phi_s=g*z_s`、invalid config rejection
- analytic profile の中心、高さ、support、great-circle rotation
- flat exact identity と uniform-offset invariance
- hydrostatic integral が各 cell の `Phi_s` から始まること
- sloping manufactured state の pressure-gradient horizontal convergence
- shallow-water flat source exact zero、non-flat sign/unit/tangency
- lat-lon periodic interpolation、duplicate/missing/non-finite/fingerprint rejection
- smoothing の constant/area mean/extrema/seam/corner property
- checkpoint restart が同じ terrain で一致し、changed fingerprint を拒否

### 9.2 convergence/regression

- DCMIP 2-0-0 の horizontal/vertical independent refinement と 6-day spurious-wind/KE history
- Williamson 5 day 10/15 reference envelope、budget、Rusanov-only validation
- linear hydrostatic response の amplitude linearity、wavelength/phase/momentum flux
- JW06 steady refinement と perturbed reference envelope
- terrain centre の panel interior/seam/corner rotation sensitivity
- mass/theta/tracer、total energy、AAM、mountain torque、pressure-work residual
- all Phase 0--5 tests、restart、GCC/Clang、ASan/UBSan、format-check

threshold は analytic error、published reference envelope、または 3-resolution observed order から
事前登録する。小さい CI fixture と production-resolution validation を区別し、前者だけの成功を
Phase 6 完了としない。結果、command、cost、known gap を `docs/validation/phase-6.md` に記録する。

## 10. コミット列

一つの commit は一つの数値責務とその test を含む。Milestone A の flat identity が通る前に外部 parser
や benchmark を追加しない。`web/` は全列で変更しない。

1. `P6.01 docs: define the fixed-orography contract`
   — ADR 0008 で immutable field、lower boundary、pressure source、I/O boundary を固定。
2. `P6.02 config: add bounded orography selection`
   — flat default、named kind、strict conditional keys、canonical-text regression。
3. `P6.03 dynamics: apply cellwise hydrostatic lower boundaries`
   — driver-owned field、column integration、flat/uniform-offset property tests。
4. `P6.04 dynamics: validate pressure gradients on sloping hybrid surfaces`
   — manufactured field、component diagnostics、horizontal/vertical refinement fixture。
5. `P6.05 test: add the DCMIP 2-0-0 resting atmosphere`
   — published terrain/state、spurious wind/KE、rotation and resolution tests。
6. `P6.06 dynamics: add the selected shallow-water terrain source`
   — separate Rusanov source component、Williamson 5 initializer、compatible rejection。
7. `P6.07 test: validate orographic responses`
   — Williamson 5 reference と bounded linear hydrostatic mountain-wave gate。
8. `P6.08 diagnostics: close terrain torque and pressure-work budgets`
   — torque/work/omega diagnostics、CSV/metadata/restart regression。
9. `P6.09 test: add the JW06 smooth-topography benchmark`
   — steady/perturbed state、UMJS14 paired regression、reference envelope。
10. `P6.10 io: import one reproducible lat-lon terrain format`
    — strict CSV、fingerprint、bilinear interpolation、one bounded smoothing pass family。
11. `P6.11 test: gate the fixed-terrain C++ core`
    — independent refinement、rotated terrain、budgets、toolchains、validation report。

```text
P6.01 -> P6.02 -> P6.03 -> P6.04 -> P6.05
                                      |
                                      v
P6.06 -> P6.07 -> P6.08 -> P6.09 -> P6.10 -> P6.11
```

## 11. 完了チェックリスト

### Entry

- [ ] Phase 5 quantitative convergence と UMJS14 envelope gap が閉じている。
- [ ] ADR 0006 と実装の state、reconstruction、SSP-RK3、CFL 契約が一致している。

### Core

- [ ] terrain は一つの immutable field で、state/checkpoint に重複保存されない。
- [ ] flat case が Phase 5 と exact に一致し、uniform offset が dynamics を変えない。
- [ ] sloping hybrid surface の pressure-gradient error が refinement で減少する。
- [ ] top/bottom relative mass flux、layer positivity、surface-pressure bounds が維持される。
- [ ] shallow-water terrain source が別 component で、Williamson 5 reference を満たす。

### Validation and data

- [ ] DCMIP 2-0-0、linear hydrostatic response、JW06 の記録済み gate を満たす。
- [ ] mountain torque、terrain work、energy/AAM residual の発生項を説明できる。
- [ ] terrain centre を seam/corner へ移しても主要 norms/budgets が同等である。
- [ ] strict CSV interpolation と smoothing が fingerprint/mean/extrema property を満たす。
- [ ] restart、旧 phases、sanitizer/toolchain gate が通る。
- [ ] `docs/validation/phase-6.md` が command、table、cost、known gap を含む。
- [ ] `web/`、FrameV2、gateway、control protocol が変更されていない。

## 12. 参照

- P. A. Ullrich et al. (2012),
  [DCMIP 2012 Test Case Document v1.7](https://public.websites.umich.edu/~cjablono/DCMIP-2012_TestCaseDocument_v1.7.pdf),
  especially test 2-0-0 for pressure-gradient error over orography.
- D. L. Williamson et al. (1992),
  [A standard test set for numerical approximations to the shallow water equations in spherical geometry](https://doi.org/10.1016/S0021-9991(05)80016-6),
  test 5 for zonal flow over an isolated mountain.
- C. Jablonowski and D. L. Williamson (2006),
  [A baroclinic instability test case for atmospheric model dynamical cores](https://doi.org/10.1256/qj.06.12),
  the smooth-surface-geopotential dry benchmark paired with the Phase 5 flat UMJS14 case.

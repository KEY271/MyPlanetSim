# Phase 5 実装計画 — 地形なし乾燥 3D 静水圧力学コア

**状態:** 実装済み（既知 gap は [Phase 5 検証報告](validation/phase-5.md) を参照）。Phase 4 の検証済み C++ column core を開始点とし、3D state と
水平・鉛直結合の判断は
[ADR 0006](adr/0006-dry-hydrostatic-state-and-coupling.md) に固定する。

## 1. 到達点と実装順序

Phase 5 では cubed-sphere と hybrid sigma-pressure column を直積し、平坦で不透過な下端を持つ
乾燥・静水圧 primitive-equation core を C++ で実装する。先に数値 state、flux、時間積分、
diagnostics、checkpoint、標準試験を C++ 単独で完了させる。その gate を通過した同じ C++ state
から表示 frame を生成し、最後に既存 visualizer を model-level map と選択 column profile へ
対応させる。

```text
Milestone A: C++ dry 3D core
  state -> flux/source -> horizontal/vertical coupling -> driver
  -> diagnostics/I/O -> transport/dynamics benchmarks -> C++ validation gate

Milestone B: visualizer update
  C++ FrameV2 -> protocol decode -> gateway -> level map/profile
  -> live end-to-end gate
```

Milestone A の全 gate が通るまで `web/`、browser-facing protocol、gateway を変更しない。

成果物は次に限定する。

- `6*N*N*K` の cell-column state と version 付き checkpoint
- surface pressure、mass-weighted horizontal momentum、potential temperature、passive tracer
- Phase 2 unique-edge Rusanov flux と Phase 4 vertical recurrence の結合
- hybrid 面上の hydrostatic geopotential と horizontal pressure-gradient acceleration
- SSP-RK3、水平・鉛直 CFL、stage ごとの invariant validation
- 3D prescribed-flow transport、静止大気、外部モード、平坦地表 baroclinic wave の C++ gate
- dry mass、theta/tracer mass、energy、axial angular momentum、continuity の診断
- C++ gate 後の `FrameV2`、gateway preset、level map、選択 column profile

## 2. スコープ境界

### 2.1 含めるもの

- 固定 `A/B`、正の固定 model-top pressure、Lorenz full/half level
- 平坦下端、`surface_geopotential_m2_s2 = 0`、上下端 mass flux 0
- 乾燥理想気体、hydrostatic balance、adiabatic potential-temperature transport
- 一つの非負 passive tracer
- collocated cell-centred horizontal wind の Cartesian tangent representation
- Rusanov + linear reconstruction + Barth--Jespersen limiter の一つの production path
- Phase 4 donor-cell/linear-minmod の鉛直 transport path
- explicit SSP-RK3 と既存の任意な水平 diffusion
- CLI、diagnostics/state CSV、checkpoint/restart
- 小規模 local run 用の full-frame binary と既存 3D/2D renderer の再利用

### 2.2 含めないもの

- 地形、水平に変化する surface geopotential、山岳 torque（Phase 6）
- moist species、condensation、radiation、boundary layer、surface exchange（Phase 7 以降）
- vertical momentum、nonhydrostatic pressure、acoustic mode
- time-varying `A/B`、vertical remap、adaptive level、isentropic coordinate
- split-explicit、semi-implicit、IMEX、multirate、barotropic subcycling
- multi-tracer registry、vertical diffusion、sponge、new filter family
- MPI/OpenMP/GPU、parallel I/O、NetCDF、compression
- browser からの temperature/tracer/wind edit、任意 config/path 指定
- server-side field slicing、remote deployment、3D volume renderer

explicit reference が小規模の C++ gate を満たせない場合、時間積分器をその場で増やさず、失敗した
case、CFL、cost を記録して別 ADR で次の方式を選ぶ。baroclinic case を通すためだけの未診断 filter
や post-step clip は追加しない。

## 3. 3D grid、index、state

### 3.1 index と prognostic state

水平 cell order は ADR 0001 の `panel-major, then j, then i` を維持する。鉛直 index `k` は
`0..K-1` で下向きに増え、volume field は `offset(cell,k)=cell*K+k` とする。column が連続するため
Phase 4 operator を `span` で呼べる。horizontal loop のためだけに field を転置して二重保持せず、
汎用 N 次元 array framework も作らない。

```text
DryHydrostaticState
  time_s
  step
  surface_pressure_pa[C]
  horizontal_momentum_mass_kg_m_s[C*K]  # Vec3, M*u
  potential_temperature_mass_k_kg_m2[C*K]
  tracer_mass_kg_m2[C*K]
```

各 cell の `ps` から既存 `AtmosphericHybridCoordinate` を用いて `p_half`、`p_full`、`delta_p`、
`M=delta_p/g`、Exner、`u`、`theta`、`q`、`T`、density、`Phi` を stage ごとに診断する。
pressure/temperature/geopotential を state や checkpoint へ二重保存しない。性能上必要なら RHS 呼出し中
だけ生存する `DerivedColumn` buffer を一つ持つ。

### 3.2 validation と checkpoint

- vector/scalar の shape が `C*K`、`ps` が `C` と一致する。
- 全値が finite、`ps` が configured range 内、全 `delta_p` と `M` が正である。
- `theta>0`、`T>=temperature_floor_k`、nonnegative case の `q>=0`。
- momentum の radial component は scale-aware rounding tolerance 内である。
- `sum_k M=(ps-p_top)/g` が各 cell と全球で閉じる。

checkpoint layout は `dry_hydrostatic_cell_column_v1` とする。restart は fingerprint、`N`、`K`、
`A/B`、field count、state size の不一致を拒否する。既存 layout は変更しない。

## 4. 設定 schema

`experiment.kind = dry_hydrostatic` を追加する。Phase 4 の `vertical.*` key は座標、temperature floor、
鉛直 transport scheme/limiter/CFL に再利用し、同じ `A/B` schema を複製しない。Phase 5 固有 key は
次に限定する。

```text
dry_hydrostatic.test_case = isothermal_rest
dry_hydrostatic.reconstruction = linear
dry_hydrostatic.limiter = barth_jespersen
dry_hydrostatic.cfl = 0.45
dry_hydrostatic.diffusion_kind = none
dry_hydrostatic.diffusion_coefficient = 0
```

test-case 固有の標準定数は named case 内に置き、意味の薄い knob をすべて config 化しない。
`A/B`、解像度、dt、diffusion のように再現性・安定性へ直接効く値は canonical config と fingerprint
へ含める。

初期 preset は次とする。

- `phase5_isothermal_rest.cfg`
- `phase5_solid_body_transport.cfg`
- `phase5_dcmip_deformational.cfg`
- `phase5_dcmip_hadley.cfg`
- `phase5_linear_wave.cfg`
- `phase5_umjs14_steady.cfg`
- `phase5_umjs14_baroclinic.cfg`
- `phase5_visualizer_rest_n4.cfg`（Milestone B で追加）

## 5. 一つの RHS での結合順

### 5.1 derived column と horizontal flux

各 SSP-RK3 stage の state から全 cell の coordinate、primitive、thermodynamic、hydrostatic field を
再生成する。`ps` が変わった stage で前 stage の pressure/mass を再利用しない。hydrostatic integration
は既存実装を使い、surface geopotential は全 cell 0 とする。

各 level で Phase 2 と同じ unique edge を一度だけ走査する。left/right の `M,u,theta,q` を
linear + Barth--Jespersen で face へ再構築し、共有 edge tangent plane で Rusanov flux を作る。

```text
U = (M, M*u, M*theta, M*q)
F(U) = (M*u_n, M*u_n*u, M*u_n*theta, M*u_n*q)
F_R = 0.5*(F_L+F_R) - 0.5*a*(U_R-U_L)
a = max(|u_n| + sqrt(gamma*Rd*T)) over left/right
gamma = cp/(cp-Rd)
```

edge length を掛け、left から引き right へ加え、cell area で割る。air mass は全球で相殺し、constant
`theta/q` の scalar flux は同じ定数と mass flux の積になることを property test にする。

### 5.2 pressure gradient と Coriolis

full level ごとに

```text
specific_volume = Rd*T/p
G = grad_eta(Phi_full) + specific_volume*grad_eta(p_full)
d(M*u)/dt_pressure = -M*G
d(M*u)/dt_coriolis = -2*M*P_tangent(Omega cross u)
```

とする。両 gradient は同じ `least_squares_gradient` を用いる。uniform rest で roundoff まで 0、
manufactured hydrostatic field で水平収束することを full driver より先に確認する。pressure gradient、
Coriolis、advection、vertical transport、diffusion は別 RHS component に保持する。

### 5.3 vertical recurrence と transport

水平 air-mass tendency `H[c,k]` を column ごとに Phase 4 の
`diagnose_vertical_mass_flux` へ渡す。得た `ps_dot` と `F` を使い、momentum の 3 Cartesian component、
theta、tracer を同じ vertical scheme で更新する。

```text
X_dot[k] = H_X[k] + F_X[k] - F_X[k+1]
```

pressure/Coriolis/diffusion は mass transport ではないため `H` recurrence に混ぜず、momentum source
として合算する。全 column の `F[0]`、`F[K]`、continuity residual を残し、bottom residual を 0 へ
上書きしない。model-coordinate の pressure velocity は `omega_relative=g*F` として診断し、幾何学的
vertical velocity `w` と誤表示しない。

### 5.4 step size と stage update

```text
dt = min(requested_dt,
         horizontal_rusanov_dt,
         vertical_transport_dt,
         ps_bound_dt,
         end_time - time)
```

SSP-RK3 の各 stage で RHS、coordinate、CFL、positivity を再評価する。stage 更新後に momentum を
cell-centre tangent plane へ射影して invariant を検査する。負の layer/temperature/tracer が生じた場合、
cell、level、stage、flux、dt を含めて停止する。

## 6. benchmark

### 6.1 3D prescribed-flow transport

力学を有効にする前に、同じ horizontal/vertical scalar coupling を prescribed velocity で検証する。
DCMIP 2012 Test 1-1 と 1-2 の velocity、coordinate、smooth tracer を使うが、Phase 5 は一 tracer の
return error と mass/shape preservation に限定し、二 tracer correlation を含む完全 conformance は
主張しない。

- `solid_body`: Phase 1 の球面回転を全層へ拡張し、鉛直 flux 0 の水平回帰を分離する。
- `dcmip_deformational`: 水平 deformation と時間変化する鉛直運動、period 後の return error。
- `dcmip_hadley`: Hadley-like overturning と time reversal による水平・鉛直 coupling error。

`N` と `K` を別々に細分化する。constant tracer は machine scale、smooth tracer は期待次数、
nonnegative tracer は新しい extrema と負値なしを gate にする。

### 6.2 静止大気と linear external mode

uniform `ps`、isothermal または stable stratification、`u=0` の水平一様 column を結合済み driver で
複数 step 進め、全 RHS component と state が roundoff scale で不変であることを確認する。

次に hydrostatic rest へ小振幅で smooth な surface-pressure/velocity perturbation を与える。振幅を
半減したときの線形性、time/space refinement、位相速度、dry mass、seam/corner error を確認し、
Rusanov wave-speed estimate と horizontal CFL を独立に検証する。

### 6.3 flat-surface baroclinic case

Phase 5 は constant surface pressure と zero surface geopotential を持つ
Ullrich--Melvin--Jablonowski--Staniforth (2014; UMJS14) を採用する。smooth surface geopotential を
含む JW06 は Phase 6 へ残す。

1. unperturbed state を進め、`ps,u,T` error、energy/AAM drift、grid imprint、seam/corner error が
   水平・鉛直 refinement で減少することを確認する。
2. standard perturbation から、公開 diagnostic time の surface pressure、低層 temperature、
   vorticity の extrema・位置・area-weighted norms を高解像度 reference/envelope と比較する。

気候統計、spectrum package、ensemble は追加せず、deterministic initial-value problem で完了とする。

## 7. diagnostics と native I/O

各 sample に少なくとも次を出す。

- global dry mass、theta mass、tracer mass と source-adjusted residual
- total dry energy と absolute axial angular momentum の値、drift、component tendency
- `ps`、`delta_p`、temperature、theta、tracer、wind speed の min/max
- top/bottom mass flux、continuity residual、最大 horizontal/vertical CFL
- hydrostatic residual、pressure-gradient cancellation/error
- pre/post projection radial momentum
- panel 別および seam/corner 近傍の error norm
- non-finite/negative/floor violation count

energy/AAM は厳密保存を主張しない。drift を pressure work、advection、Coriolis、diffusion、time
integration に分け、refinement で減るか明示した diffusion と整合することを gate にする。

human mode は `dry_hydrostatic.csv`、`dry_hydrostatic_diagnostics.csv`、必要時の checkpoint を出す。
state CSV は `panel,i,j,k` と unit suffix 付き列を持つ。大規模 format は導入しない。restart は連続 run
と final state、step/time、budget、diagnostic schedule が一致することを検査する。

## 8. Milestone A — C++ 検証 gate

### 8.1 unit/property

- state index/shape/flatten/checkpoint round-trip と invalid input rejection
- 各 cell の hybrid pressure、mass sum、hydrostatic integration が Phase 4 と一致
- vector reconstruction/flux の seam crossing、edge antisymmetry、tangency
- horizontal global mass cancellation と constant theta/tracer consistency
- arbitrary horizontal `H` に対する全 column の Phase 4 recurrence
- uniform rest の全 RHS component が roundoff scale
- pressure-gradient manufactured field の水平収束
- horizontal/vertical CFL と `ps` bound が危険な step を拒否

### 8.2 convergence/regression/toolchain

- DCMIP 1-1/1-2 subset の horizontal/vertical independent refinement
- linear wave の amplitude linearity、space/time convergence、seam/corner error
- UMJS14 steady state の refinement と rotated-grid sensitivity
- UMJS14 perturbed wave の公開 time diagnostic envelope
- dry/theta/tracer mass、energy/AAM、continuity、hydrostatic budget
- restart/no-restart final state と diagnostics の一致
- `N`, `K`, step count に対する cell-layer-updates/s と peak memory の report
- 全 Phase 0--4 tests、dev/release、GCC/Clang、ASan/UBSan、format-check

tolerance は解析解、roundoff scale、observed order、または公開 reference envelope の出典を test/report
に記録する。`docs/validation/phase-5.md` の C++ 節に再現 command、table、budget、known gap を記録し、
ADR 0006 を accepted にしてから Milestone B を開始する。

## 9. Milestone B — C++ `FrameV2`

`FrameV2` は checkpoint ではなく表示用 snapshot である。C++ が validated state から primitive/
diagnostic field を生成し、browser は pressure/Exner/hydrostatic equation を再実装しない。既存 shallow-
water `FrameV1` と magic を維持し、decoder は magic/schema で分岐する。

```text
FrameV2 metadata
  schema_version = 2
  model_kind = dry_hydrostatic
  mapping, N, K, flatten_order, time, step, config_fingerprint

payload (little-endian Float64)
  surface_pressure_pa[C]
  pressure_pa[C*K]
  potential_temperature_k[C*K]
  temperature_k[C*K]
  wind_x_m_s[C*K]
  wind_y_m_s[C*K]
  wind_z_m_s[C*K]
  tracer_mixing_ratio[C*K]
```

volume order は solver と同じ `cell-major, then k` とする。wind speed は UI の algebraic derived field
としてよい。geopotential、vertical mass flux、tendency component は初期 UI に不要なので含めない。

full frame を採用し、field/level range API、compression、streaming decoder を先回りして作らない。
interactive allowlist は原則 `N<=24`, `K<=30` とし、gateway と C++ は frame count に加え run 全体の
published binary bytes を既定 `256 MiB` 以下にする。fixed payload は
`8*C*(1+7*K)` bytes（`N=24,K=30` で 5,833,728 bytes/frame）なので、request/runtime の両方で
overflow-safe に計算する。request の planned estimate を preflight し、CFL により step 数が増える場合は
C++ と gateway が累積実 byte 数を再検査する。実測で不足した場合だけ partial retrieval を別設計する。

C++ test は exact byte length、field offset、round-trip、truncated/wrong version/shape、non-finite、
fingerprint mismatch、atomic publication を検査する。frame 0 は初期化後・積分前の state と一致する。

## 10. protocol と gateway

既存 `RunRequestV1` の run/grid override と `/api/v1` lifecycle を維持する。Phase 5 UI は fixed preset
を選び `initialCondition.edits=[]` とする。新しい edit kind は追加しない。

capabilities は既存 field を壊さず、次の preset descriptor を additive に返す。

```text
presetDetails[]
  id
  modelKind = shallow_water | dry_hydrostatic
  frameSchemaVersion = 1 | 2
  levels = null | K
  supportedEdits[]
  maximumCellsPerPanel
```

gateway は preset ごとの edit/N limit を検査し、dry-hydrostatic preset への `gaussian_depth` を拒否する。
`frame.ready` へ `frameSchemaVersion` を加え、gateway は C++ `--describe-control` と preset descriptor の
compatibility を起動時に検査する。既存 loopback/token、Origin/Host、shell false、path containment、
cancel、atomic publish、bundle を維持する。累積 byte 超過は graceful cancel と明示 error にする。

## 11. visualizer

`VisualDatasetV1` と shallow-water UI を維持し `VisualFrameV2` を追加する。既存 `Globe3D` と `Map2D`
へ `levelSlice(frame,fieldId,k)` が返す cell scalar dataset を渡し、geometry、projection、grid overlay、
color scale を書き直さない。current field、level、time、selected cell は一つの state で共有する。

dry-hydrostatic UI の追加は次だけとする。

- model kind/preset selector
- model-level selector
- surface pressure、pressure、theta、temperature、tracer、wind speed の field selector
- 3D globe と 2D map の同じ level slice
- click した cell の、現在 field 対 pressure の column profile
- selected `panel,i,j`、level、pressure、raw C++ value、time/step の inspector

profile の pressure axis は downward-increasing log scale とし、V2 の `pressure_pa` を使う。新しい chart
dependency を入れず bounded SVG/canvas 一つで描く。surface field では level/profile を明示的に無効化し、
架空の vertical copy を作らない。Phase 3 の Gaussian editor は shallow-water preset でのみ表示し、
dry 3D mode の click は selection とする。未受信 time/level の補間はしない。

V2 は `256 MiB` run budget 内で全受信 frame を保持し、geometry は `N` ごとに cache、frame/level change
は color buffer 更新だけにする。測定なしに worker pool、virtual store、GPU volume rendering を追加しない。

## 12. Milestone B 検証 gate

- FrameV1/V2 magic dispatch、shape、offset、unit、flatten order
- V2 level slice と selected column が同じ C++ fixture の `(cell,k)` に一致
- invalid `N/K`、truncation、extra bytes、non-finite、fingerprint mismatch rejection
- preset descriptor、per-preset N/edit restriction、frame/byte budget
- gateway cancel/reconnect/bundle/path/atomic publication の V1/V2 回帰
- field/level/time/selection が globe、map、profile で同期
- top/bottom level、seam/corner cell、surface/constant/signed field の表示
- keyboard、focus、reduced motion、narrow viewport、WebGL failure の既存 gate
- frame/level change で geometry/node が再生成されない構造 test

live E2E は小さい `N=4`, `K<=16` rest preset で capabilities、empty-edit submit、frame 0/後続 frame、
level/map/profile の値一致、complete を確認する。続けて shallow-water V1 edit run と V2 cancel を通し、
orphan/partial frame がないことを確認する。`docs/validation/phase-5.md` は C++ gate と visualizer gate を
別見出しにし、後者だけの成功を Phase 5 完了としない。

## 13. コミット列

Milestone A の commit は C++ 数値責務と対応 test を含み `web/` を変更しない。Milestone B は C++
FrameV2、protocol、gateway、UI の順に進める。generic field/tracer/chart framework の空 shell は作らない。

1. `P5.01 docs: define the dry hydrostatic coupling contract`
   — ADR 0006 で state、mass flux、pressure gradient、time integration、flat boundary を固定。
2. `P5.02 config: add dry hydrostatic experiments`
   — experiment kind、既存 `vertical.*` 再利用、strict config tests。
3. `P5.03 dynamics: add the dry hydrostatic 3D state`
   — shape/index、primitive diagnosis、validation、flatten tests。
4. `P5.04 dynamics: transport layered mass and scalars on the cubed sphere`
   — reconstruction、Rusanov mass/momentum/theta/tracer flux、property tests。
5. `P5.05 dynamics: assemble dry hydrostatic pressure and rotation sources`
   — hydrostatic field、pressure force、Coriolis、uniform/manufactured tests。
6. `P5.06 dynamics: couple layered convergence to vertical mass flux`
   — Phase 4 recurrence、shared vertical transport、continuity tests。
7. `P5.07 dynamics: integrate the dry hydrostatic core`
   — component RHS、combined CFL、SSP-RK3、observer/cancel/restart hooks。
8. `P5.08 test: validate coupled 3D tracer transport`
   — DCMIP 1-1/1-2 one-tracer subset と independent refinement。
9. `P5.09 test: add flat-surface dry dynamics benchmarks`
   — rest、linear mode、UMJS14 steady/perturbed cases。
10. `P5.10 io: close dry hydrostatic budgets and restart`
    — diagnostics、CSV、checkpoint、CLI、presets、file regression。
11. `P5.11 test: gate the dry hydrostatic C++ core`
    — convergence、budgets、toolchains、cost、validation report。ここまで `web/` を変更しない。
12. `P5.12 io: publish validated dry hydrostatic frames`
    — FrameV2 writer/reader、observer、size/atomic tests。
13. `P5.13 web: decode dry hydrostatic frames and level slices`
    — VisualFrameV2、magic dispatch、level/column adapter、V1 regression。
14. `P5.14 gateway: serve bounded dry hydrostatic runs`
    — preset descriptor/restriction、byte budget、V2 lifecycle。
15. `P5.15 web: visualize model levels and selected columns`
    — controls、既存 globe/map slice、shared selection、pressure profile。
16. `P5.16 test: validate the Phase 5 native visualization path`
    — V2 E2E、V1 regression、cancel/memory/renderer gate と report 完成。

```text
P5.01 -> ... -> P5.11  [Milestone A / C++ gate]
                    |
                    v
P5.12 -> P5.13 -> P5.14 -> P5.15 -> P5.16  [Milestone B]
```

## 14. 完了チェックリスト

### C++ core

- [x] state/mass/index/staggering が ADR 0005/0006 と一致する。
- [x] horizontal unique-edge flux と vertical recurrence が一つの dry-mass budget を閉じる。
- [x] constant theta/tracer、positivity、tangency、上下端 flux が全 stage で成立する。
- [ ] pressure gradient、linear wave、3D transport が設計どおり収束する。
      （bounded initialization/regression まで。定量的な refinement は既知 gap。）
- [ ] flat rest と UMJS14 steady/perturbed case が記録した gate を満たす。
      （rest は満たす。UMJS14 の reference envelope 比較は既知 gap。）
- [x] mass/theta/tracer、energy/AAM、continuity/hydrostatic residual を説明できる。
- [x] checkpoint/restart、CLI/CSV、旧 test、sanitizer/toolchain gate が通る。
- [x] C++ validation report 完了まで `web/` が変更されていない。

### Visualizer

- [x] FrameV2 が C++ state から生成され、V1 と独立に version/shape/size 検証される。
- [x] gateway が preset/edit/N/frame/byte budget と既存 local security を守る。
- [x] globe/map/profile の field/level/time/cell が同じ C++ frame に一致する。
- [x] browser が pressure/hydrostatic/transport を再計算しない。
- [x] shallow-water V1 UI と live edit が回帰しない。
- [x] native V2 E2E と cancel が再現可能である。

## 15. 参照

- A. J. Simmons and D. M. Burridge (1981),
  [An Energy and Angular-Momentum Conserving Vertical Finite-Difference Scheme and Hybrid Vertical Coordinates](https://doi.org/10.1175/1520-0493(1981)109%3C0758:AEAAMC%3E2.0.CO;2).
- DCMIP (2012),
  [Dynamical Core Model Intercomparison Project test case document](https://public.websites.umich.edu/~cjablono/DCMIP-2012_TestCaseDocument_v1.7.pdf).
- P. A. Ullrich, T. Melvin, C. Jablonowski, and A. Staniforth (2014),
  [A proposed baroclinic wave test case for deep- and shallow-atmosphere dynamical cores](https://doi.org/10.1002/qj.2241).
- C. Jablonowski and D. L. Williamson (2006),
  [A baroclinic instability test case for atmospheric model dynamical cores](https://doi.org/10.1256/qj.06.12),
  retained for Phase 6 because its balanced state uses smooth surface geopotential.

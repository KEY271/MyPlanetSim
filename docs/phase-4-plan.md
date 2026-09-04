# Phase 4 実装計画 — 鉛直 1D と hybrid sigma-pressure

**状態:** 実装・検証済み。Phase 3 の完了状態を開始点とし、鉛直座標と column state の
判断は [ADR 0005](adr/0005-hybrid-vertical-coordinate-and-column-state.md) に固定する。

## 1. 到達点

Phase 4 では、Phase 5 の乾燥 3D 静水圧力学へ組み込める鉛直 1D column core を C++ で
実装し、水平格子から独立に座標、熱力学、静水圧、鉛直質量 flux、保存型 scalar 輸送を
検証する。Phase 4 は C++ の数値仕様・出力・検証 gate までとし、cubed-sphere と結合した
user-facing visualizer は 3D state が確定する Phase 5 で更新する。

成果物は次を含む。

- version 付き `A/B` hybrid coordinate と full/half-level geometry
- `delta_p/g`、Exner、温位・温度変換、静水圧 geopotential
- surface-pressure tendency と水平質量 tendency からの鉛直 mass flux 診断
- mass-weighted potential temperature と passive tracer の保存型鉛直輸送
- donor-cell 基準方式と bounded linear reconstruction
- isothermal、dry-adiabatic、moving-`ps`、manufactured transport の column case
- column diagnostics、checkpoint、CSV、CLI preset
- Phase 4 の再現手順と検証報告

Phase 4 の完了は C++ 単独の unit、convergence、conservation gate で判定する。column 専用の
Web tab、binary frame、live protocol は作らない。検証用 profile は CSV と数値 test で確認し、
Phase 5 で一体化した globe/map/column profile UI を設計する。

## 2. 実装順序の境界

### 2.1 必須の順序

Phase 4 は C++ column core の一つの milestone とする。

```text
Phase 4 C++ column core
  coordinate -> thermodynamics -> hydrostatic -> state -> mass flux
  -> transport -> driver -> diagnostics/I/O -> C++ validation
```

Phase 4 では `web/` を変更しない。これにより Web schema や表示都合が C++ state、staggering、
単位、数値方式を決めることを防ぐ。Phase 5 では結合済み C++ 3D state を authoritative data とし、
ブラウザ内に別の静水圧積分や鉛直輸送を実装しない。

### 2.2 Phase 4 に含めないもの

- cubed-sphere と鉛直 column の直積 field、水平 flux との実結合（Phase 5）
- 水平風、Coriolis、圧力傾度、運動エネルギー、全エネルギー closure（Phase 5）
- 地形、傾いた coordinate surface、surface geopotential の水平変化（Phase 6）
- moist thermodynamics、condensate、放射、境界層、対流調節
- time-varying `A/B`、鉛直 AMR、isentropic coordinate、conservative remap
- semi-implicit 鉛直音波/重力波 solver
- Web protocol、binary live frame、gateway、visualizer の変更（Phase 5）

Phase 4 の column driver へ入る水平 tendency は manufactured forcing であり、水平力学の代用品
ではない。総エネルギーが閉じると主張せず、dry mass、potential-temperature mass、tracer mass、
hydrostatic residual を Phase 4 の budget とする。

## 3. 鉛直配置と state

### 3.1 index と pressure

`docs/conventions.md` のとおり `k` は下向きに増える。

```text
half level 0          p_top
full level 0
half level 1
...
full level nz-1
half level nz         ps
```

固定係数に対し

```text
p_half[k] = A_pa[k] + B[k] * ps
delta_p[k] = p_half[k+1] - p_half[k]
air_mass[k] = delta_p[k] / g
```

とする。production column は `A[0]=p_top`, `B[0]=0`, `A[nz]=0`, `B[nz]=1` を必須とし、
`sum(air_mass)=(ps-p_top)/g` が telescoping により成立する。

### 3.2 full-level pressure

`Pi=(p/p0)^kappa` を用い、full-level Exner を log-pressure 区間平均で定義する。

```text
Pi_full[k] = (Pi_half[k+1] - Pi_half[k])
             / (kappa * log(p_half[k+1] / p_half[k]))
p_full[k] = p0 * Pi_full[k]^(1/kappa)
```

近接圧力では `log1p`/`expm1` 相当の安定式または解析的極限を用いる。全 full level が対応する
half level の厳密な内側にあり、有限・正であることを検証する。算術平均圧力を暗黙に混用しない。

### 3.3 `VerticalColumnState`

最小 state は次とする。

```text
VerticalColumnStateV1
  time_s
  step
  surface_pressure_pa
  potential_temperature_mass[nz]  # air_mass * theta
  tracer_mass[nz]                 # air_mass * q
```

`air_mass`、`theta`、`q`、pressure、temperature、geopotential は coordinate と state から診断する。
同じ量を prognostic と diagnostic の両方に保存しない。全値は finite、`ps` は設定範囲内、
`air_mass>0`、`theta>=temperature_floor_k/Pi_full`、非負 tracer case では `tracer_mass>=0` を
stage 境界で検査する。

Checkpoint layout ID は `vertical_column_hybrid_v1` とし、既存 ODE、transport、shallow-water
layout を変更しない。config fingerprint に `A/B` 全要素、planet thermodynamic constants、floor、
scheme を含める。

## 4. 座標係数と設定 schema

### 4.1 係数 validation

低レベルの `HybridPressureCoefficients` は pure-pressure limit も表現できるよう、長さ、有限性、
指定 `ps` 範囲での単調性と layer thickness を検査する。その上に production 用
`AtmosphericHybridCoordinate` factory を置き、次を拒否する。

- `nz < 1`、`A/B` の長さ不一致、非有限値
- `p_top<=0`、endpoint 規約違反
- `B` が `[0,1]` 外、または下向きに減少する profile
- configured `ps_min` または `ps_max` で非単調な interface pressure
- いずれかの layer が `minimum_pressure_thickness_pa` 未満
- `ps_min<=p_top`、`ps_reference` が許容範囲外

`A` 単独の単調性は要求しない。実際に意味を持つ `A+B*ps` を検査する。pressure は `ps` の
affine function なので、許容区間の両端を検査すれば全区間を cover できる。

### 4.2 config

`experiment.kind = vertical_column` を追加する。初期 schema は次とする。

```text
experiment.kind = vertical_column
vertical.test_case = isothermal | dry_adiabatic | moving_surface_pressure | manufactured_transport
vertical.levels = <positive integer>
vertical.a_half_pa = <comma-separated nz+1 values>
vertical.b_half = <comma-separated nz+1 values>
vertical.surface_pressure_pa = <Pa>
vertical.minimum_surface_pressure_pa = <Pa>
vertical.maximum_surface_pressure_pa = <Pa>
vertical.minimum_pressure_thickness_pa = <Pa>
vertical.surface_geopotential_m2_s2 = <m2 s-2>
vertical.initial_temperature_k = <K>
vertical.initial_potential_temperature_k = <K>
vertical.temperature_floor_k = <K>
vertical.transport_scheme = donor_cell | linear
vertical.limiter = none | minmod
vertical.cfl = <0,1]
vertical.forcing_amplitude = <case-specific finite value>
diagnostics.interval_steps = <positive integer>
```

list parser は空要素、末尾 comma、非有限値、unit suffix の混在を拒否し、canonical writer は
binary64 の round-trip に十分な桁数で出力する。case 固有の不要 key も拒否し、暗黙 default に
よる別 experiment の成立を避ける。

標準 preset は小さい `nz=8` smoke、収束用 `nz=16/32/64/128` を持つ。production hybrid profile
をコード内の名前だけで生成せず、各 preset に `A/B` 数値列を明記して再現可能にする。test helper
には pressure-limit と positive-top sigma-like profile の generator を置いてよい。

## 5. 乾燥熱力学と静水圧

### 5.1 pure functions

`dry_thermodynamics` は allocation しない scalar pure function とする。

```text
kappa = Rd / cp
Pi(p) = (p/p0)^kappa
T(theta,p) = theta*Pi(p)
theta(T,p) = T/Pi(p)
rho(p,T) = p/(Rd*T)
```

入力の pressure、temperature、`Rd`、`cp`、`p0` は正、`cp>Rd` を要求する。正常な正値入力で
overflow/underflow しない範囲を config validation で固定し、NaN を default 値へ置換しない。

### 5.2 hydrostatic integration

surface interface から上向きに

```text
Phi_half[nz] = Phi_surface
Phi_half[k] = Phi_half[k+1]
              + cp*theta[k]*(Pi_half[k+1]-Pi_half[k])
```

を積分する。full-level geopotential は同じ layer integral を lower interface から `Pi_full` まで
適用する。`height_m=Phi/g` は表示・診断用 derived quantity であり、geopotential と混同しない。

hydrostatic residual は再積分値との差ではなく、各 layer の離散 balance

```text
Phi_half[k]-Phi_half[k+1]
  + cp*theta[k]*(Pi_half[k]-Pi_half[k+1])
```

を独立に評価する。isothermal と constant-theta dry-adiabatic column は layer ごとの解析積分と
比較する。

## 6. surface pressure と鉛直 mass flux

### 6.1 continuity contract

Phase 5 が各 horizontal cell から渡す量を先に C++ API として固定する。

```text
H[k]     horizontal air-mass tendency, kg m^-2 s^-1
H_theta  horizontal theta-mass tendency, K kg m^-2 s^-1
H_q      horizontal tracer-mass tendency, kg m^-2 s^-1
```

`H>0` は net convergence とする。impermeable boundary の下で

```text
ps_dot = g*sum(H)
target_mass_dot[k] = delta_B[k]*ps_dot/g
F[0] = 0
F[k+1] = F[k] + H[k] - target_mass_dot[k]
```

により正下向きの coordinate-relative mass flux `F` を診断する。`F[nz]`、telescoping mass
residual、`ps_dot-g*sum(H)` を diagnostics に残す。許容差超過は failure とし、bottom flux を
ゼロへ上書きしない。

### 6.2 scalar flux

mass-weighted scalar `M*q` は

```text
d(M*q)[k]/dt = H_q[k] + F_q[k] - F_q[k+1]
```

で更新する。`F_q=F*q_face` とし、theta と tracer が同一の `F` と reconstruction policy を使う。
最初に donor-cell を実装し、その後 nonuniform pressure/mass coordinate 上の piecewise-linear
reconstruction と minmod limiter を追加する。

limiter は nonnegative tracer と positive temperature を保証する。更新後 clip は保存量を変えるため
禁止し、許容 CFL を満たしても positivity を失う場合は layer、interface、flux、dt を含む error で
停止する。

### 6.3 time step

各 layer の outward mass rate から

```text
outward_H[k] = max(-H[k], 0)
dt_out[k] = cfl * air_mass[k]
            / (max(-F[k],0) + max(F[k+1],0) + outward_H[k])
```

を評価し、有限な最小値と `run.time_step_s` の小さい方を stage 共通上限とする。manufactured
forcing が mass を除く場合は `ps_min` と layer mass floor までの距離も上限に含める。最終 step は
`end_time_s` へ揃える。SSP-RK3 の各 stage で coordinate、forcing、flux、positivity を再評価する。

## 7. column driver と benchmark

### 7.1 static hydrostatic cases

- `isothermal`: 指定 `T0` から各 full level の `theta=T0/Pi_full` を作る。
- `dry_adiabatic`: 指定 constant `theta0` から `T=theta0*Pi_full` を作る。

両 case は無 forcing で state を保持し、hydrostatic profile と analytic geopotential を出力する。
静水圧積分自体は state tendency を生成しないため、「静止柱が加速しない」は state bitwise 不変、
zero mass flux、zero budget drift として検査する。鉛直 momentum は Phase 4 に導入しない。

### 7.2 moving surface pressure

analytic な `H[k,t]` を与えて `ps` が許容範囲内で周期変化し、`H_theta=theta0*H`、
`H_q=q0*H` とする。診断された `F` を通して constant theta/tracer が丸め誤差 scale で保持され、
layer mass sum、上下 boundary、周期終了時 state が一致することを確認する。

### 7.3 manufactured transport

滑らかな非定数 scalar と解析的な `H/H_q` を組にした forced problem を使い、`nz` と `dt` を
別々に細分化する。donor-cell は一次、unlimited linear は滑らかな解で二次を期待し、minmod case
は extrema で boundedness を優先して次数 gate を分ける。

## 8. diagnostics と I/O

### 8.1 diagnostics

共通 sample に次を出す。

- total dry mass と `(ps-p_top)/g` の差
- total potential-temperature mass、total tracer mass と各 budget residual
- minimum/maximum `delta_p`、air mass、theta、temperature、tracer
- top/surface flux、continuity residual、maximum vertical CFL
- hydrostatic `L1/L2/Linf` residual
- state、pressure、diagnostic field の non-finite count

forced case は initial 値だけとの drift ではなく、time-integrated source/flux を差し引いた budget
residual を報告する。reduction order を固定し、既存 diagnostics CSV を変更せず `column_diagnostics.csv`
を追加する。

### 8.2 human output と checkpoint

`experiment.kind=vertical_column` を既存 CLI dispatch に追加する。human mode は
`column_profile.csv`、`column_diagnostics.csv`、checkpoint を出力し、Phase 0–3 の file name と列を
変更しない。profile CSV は level location を明示するため layer と interface を別 file または
`location,index` 列で区別し、pressure・geopotential・unit を header 名へ含める。

### 8.3 Phase 5 への出力 handoff

Phase 4 は human-readable な column profile/diagnostics CSV と checkpoint までを所有する。
column 専用 binary frame や live control schema は固定しない。Phase 5 で 3D state、field 数、
staggering、必要な部分取得方法、frame size 上限を確定してから、surface map、model level、選択 column
profile を一つの run/frame contract として設計する。

## 9. C++ 検証 gate

### 9.1 coordinate と thermodynamics

- pure-pressure fixture (`B=0`) が `p_half=A` と既知の `delta_p/g` に一致する。
- positive-top sigma-like fixture `A=p_top*(1-B)` が
  `p=p_top+B*(ps-p_top)` に一致する。
- production endpoints で mass sum が `(ps-p_top)/g` に丸め誤差 scale で一致する。
- `ps_min` と `ps_max` の間で全 interface が正・単調、full pressure が layer 内部にある。
- `theta -> T -> theta`、`p -> Pi` が representative/extreme valid range で round-trip する。
- invalid length、endpoint、nonmonotone、thin layer、zero top、NaN/Inf を constructor/config が拒否する。

### 9.2 hydrostatic

- isothermal column は `Delta Phi=Rd*T*log(p_lower/p_upper)` と一致する。
- constant-theta column は `Delta Phi=cp*theta*(Pi_lower-Pi_upper)` と一致する。
- arbitrary surface geopotential は全 profile を同じ定数だけ shift する。
- analytic smooth stratification は `nz=16,32,64,128` で最低二次収束を示す。
- pressure range が狭い layer でも cancellation による NaN/負 full pressure を生じない。

### 9.3 mass flux と transport

- arbitrary finite `H` に対し recurrence、top/bottom flux、mass tendency が machine scale で閉じる。
- moving-`ps` constant theta/tracer case の最大誤差が
  `256 epsilon * max(1, field scale)` を超えない。
- donor-cell が dry mass と scalar mass を保存し、nonnegative tracer を保つ。
- linear unlimited manufactured solution の `L1/L2` observed order が原則 `>=1.7`。
- minmod case は新 extrema を作らず、全 layer mass/temperature が floor より上にある。
- time-step refinement で SSP-RK3 の observed order が原則 `>=2.7`。
- restart 前後で final state と全 budget が許容誤差内で一致する。

### 9.4 regression と sanitizer

- 既存 Phase 0–3 C++ test を全て維持する。
- Debug、GCC/Clang、ASan/UBSan、format-check を通す。
- column CSV/checkpoint は round-trip し、wrong layout、wrong shape、non-finite、fingerprint mismatch を
  拒否する。

この節の gate と `configs/phase4_*` の clean run が全て通った時点を Phase 4 完了とする。

## 10. コミット規約

各コミットは一つの C++ 数値責務と対応 test を含む。Phase 4 では `web/`、browser-facing
protocol、gateway を変更しない。既存 Phase 3 の C++/Web test は回帰 gate として維持する。

generated column output、checkpoint、test report は commit しない。reference tolerance は
式・conditioning・truncation error から決め、失敗後に根拠なく緩和しない。

想定する新規 module 境界は次とする。実装中に名前を変更しても責務は混ぜない。

```text
include/myplanetsim/vertical/hybrid_pressure_coordinate.hpp
include/myplanetsim/vertical/hydrostatic_column.hpp
include/myplanetsim/vertical/vertical_transport.hpp
include/myplanetsim/thermodynamics/dry_thermodynamics.hpp
include/myplanetsim/dynamics/vertical_column_state.hpp
include/myplanetsim/dynamics/vertical_column_driver.hpp
include/myplanetsim/diagnostics/vertical_column_diagnostics.hpp
```

対応する `src/` 実装と `tests/test_*` を同じ commit へ置く。

## 11. コミット列

### P4.01 — vertical coordinate ADR

```text
docs: define the hybrid vertical column contract
```

ADR 0005、index、pressure、mass、state、mass-flux sign、remap 非採用を固定する。

### P4.02 — vertical experiment config

```text
config: add strict vertical column experiments
```

`vertical_column` kind、`A/B` list parser、case/scheme/floor/pressure range validation と preset を追加。

### P4.03 — hybrid coordinate geometry

```text
vertical: construct hybrid pressure columns
```

half/full pressure、`delta_p`、air mass、stable log-pressure mean と invalid profile tests を追加。

### P4.04 — dry thermodynamics

```text
thermodynamics: add dry Exner conversions
```

Exner、theta/temperature、density pure functions と round-trip/range tests を追加。

### P4.05 — hydrostatic column

```text
vertical: integrate hydrostatic geopotential
```

interface/full geopotential、height、residual、isothermal/adiabatic analytic tests を追加。

### P4.06 — column state and checkpoint

```text
dynamics: add a versioned vertical column state
```

mass-weighted state、invariant validation、flatten/unflatten、checkpoint layout を追加。

### P4.07 — vertical mass-flux diagnosis

```text
vertical: diagnose hybrid-coordinate mass flux
```

`ps_dot`、target layer mass tendency、interface recurrence、boundary/budget diagnostics を追加。

### P4.08 — conservative donor-cell transport

```text
vertical: transport column scalars conservatively
```

theta/tracer の共有 donor-cell flux、RHS、CFL、mass/constant/positivity tests を追加。

### P4.09 — bounded linear transport

```text
vertical: reconstruct bounded column profiles
```

nonuniform coordinate 上の linear reconstruction、minmod limiter、space convergence gate を追加。

### P4.10 — column forcing and driver

```text
dynamics: integrate forced vertical columns
```

SSP-RK3 stage coupling、static/moving-`ps`/manufactured case と restart を追加。

### P4.11 — column diagnostics

```text
diagnostics: close vertical column budgets
```

dry/scalar mass、source-integrated residual、hydrostatic residual、range/CFL samples を追加。

### P4.12 — column CLI and I/O

```text
io: publish vertical column results
```

CLI dispatch、profile/diagnostics CSV、checkpoint/restart と configs を追加する。既存 V1
shallow-water machine run と `FrameV1` は変更しない。

### P4.13 — C++ Phase 4 gate と検証報告

```text
test: gate the vertical column core
```

coordinate limits、analytic hydrostatic、space/time convergence、moving-`ps`、restart、sanitizer の
Phase 4 gate を固める。clean build/run、数値 tolerance、収束表、budget、既知 gap を
`docs/validation/phase-4.md` に記録する。

## 12. 依存関係

```text
P4.01 -> P4.02 -> P4.03 -> P4.04 -> P4.05 -> P4.06
                                                   |
                                                   v
P4.07 -> P4.08 -> P4.09 -> P4.10 -> P4.11 -> P4.12 -> P4.13
```

実装は commit 番号順を基本とする。P4.04 の scalar thermodynamics tests は P4.03 と並行可能だが、
履歴上は coordinate convention 確定後に置く。visualizer work は Phase 5 の 3D state 検証後に行う。

## 13. 完了チェックリスト

- [x] `A/B`、full/half level、pressure/mass convention が ADR と code で一致する。
- [x] arbitrary valid `ps` で pressure が単調、layer mass が正、mass sum が閉じる。
- [x] isothermal と dry-adiabatic hydrostatic profile が解析式に一致する。
- [x] moving-`ps` で上下 boundary flux、continuity、constant scalar が整合する。
- [x] donor-cell/linear transport の conservation、boundedness、設計収束が確認できる。
- [x] C++ checkpoint/restart と profile/diagnostics CSV が通る。
- [x] Phase 0–3 C++ gate、GCC/Clang、ASan/UBSan、format が回帰しない。
- [x] `web/`、Phase 3 protocol、shallow-water UI に Phase 4 固有変更が入っていない。
- [x] Phase 4 検証報告が clean checkout から再生成できる。

# Phase 8 実装計画 — 惑星設定、陸海 surface、軌道強制

**状態:** bounded implementation 完了。Phase 5--7 の production gate と ADR 0008 の acceptance は
未完了だが、2026-09-05 のユーザー判断により既知 gap を明記したまま予定コミット列を実装した。
これは production validation の完了を意味せず、結果と未完了項目は
[Phase 8 検証報告](validation/phase-8.md)に記録する。

## 1. 到達点と設計原則

Phase 8 は、乾燥・静水圧 core をコード変更なしで複数の地球型惑星条件へ切り替えられるようにする。
惑星定数を変えるだけでなく、恒星方向、解析的な惑星 benchmark forcing、固定地表、陸海の熱慣性を
一つの再現可能な run contract にする。

海洋は流体として解かない。海面下地形、海洋深度、海流、海洋熱輸送、混合層の鉛直構造は持たず、
陸と異なる面積熱容量を持つ、一セル一温度のゼロ次元 surface reservoir とする。海岸セルを二値化せず、
全水平セルに `0 <= f_land <= 1` を保持する。純海洋は `f_land=0`、純陸地は `f_land=1`、海岸を含む
mixed cell は必ず `0<f_land<1` のまま計算する。

```text
Phase 5--7 production gaps / ADR 0008
                   |
                   v
Milestone A: planet and orbit contract
  SI constants -> rotating-frame star vector -> dimensionless metadata
                   |
                   v
Milestone B: immutable surface boundary
  uniform / Earth geography -> Phi_s + f_land + provenance
                   |
                   v
Milestone C: two deliberately separate forcing paths
  analytic benchmark forcing | surface energy balance
                   |
                   v
Milestone D: presets and bounded validation
  Earth regression / Earth geography / slow / rapid / synchronous
```

解析的 Newtonian forcing は既知の大気 benchmark を再現するための経路、surface energy balance は
陸海熱容量と日変化を検証するための経路である。後者を前者の published reference と直接比較しない。

## 2. スコープ境界

### 2.1 含めるもの

- 既存の `R, Omega, g, Rd, cp, p0` と整合する軌道・恒星パラメータ
- 円・楕円軌道、軸傾斜、初期軌道位相から得る planet-fixed 恒星方向と `mu0`
- axisymmetric と substellar の二種類に限定した解析的 Newtonian forcing
- 時間不変な `surface_geopotential_m2_s2[C]` と `land_fraction[C]`
- `uniform` geography と、出典・生成物を固定した `earth` geography
- 一セル一つの surface temperature と面積熱容量だけを持つ陸海 slab
- surface--lowest-atmospheric-layer の sensible heat exchange と全球 energy budget
- Earth-like regression、Earth geography、slow rotator、rapid rotator、tidally locked preset
- Rossby 数、Burger 数、変形半径、回転・放射・surface 時間尺度比の metadata/diagnostics
- old config、old checkpoint、`physics.kind=none|held_suarez` の exact compatibility

### 2.2 含めないもの

- 海洋深度、bathymetry、海流、水平 ocean heat transport、鉛直混合、upwelling
- mixed cell 内の陸温度と海温度の分離。mixed cell も一つの `T_s` だけを持つ
- 陸海別 albedo、emissivity、roughness、drag、evaporation。Phase 8 で異なるのは熱容量だけ
- 動く海岸線、海面高度、洪水、湖・海氷・積雪・植生、地下熱拡散
- moist physics、雲、凝結、gray/line-by-line radiation、温室効果の conformance claim
- tides、precession、nutation、複数恒星、衛星、軌道の時間発展
- nonhydrostatic dynamics、top sponge、gravity-wave drag
- THAI conformance。THAI は moist/radiative physics を導入する Phase 9 以後へ残す
- arbitrary process registry、parameter-sweep DSL、online data download、GDAL/NetCDF runtime dependency
- FrameV3 と Web UI。Phase 8 の surface field は standalone CSV と checkpoint で先に検証する

## 3. planet、orbit、恒星方向 contract

### 3.1 設定

既存の `planet.*` は SI 値を明示する現在の契約を維持し、`planet.preset` のような値の二重管理は
追加しない。preset は version 管理された `.cfg` ファイルであり、展開後の全値が canonical config と
metadata に残る。

`substellar` forcing または surface energy balance を選んだときだけ、次を必須にする。
axisymmetric forcing は恒星位置を使わないため orbit key を要求しない。

```text
orbit.period_s
orbit.eccentricity
orbit.obliquity_rad
orbit.longitude_of_periapsis_rad
orbit.initial_mean_anomaly_rad
orbit.initial_substellar_longitude_rad
star.flux_at_semimajor_axis_w_m2
```

`period_s>0`、`0<=eccentricity<1`、有限な角度、正の恒星 flux を要求する。時刻は
`run.start_time_s` と同じ epoch からの SI seconds である。`planet.rotation_rate_rad_s` は慣性系に対する
符号付き sidereal angular velocity とし、既存の負回転許容を維持する。

### 3.2 幾何

平均運動 `n=2*pi/P_orbit` とし、平均近点角から Kepler 方程式を決定論的に解く。公転距離 `r` により

```text
S(t) = S_a * (a/r(t))^2
mu0(c,t) = max(0, dot(cell_center[c], star_direction_planet_fixed(t)))
```

を得る。`t=0`、zero obliquity、zero eccentricity では substellar longitude が
`orbit.initial_substellar_longitude_rad` と一致する。prograde の恒星日候補は
`2*pi/abs(Omega-n)` とし、`Omega=n` は有限値を捏造せず `synchronous=true`、stellar day undefined と
metadata に書く。retrograde も符号付き角速度から同じ座標変換で扱う。

角度の曖昧さを避けるため、`longitude_of_periapsis` は北半球の春分方向から prograde に測る。
true anomaly `nu` に対する planet-to-star の solar longitude を
`L_star=nu+longitude_of_periapsis+pi` とし、equatorial inertial vector を
`(cos(L_star), cos(obliquity)*sin(L_star), sin(obliquity)*sin(L_star))` とする。`t=0` の z 軸回転位相は
この vector の right ascension と `initial_substellar_longitude` から一意に決める。Kepler solve は
safeguarded Newton/bisection、固定 iteration 上限、scale-aware residual を使い、未収束を明示的に失敗
させる。

幾何 unit test は少なくとも次を固定する。

- `mu0` は `[0,1]`、night side は exact zero
- 球面積分した瞬時入射 power は `pi*R^2*S` に quadrature error 内で一致
- zero obliquity/eccentricity の非同期回転では substellar longitude が `n-Omega` で移動
- synchronous case では恒星方向が planet-fixed で不変
- eccentric orbit の periapsis/apoapsis flux 比が解析値と一致

## 4. 固定 surface geography contract

### 4.1 run-owned field

ADR 0008 の `SurfaceOrography` を包含する、不変な run-owned object を導入する。

```text
SurfaceBoundary
  surface_geopotential_m2_s2[C]
  land_fraction[C]
  source_id
  source_fingerprint

0 <= land_fraction[c] <= 1
C_surface[c] = f_land[c] * C_land
             + (1 - f_land[c]) * C_ocean
```

field 順序は既存どおり `panel-major, then j, then i` とする。`land_fraction` を閾値で 0/1 に丸めない。
`C_land` と `C_ocean` は正の面積熱容量 `[J m-2 K-1]` であり、密度、比熱、深さへ分解しない。
したがって config に ocean depth は存在しない。`C_ocean>C_land` は標準 preset の選択だが、共通 parser
は両方が正で有限であることだけを要求する。

surface energy balance では次を明示する。

```text
surface.land_heat_capacity_j_m2_k
surface.ocean_heat_capacity_j_m2_k
surface.initial_temperature_k
surface.albedo
surface.emissivity
surface.air_exchange_coefficient_w_m2_k
surface.internal_heat_flux_w_m2
```

Phase 8 の登録初期値はそれぞれ `2.0e6 J m-2 K-1`、`4.0e7 J m-2 K-1`、`288 K`、
`0.3`、`1.0`、`10 W m-2 K-1`、`0 W m-2` とする。値は物理層の深さを意味せず、まず熱慣性比 20 の
応答差を測るための idealized preset 値である。production 結果を見た後に暗黙に調整せず、変更時は
config と validation report に残す。

`surface.geography=uniform` は全セルで一つの `surface.uniform_land_fraction` を使い、既存
`orography.*` と組み合わせられる。`surface.geography=earth` は elevation と land fraction の両方を
Earth product から供給し、別の non-flat `orography.kind` との併用を拒否する。key が欠落した既存 config
は surface 無効として扱い、canonical text と fingerprint を変えない。

### 4.2 Earth geography

Earth は runtime download を行わない。repository に小さい派生 product、生成 script、manifest を置く。
raw source は次に固定する。

- elevation: NOAA NCEI ETOPO 2022 60 arc-second Ice Surface global relief
- land/water polygon: NOAA GSHHG 2.3.7

generator は raw file の version、取得元、取得日、SHA-256、license を manifest に記録する。water sample
の bathymetry は捨てて高さを exact zero とする。land polygon 内では ETOPO surface elevation を使うため、
海面下にある陸地の負 elevation は保持できる。ice sheet は Ice Surface elevation とする。

派生 product は 1 degree x 1 degree の rectilinear cell average に固定し、各 source cell は次を持つ。

```text
longitude_deg,latitude_deg,mean_surface_height_m,land_fraction
```

`land_fraction` は polygon と source-cell 面積の交差比、`mean_surface_height_m` は水面を 0 とした
source-cell 全面積平均である。target cubed-sphere cell へは固定次数の面積 quadrature で平均し、
quadrature order、source resolution、runtime FNV-1a fingerprint を canonical metadata に残す。
補間値ではなく target-cell average を作るため、海岸 cell は `0<f_land<1` になる。高周波地形の
smoothing は ADR 0008 の一種類の bounded filter だけを再利用し、pass count を preset に明記する。

Earth product の gate は次とする。

- 全値 finite、`0<=f_land<=1`、少なくとも一つずつ pure land、pure ocean、mixed coast cell がある
- `sum(A*f_land)` が派生 source product の全球 land area と登録許容差内で一致
- constant fraction と constant height を remap すると exact constant
- water-only target cell の surface height は exact zero
- seam/corner を含む回転テストで global mean と histogram が登録許容差内で不変
- Earth preset の surface pressure/layer thickness が地形込みで正

## 5. forcing を二経路に分離する

### 5.1 analytic planetary benchmark

`physics.kind=planetary_newtonian` を追加し、`forcing.geometry` を次の二値に限定する。

```text
forcing.geometry = axisymmetric | substellar
```

`axisymmetric` は Held--Suarez と同じ profile を任意の `Omega` で使う rotation-sensitivity、
`substellar` は Merlis--Schneider/Heng 型 tidally locked benchmark に使う。`chi` を
axisymmetric では `-sin(phi)^2`、substellar では符号を保った `dot(r_hat,s_hat)` として、共通式を

```text
T_eq = max(200 K,
           (315 K + 60 K*chi
                  - 10 K*log(p/p0)*cos(phi)^2) * (p/p0)^kappa)

k_T = 1/(40 day) + (1/(4 day)-1/(40 day))
                     * max(0,(sigma-0.7)/0.3) * cos(phi)^4
k_v = 1/(1 day) * max(0,(sigma-0.7)/0.3)

dT/dt = -k_T*(T-T_eq)
du/dt = -k_v*u
```

と固定する。`day=86400 s` であり rotation/orbit period に読み替えない。substellar の `chi` は night
side でもゼロ clip しない。これは published longitudinal dipole の負側を保持するためである。
`held_suarez` の意味は変更せず、axisymmetric kernel が同じ入力で既存 kernel と rounding 内で一致する
property test を置く。substellar direction は固定 longitude を config に直書きせず、同じ orbit kernel
から作る。

この経路は surface temperature を持たず、published analytic benchmark の比較専用である。
`surface.geography=uniform`、flat terrain だけを受理する。slow/rapid rotator は同じ axisymmetric forcing で
`Omega` だけを一要因変更し、tidally locked は `Omega=n` と substellar forcing を必須にする。

### 5.2 zero-depth surface energy balance

`physics.kind=surface_energy_balance` のときだけ `surface_temperature_k[C]` を予報する。cell `c` で

```text
Q_sw = (1-alpha) * S(t) * mu0
H = lambda_air * (T_s - T_air,bottom)

C_surface * dT_s/dt = Q_sw + Q_internal
                         - epsilon*sigma_SB*T_s^4 - H

d(M*theta)_bottom/dt += H / (cp * Pi_bottom)
```

とする。`H>0` は surface から大気への upward heat flux である。`alpha`、`epsilon`、`lambda_air`、
`Q_internal` は全球共通で、陸海差を付けない。Stefan--Boltzmann 定数はコード定数とし設定化しない。
大気 source は lowest model layer だけへ入り、mass、momentum、tracer、surface pressure を直接変えない。
surface mode の低層 Rayleigh drag は上記 `k_v` profile を再利用して全球一様に適用し、`f_land` へ
依存させない。drag work は surface reservoir へ戻さず、外部 energy sink として別に報告する。

これは短波の大気吸収と downwelling longwave を省いた idealized energy balance であり、現実の Earth
surface temperature や THAI climate を再現する放射 scheme とは呼ばない。Earth geography preset は
地形・海岸・日変化・熱慣性の統合試験であり、観測 climate の validation ではない。

surface RHS と大気 RHS は同じ stage state/time から SSP-RK3 の各 stage で再評価する。radiative cooling
が temperature positivity を破る場合は step を半減して retry し、temperature clip で energy error を
隠さない。surface-only scalar problem の time-step refinement で三次精度を確認する。

### 5.3 state、checkpoint、出力

既存 `DryHydrostaticState` の atmospheric arrays と flatten 順は変えない。surface mode だけが使う
optional `surface_temperature_k[C]` を末尾に持たせ、空 vector は従来 state とする。

- 従来 run: layout `dry_hydrostatic_cell_column_v1`、payload/FrameV2/canonical config を exact 維持
- surface run: 新 layout `dry_hydrostatic_surface_v1`、末尾に `T_s[C]`
- fixed `Phi_s`、`f_land`、`C_surface` は checkpoint に複製せず config/data fingerprint から再構築
- standalone output は `surface_state.csv` と `surface_diagnostics.csv`
- FrameV2 と gateway allowlist は変更せず、surface mode の machine-control run は Phase 8 では拒否

`surface_state.csv` は final cell field とし、`panel,i,j,longitude_deg,latitude_deg,height_m,land_fraction,`
`heat_capacity_j_m2_k,surface_temperature_k` を deterministic order で書く。

## 6. energy budget と無次元 metadata

surface reservoir energy anomaly は任意だが固定した reference `T_ref` に対して

```text
E_surface = sum_c A_c * C_surface[c] * (T_s[c] - T_ref)
```

とする。各 diagnostic interval で absorbed stellar power、internal heat、outgoing longwave、sensible
exchange、surface-energy change、大気への discrete thermal contribution、combined residual を分けて
出す。`H` の surface loss と大気 source は同じ stage-weighted flux から集計し、別式で再計算しない。
大気 total-energy への寄与は Phase 7 と同じ directional-derivative test で検証する。

dimensionless metadata は暗黙の代表 scale を使わない。Phase 8 preset は次を明示する。

```text
scales.characteristic_wind_m_s = U
scales.characteristic_length_m = L
scales.buoyancy_frequency_s_1 = N
scales.radiative_time_s = tau_rad
```

これと `T_ref` から pressure scale height `H_scale=Rd*T_ref/g`、reference latitude 45 deg の
`f_ref=sqrt(2)*abs(Omega)` を用いて次を記録する。

```text
Ro = U / (f_ref*L)
L_D = N*H_scale / f_ref
Bu = (L_D/L)^2
thermal_Ro = Rd*DeltaT / (Omega^2*R^2)
tau_rad_over_rotation = tau_rad / rotation_period
tau_surface(c) = C_surface[c] / (4*epsilon*sigma_SB*T_ref^3)
```

`Omega=0`、`N=0`、synchronous stellar day のような特異ケースは `inf`/NaN を数値 field に書かず、
値と `defined=false` を併記する。run 後には同じ定義のほか mass-weighted `U_rms` を使った diagnosed
Rossby number も別名で出し、input scale と観測値を混同しない。

## 7. preset と validation matrix

### 7.1 named preset

1. `phase8_earth_like.cfg`
   — Phase 7 Held--Suarez config の名前付き複製。surface/orbit key を追加せず、state と統計を exact regression。
2. `phase8_earth_geography.cfg`
   — Earth constants、Earth `Phi_s/f_land`、surface energy balance、季節なしの短い diurnal integration。
3. `phase8_slow_rotator.cfg`
   — flat uniform surface、axisymmetric analytic forcing、Earth-like case の `Omega/10`。
4. `phase8_rapid_rotator.cfg`
   — 同じ forcing で Earth-like case の `4*Omega`。
5. `phase8_tidally_locked.cfg`
   — flat aquaplanet (`f_land=0`)、`Omega=Omega_earth/365=n`、zero obliquity/eccentricity、
   substellar analytic forcing。

倍率や benchmark 定数は ADR 0010 で published comparison と pilot cost を見ずに事前登録する。
Earth geography と analytic Earth-like を同名にせず、現実地形を使うことと Earth-like benchmark である
ことを区別する。

### 7.2 unit/property/CI gate

- orbit/star geometry の解析値と global insolation integral
- parser の conditional required key、range、invalid pairing、old-config identity
- `f_land=0,1,0.25,0.5,0.75` の capacity interpolation と、閾値化されない mixed cell
- `C_ocean>C_land` のとき同じ step forcing に対する ocean の `|dT/dt|` が小さい
- no-forcing equilibrium、surface-only radiation、surface--air exchange の符号と time convergence
- sensible exchange だけなら combined surface+dynamics source budget が tolerance 内で閉じる
- old dry state/checkpoint/FrameV2 が byte-for-byte 不変、新 surface layout の restart が uninterrupted と一致
- Earth data の range、land area、height、seam/corner、fingerprint tests
- GCC/Clang、ASan/UBSan、format-check と全 Phase 0--7 regression

### 7.3 production gate

- Phase 7 Earth-like preset が同じ seed/resolution/time step で既存 baseline を変えない
- 無次元相似 case を無次元の長さ倍率 `lambda`、時間倍率 `tau`、圧力倍率 `P` で変換し、
  `R'=lambda R`、`Omega'=Omega/tau`、`g'=lambda g/tau^2`、velocity=`lambda/tau`、
  temperature-energy scale=`lambda^2/tau^2`、forcing time=`tau` とした normalized solution が一致
- slow/rapid rotator は一要因差として jet 数、緯度、Rossby/Burger 数、mass/energy budget を比較
- `Omega -> 0` case が Coriolis tendency zero を保ち、未定義 metadata を安全に出す
- tidally locked case は day/night temperature と lower/upper wind pattern を公開 benchmark envelope と比較
- Earth geography run は mixed coast、mountain pressure-gradient、diurnal land/ocean amplitude/phase、budget を報告
- CFL、surface/layer temperature、surface pressure、layer thickness、wind の preregistered safety envelope を
  外れた run は明示的に失敗し、clip や silent dt underflow を許さない

## 8. 変更境界

予定する新規責務は次とする。実装時の最終名は ADR 0010 で固定する。

```text
include/myplanetsim/core/orbit_parameters.hpp
src/core/orbit_parameters.cpp
include/myplanetsim/dynamics/surface_boundary.hpp
src/dynamics/surface_boundary.cpp
include/myplanetsim/physics/planetary_newtonian.hpp
src/physics/planetary_newtonian.cpp
include/myplanetsim/physics/surface_energy_balance.hpp
src/physics/surface_energy_balance.cpp
include/myplanetsim/diagnostics/planetary_scales.hpp
src/diagnostics/planetary_scales.cpp
data/earth/manifest.txt
data/earth/earth_surface_derived.csv
tools/build_earth_surface.py
tests/test_orbit_parameters.cpp
tests/test_surface_boundary.cpp
tests/test_surface_energy_balance.cpp
tests/test_planetary_newtonian.cpp
tests/test_planetary_scales.cpp
docs/adr/0010-planetary-forcing-and-fractional-surface.md
docs/validation/phase-8.md
```

既存変更は config parser/writer、surface orography adapter、dry state の optional tail、dry driver の complete
stage RHS、checkpoint layout selection、standalone CLI/CSV、CMake/test registration に留める。Earth generator は
offline tool であり、native executable に Python/GDAL dependency を加えない。

## 9. コミット列

1. `P8.01 docs: fix planet forcing and fractional-surface contracts`
   — ADR 0010、式、pairing、preset constants、reference envelope、Earth data provenance を固定。
2. `P8.02 core: add orbital geometry and planetary-scale metadata`
   — signed rotation、Kepler solver、star vector、stellar day、dimensionless diagnostics。
3. `P8.03 surface: add immutable land-fraction boundary`
   — `SurfaceBoundary`、uniform fraction、existing orography adapter、mixed-cell properties。
4. `P8.04 data: derive the reproducible Earth surface product`
   — ETOPO/GSHHG manifest、offline generator、Earth remap/range/area tests。
5. `P8.05 physics: add analytic planetary benchmark forcing`
   — axisymmetric/substellar kernels、conserved tendencies、slow/rapid/synchronous configs。
6. `P8.06 physics: add zero-depth surface energy balance`
   — `T_s`、capacity mixture、radiation/sensible RHS、source budgets、source-only convergence。
7. `P8.07 io: checkpoint and report the optional surface state`
   — new layout、restart equality、surface CSV、old layout/FrameV2 identity。
8. `P8.08 test: gate the coupled planetary core`
   — CI short runs、Earth geography health、similarity test、toolchains/sanitizers。
9. `P8.09 validation: run the bounded planetary matrix`
   — Earth regression、slow/rapid/tidally locked published comparison、known gaps/costs。

```text
entry -> P8.01 -> P8.02 -> P8.03 -> P8.04
                    |                  |
                    +-> P8.05          +-> P8.06 -> P8.07 -> P8.08 -> P8.09
```

## 10. 完了チェックリスト

### Compatibility and provenance

- [ ] Phase 5--7 production gap と ADR 0008 の扱いが validation report に明記されている。
- [ ] old config、Held--Suarez state/checkpoint/FrameV2、統計が exact regression する。
- [ ] Earth raw/derived data の出典、version、license、hash、生成 command が再現できる。

### Surface and forcing

- [ ] 全 cell の `0<=f_land<=1` と coast の `0<f_land<1` が保持され、二値 threshold がない。
- [ ] 海は depth/state/transport を持たず、陸海差は面積熱容量だけである。
- [ ] surface--atmosphere coupling が各 SSP-RK3 stage に一度だけ入り、budget が閉じる。
- [ ] Earth geography が地形と land fraction を同じ provenance から供給する。
- [ ] analytic benchmark と surface model の結果を混同していない。

### Planet validation

- [ ] Earth-like regression、similarity、`Omega -> 0` が登録済み gate を満たす。
- [ ] slow/rapid/tidally locked の設定と再現 command、平均/spread、reference comparison がある。
- [ ] Rossby/Burger/deformation radius/time-scale metadata の定義と undefined case が明示される。
- [ ] restart、CFL、positivity、GCC/Clang、ASan/UBSan、全 regression が通る。

## 11. 参照

- NOAA National Centers for Environmental Information (2022),
  [ETOPO 2022 Global Relief Model](https://www.ncei.noaa.gov/products/etopo-global-relief-model):
  Earth surface elevation の固定 raw source。
- NOAA NCEI,
  [Global Self-consistent, Hierarchical, High-resolution Geography Database](https://www.ngdc.noaa.gov/mgg/shorelines/shorelines.html):
  fractional land/water area を作る固定 shoreline polygon source。
- K. Heng, K. Menou, and P. J. Phillipps (2011),
  [Atmospheric circulation of tidally locked exoplanets: a suite of benchmark tests for dynamical solvers](https://doi.org/10.1111/j.1365-2966.2011.18315.x):
  dry tidally locked benchmark と inter-model comparison。
- T. M. Merlis and T. Schneider (2010),
  [Atmospheric Dynamics of Earth-Like Tidally Locked Aquaplanets](https://doi.org/10.3894/JAMES.2010.2.13):
  substellar longitudinal forcing と tidally locked circulation benchmark。
- D. E. Sergeev et al. (2023),
  [Simulations of idealised 3D atmospheric flows on terrestrial planets using LFRic-Atmosphere](https://doi.org/10.5194/gmd-16-5601-2023):
  Held--Suarez、Earth-like、tidally locked Earth の modern cubed-sphere comparison。

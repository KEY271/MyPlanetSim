# Phase 7 実装計画 — 理想化乾燥物理と気候統計

**状態:** 計画済み・未着手。Phase 7 のコード変更は、
[Phase 5 検証報告](validation/phase-5.md) と
[Phase 6 検証報告](validation/phase-6.md) の known gap を閉じ、
[ADR 0006](adr/0006-dry-hydrostatic-state-and-coupling.md) の時間積分契約と実装を一致させてから
開始する。結果は `docs/validation/phase-7.md` に記録する。

## 1. 到達点と開始条件

Phase 7 は、検証済みの乾燥・静水圧 3D core に Held--Suarez (1994) の Newtonian
temperature relaxation と低層 Rayleigh drag を一つだけ加え、1200 日積分のうち最初の 200 日を
spin-up として除外した気候統計を生成する。目的は汎用 physics package を作ることではなく、力学 core
が長時間の強制・散逸下でも既知の統計的循環を再現するかを測ることである。

```text
Phase 5/6 production validation + ADR 0006 implementation audit
                              |
                              v
Milestone A: one forcing contract
  Held--Suarez kernel -> conserved-variable tendencies -> source budgets
                              |
                              v
Milestone B: coupled short-run gate
  every SSP-RK3 stage -> invariants -> restart before statistics window
                              |
                              v
Milestone C: bounded climate report
  online zonal/time moments -> three seeds -> one-factor sensitivities
  -> Phase 7 validation report
```

開始前に次を満たす。

- Phase 5 の quantitative 3D convergence と UMJS14 reference envelope を完了する。
- Phase 6 の DCMIP 2-0-0、Williamson 5、linear mountain wave、JW06 の production gate を完了し、
  ADR 0008 を accepted にする。
- `DryHydrostaticDriver` が ADR 0006 どおり complete RHS を SSP-RK3 の各 stage で再評価することを
  test で固定する。現在の単段更新を Held--Suarez forcing で覆わない。
- unforced/flat run の mass、energy、AAM、restart の基準値と許容値を先に固定する。

これらは Phase 7 の成果に数えない。前段の不整合を relaxation や drag で安定化して完了扱いしない。

Phase 7 の成果物は次に限定する。

- `none` と `held_suarez` の二値だけを持つ idealized dry-physics selection
- 論文定数を固定した一つの Held--Suarez forcing kernel
- potential-temperature mass と horizontal momentum への source tendency
- thermal relaxation、Rayleigh work、model-energy residual の budget
- flat Earth-like Held--Suarez initializer と seed 付き微小 temperature perturbation
- cell-area/time-weighted の帯状平均、eddy moments、Hadley mass streamfunction
- 一つの budget 時系列 CSV、一つの最終 climate-statistics CSV
- C++ unit/property/short-run regression と production validation report

## 2. スコープ境界

### 2.1 含めるもの

- flat surface、perpetual equinox、乾燥大気の標準 Held--Suarez case
- full-level の `sigma = p_full/ps` による温度緩和率と低層 drag
- 既存の `M*theta`、`M*u` state へ直接変換した source
- dynamics と同じ SSP-RK3 stage で毎回評価する unsplit coupling
- `run.random_seed` だけで決まる、実装間で再現可能な微小 temperature perturbation
- cubed-sphere cell centre を equal-area latitude bin へ集約する online moments
- 1200 日 run、200 日 spin-up、1000 日統計窓
- baseline 3 seed と、時間刻み・解像度・水平拡散を一要因ずつ変える production comparison

### 2.2 含めないもの

- radiation solver、surface energy balance、soil/ocean、seasonal/diurnal cycle
- moist species、condensation、cloud、latent/sensible heat flux
- boundary-layer diffusion、vertical mixing、dry convective adjustment
- gravity-wave drag、top sponge、stochastic physics
- terrain と Held--Suarez forcing の組合せ
- configurable な equilibrium-temperature profile や timescale 群
- generic process registry、plugin ABI、operator-splitting framework
- physics coupling interval、subcycling、implicit/exact relaxation integrator
- arbitrary ensemble scheduler、parameter-sweep DSL、database、NetCDF
- full 4D history、spectral transform、EOF、formal uncertainty quantification
- Web UI、FrameV2、gateway、control protocol の変更
- Phase 8 の slow/rapid rotator、tidally locked、惑星別 forcing preset

Phase 7 では coupling を各 SSP-RK3 stage に固定する。変更可能な coupling interval を追加してから感度を
調べるのではなく、まず time-step halving で結合誤差を測る。別 interval が必要だという測定結果が出た場合に
だけ、後続 ADR で split/subcycle を比較する。

## 3. forcing contract

### 3.1 named physics selection

設定には一つの additive selection だけを加える。

```text
physics.kind = held_suarez
```

key がなければ `none` とし、既存 Phase 0--6 config の canonical text と fingerprint を変えない。
`held_suarez` は `experiment.kind=dry_hydrostatic`、`dry_hydrostatic.test_case=held_suarez`、flat
orography、登録済み Earth-like planet constants との組合せだけを受理する。論文定数を個別 config key に
展開しない。

### 3.2 equations and fixed constants

cell `c`、level `k` で `phi=asin(z_c)`、`sigma=p_full/ps`、`kappa=Rd/cp` とする。
標準 forcing は次の形を用いる。

```text
T_eq = max(T_strat,
           (T0 - Delta_T_y*sin(phi)^2
               - Delta_theta_z*log(p/p0)*cos(phi)^2) * (p/p0)^kappa)

k_T = k_a + (k_s-k_a)
            * max(0, (sigma-sigma_b)/(1-sigma_b)) * cos(phi)^4
k_v = k_f * max(0, (sigma-sigma_b)/(1-sigma_b))

dT/dt|HS = -k_T * (T-T_eq)
du/dt|HS = -k_v * u
```

定数は Held--Suarez の標準値 `T_strat=200 K`、`T0=315 K`、`Delta_T_y=60 K`、
`Delta_theta_z=10 K`、`tau_a=40 day`、`tau_s=4 day`、`tau_f=1 day`、
`sigma_b=0.7` に固定する。`day` は 86400 SI seconds であり、惑星の rotation period へ
読み替えない。

physics kernel は derived primitive と planet constants を読み、次だけを返す。

```text
HeldSuarezTendency
  potential_temperature_mass_k_kg_m2_s[C*K]
  horizontal_momentum_mass_kg_m_s2[C*K]
  diagnostic source rates
```

`M`、`ps`、tracer mass は変更しない。temperature tendency は同じ full-level Exner function で
`M*theta` tendency へ変換する。drag は Cartesian tangent velocity に scalar rate を掛け、source 後も
tangent であることを property test にする。kernel は state を直接 mutate しない。

### 3.3 initializer and seed

原論文は initial temperature を固定していないため、Heng et al. の再現計算に合わせて flat hydrostatic
column、zero wind、uniform 264 K を基準とし、その選択を ADR に明記する。baroclinic eddy を励起する
temperature perturbation は amplitude と生成法を ADR に固定し、自由 knob にしない。
`run.random_seed` と integer `(cell,k)` から明示した整数 hash で値を作り、標準ライブラリの
distribution 実装差に依存させない。perturbation は level ごとに area-weighted mean を除き、同じ seed
では GCC/Clang 間で同じ初期 field になることを test する。

## 4. physics--dynamics coupling

### 4.1 one unsplit path

各 SSP-RK3 stage で次を行う。

```text
stage state
  -> existing dynamics diagnostics and tendencies
  -> Held--Suarez tendency from the same stage state
  -> add M*theta and M*u sources once
  -> stage update, tangent projection, invariant validation
```

physics は mass continuity、vertical mass flux、pressure gradient、CFL estimate を変更しない。
`none` では source arrays を構築せず、RHS、stable dt、stage state、checkpoint、出力が Phase 6 と exact
に一致する。forcing の最短 timescale は 1 day であり、既存 dynamics time step より十分長い。Phase 7
では exact exponential update や split solver を追加しない。

### 4.2 source budgets

各 diagnostic sample で少なくとも次を出す。

- relaxation による `M*theta` source の global integral
- Rayleigh drag の zonal/meridional momentum tendency と kinetic-energy sink
- current discrete total-energy definition に対する physics source contribution
- measured total-energy change、physics contribution、既存 numerical/dynamics contribution の residual
- mass、tracer mass、minimum temperature、maximum wind、non-finite count

source contribution は実際に `M*theta` と `M*u` へ加えた discrete tendency から計算し、別の連続式を
「期待値」として二重管理しない。Newtonian relaxation は外部 energy source/sink、Rayleigh drag は外部
kinetic-energy sink なので total energy の保存は要求しない。要求するのは符号、単位、source attribution、
時間積分した budget の閉じ方が説明できることである。

`physics_diagnostics.csv` は `diagnostics.interval_steps` ごとに書き、瞬時 rate と累積 contribution を
分ける。production acceptance に使う累積 budget は統計と同じ day 200 でゼロから開始し、spin-up 中は
瞬時値だけを健全性確認に使う。無制限な per-step field output は作らない。

## 5. climate statistics

### 5.1 online accumulation

統計は full state history を保存せず、accepted step ごとに実時間幅 `dt` を重みにして online に更新する。
最初の 200 日は accumulation から除外し、day 200--1200 の 1000 日だけを使う。latitude bins は
`mu=sin(phi)` を等間隔に `2*N` 分割して面積をほぼ等しくし、cell は centre が属する bin に一度だけ
入れる。panel/seam/corner を別経路にしない。

各 accepted step で bin、native model level ごとの area-weighted zonal moment を先に作り、その値を
`dt` で時間積分する。次の十分統計量だけを蓄積する。

- mean pressure、temperature、zonal wind、meridional wind
- `u^2`、`v^2`、`T^2`、`u*v`、`v*T`
- layer mass と meridional mass transport

prime は各時刻の zonal mean からの偏差と定義し、その covariance/variance を time-average する。ここから
time/zonal mean、`u'v'`、`v'T'`、eddy kinetic energy、eddy temperature variance、Hadley mass
streamfunction を導出する。風は planet rotation axis に対する east/north basis へ一度だけ射影する。
Phase 7 preset の native hybrid levels をそのまま使い、任意 pressure grid への補間器を作らない。

### 5.2 output contract

最終 `climate_statistics.csv` は一行を `(latitude_bin, level)` とし、少なくとも次を含む。

```text
latitude_lower_deg,latitude_upper_deg,level,mean_pressure_pa,
mean_temperature_k,mean_zonal_wind_m_s,mean_meridional_wind_m_s,
eddy_momentum_flux_m2_s2,eddy_heat_flux_k_m_s,
eddy_kinetic_energy_m2_s2,temperature_variance_k2,
mass_streamfunction_kg_s,accumulated_time_s
```

順序、units、normalization、prime の定義を ADR と列名に固定する。CSV と既存 run
metadata/config fingerprint だけで結果を再現できるようにし、画像を acceptance artifact にしない。

### 5.3 restart boundary

atmospheric checkpoint layout は変更しない。restart equality は spin-up 中の登録済み checkpoint から
day 1200 まで再開し、uninterrupted run と最終 state、budget、day 200--1200 statistics が一致することで
検証する。統計窓の途中から accumulator を再開するための別 checkpoint schema は Phase 7 に含めない。
その必要性が実運用で確認された場合は、state checkpoint と analysis state を別 contract として設計する。

## 6. benchmark と acceptance

### 6.1 unit/property gate

- `T_eq`、`k_T`、`k_v` が equator、pole、`sigma_b` の上下、stratospheric floor で論文式に一致する。
- `T=T_eq` かつ `u=0` で source が exact zero になる。
- drag work は非正で、upper atmosphere では exact zero になる。
- physics は `ps`、air mass、tracer mass を変えず、momentum tangent property を保つ。
- `physics.kind=none` は既存 dry RHS/stage/checkpoint/canonical config と exact identity になる。
- 同じ seed の initializer と statistics reduction が toolchain 間で再現する。
- constant/zonal analytic fields の bin mean、covariance、streamfunction が手計算値に一致する。
- online accumulation を複数 chunk に分けても一括処理と許容誤差内で一致する。

### 6.2 CI short-run gate

CI は climate conformance を主張せず、小さい `N/K` と短期間で次だけを gate する。

- forcing ありで temperature が equilibrium profile 側へ進み、低層 wind が減衰する。
- SSP-RK3 time-step refinement で source-only manufactured solution が設計次数を示す。
- mass/tracer invariants、temperature floor、surface-pressure bounds、NaN/Inf absence を維持する。
- budget CSV と climate CSV の row count、units、finite、deterministic ordering を検査する。
- spin-up 前 checkpoint からの restart が uninterrupted 結果と一致する。
- all Phase 0--6 tests、GCC/Clang、ASan/UBSan、format-check が通る。

### 6.3 production climate gate

標準 run は 1200 日、最初の 200 日を破棄し、残り 1000 日を平均する。実装前の ADR で
Held--Suarez 原論文と複数の公開実装から比較量、pressure/latitude location、reference envelope を
引用付きで登録し、結果を見た後に threshold を緩和しない。少なくとも次を比較する。

- time/zonal mean temperature と zonal wind の構造、extrema、位置
- Hadley cell の符号、強度、緯度・鉛直範囲
- eddy momentum flux、eddy heat flux、eddy kinetic energy、temperature variance の peak と位置
- hemispheric symmetry と、統計窓前半/後半の stationarity
- forcing/drag の累積 contribution と model-energy residual

production matrix は組合せ sweep にせず、次の 6 run に限定する。

1. 登録済み baseline resolution/time step/diffusion の seed 0。
2. 同じ baseline の seed 1。
3. 同じ baseline の seed 2。
4. seed 0 で time step を 1/2。
5. seed 0 で水平・鉛直を一段だけ高解像度化。
6. seed 0 で水平 diffusion を一段だけ変更。

P7.01 では baseline 候補と採否基準を登録する。P7.06 の short/pilot run で wall time、安定性、最低限の
eddy activity を確認し、1200 日の結果を見る前に baseline resolution を一つに固定する。追加の factorial
sweep は行わない。3 seed の spread は内部変動として report し、
感度差を単一 snapshot の差と混同しない。

## 7. 変更境界

予定する新規責務は次に限定する。実装時に実際の命名を ADR で確定する。

```text
include/myplanetsim/physics/held_suarez.hpp
src/physics/held_suarez.cpp
include/myplanetsim/diagnostics/climate_statistics.hpp
src/diagnostics/climate_statistics.cpp
tests/test_held_suarez.cpp
tests/test_climate_statistics.cpp
configs/phase7_held_suarez.cfg
docs/adr/0009-held-suarez-coupling-and-statistics.md
docs/validation/phase-7.md
```

既存変更は config parser/writer、dry benchmark initializer、dry driver の stage RHS、CLI の二つの CSV
writer、CMake/test registration に留める。`DryHydrostaticState`、checkpoint payload、FrameV2、`web/` は
変更しない。将来候補の空 interface や directory は作らない。

## 8. コミット列

一つの commit は一つの責務と test を含む。entry gate が閉じる前に P7.02 以降へ進まない。

1. `P7.01 docs: fix the Held--Suarez contract and reference envelopes`
   — ADR 0009 で式、定数、stage coupling、統計定義、production cost と閾値を固定。
2. `P7.02 config: add one idealized dry-physics selection`
   — additive `none|held_suarez`、strict pairing、initializer、seed reproducibility。
3. `P7.03 physics: implement Held--Suarez tendencies`
   — temperature relaxation、Rayleigh drag、conserved-variable conversion、unit/property tests。
4. `P7.04 dynamics: couple forcing at every SSP-RK3 stage`
   — complete RHS integration、none identity、source budgets、source-only time convergence。
5. `P7.05 diagnostics: accumulate bounded climate statistics`
   — equal-area latitude bins、online moments、eddy quantities、streamfunction、CSV contract。
6. `P7.06 test: gate the forced dry core`
   — short coupled run、invariants、budget attribution、pre-spin-up restart、toolchains。
7. `P7.07 validation: run the bounded climate matrix`
   — 3 seed + 3 one-factor comparisons、published envelopes、cost/commands/known gaps。

```text
entry gate -> P7.01 -> P7.02 -> P7.03 -> P7.04 -> P7.05 -> P7.06 -> P7.07
```

## 9. 完了チェックリスト

### Entry

- [ ] Phase 5/6 の production validation gap が閉じている。
- [ ] ADR 0006 の SSP-RK3/stage/CFL 契約と dry driver 実装が一致している。
- [ ] ADR 0008 が accepted で、flat dry baseline が再生成できる。

### Physics and coupling

- [ ] standard Held--Suarez formula/constants が unit test と ADR に固定されている。
- [ ] source が `M*theta` と `M*u` に一度だけ入り、mass/tracer/terrain を変更しない。
- [ ] `none` が Phase 6 の exact identity である。
- [ ] thermal/drag/model-energy budget の符号、単位、residual を説明できる。
- [ ] source-only と coupled run の time-step refinement が登録済み gate を満たす。

### Statistics and validation

- [ ] day 200--1200 の面積・時間重み付き統計が一つの CSV から再生成できる。
- [ ] mean `T/u`、Hadley cell、四つの eddy statistics が公開 envelope と比較されている。
- [ ] 3 seed の平均と spread、3 one-factor sensitivity が報告されている。
- [ ] spin-up 前 restart が最終 state、budget、統計を変えない。
- [ ] CI gate と production gate を区別し、command、wall time、threshold、known gap を記録している。
- [ ] `DryHydrostaticState`、checkpoint payload、FrameV2、`web/` が変更されていない。

## 10. 参照

- I. M. Held and M. J. Suarez (1994),
  [A Proposal for the Intercomparison of the Dynamical Cores of Atmospheric General Circulation Models](https://www.gfdl.noaa.gov/bibliography/related_files/ih9401.pdf):
  standard Newtonian relaxation、Rayleigh friction、1200-day/200-day spin-up benchmark。
- H. Wan et al. (2008),
  [Ensemble Held--Suarez Test with a Spectral Transform Model: Variability, Sensitivity, and Convergence](https://doi.org/10.1175/2007MWR2044.1):
  1000-day mean の内部変動、ensemble、解像度感度を acceptance の解釈に使う。
- K. Heng, K. Menou, and P. J. Phillipps (2011),
  [Atmospheric circulation of tidally locked exoplanets: a suite of benchmark tests for dynamical solvers](https://doi.org/10.1111/j.1365-2966.2011.18315.x):
  Earth Held--Suarez の独立実装比較。tidally locked cases 自体は Phase 8 へ残す。
- E. Sergeev et al. (2023),
  [Simulations of idealised 3D atmospheric flows on terrestrial planets using LFRic-Atmosphere](https://doi.org/10.5194/gmd-16-5601-2023):
  cubed-sphere 系での standard formula、長時間平均、現代的な reference comparison。

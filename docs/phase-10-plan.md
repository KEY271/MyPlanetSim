# Phase 10 実装計画 — 重力波 semi-implicit 時間積分

**状態:** P10.15 進行中（2026-09-07）。P10.01--P10.14 は完了し、残る完了条件は 1200 日 run のみ。
水平 Lamb 波と高速な内部重力波だけを semi-implicit に扱う。`N=12`, `K=20` の標準的な乾燥大気で、
現在の 200 s から通常 1200--1800 s、目標 1800 s の大刻みへ移ることを狙う。

## 1. 到達点と設計原則

Phase 9 後の乾燥静水圧コアは、移流 Rusanov 散逸と Lamb CFL 速度を分離しているが、時間積分は
全項 SSP-RK3 のままである。`phase7_held_suarez` と同じ `N=12`, `K=20` の初期状態では、
per-cell Lamb CFL 上限は約 240.6 s であり、preset は 200 s を使う。1200 日積分は 518,400 step を
必要とする。1800 s を常用できれば 57,600 step となり、step 数は 1/9 になる。

Phase 10 の到達点を次に固定する。

1. explicit SSP-RK3 を回帰経路として維持したまま、乾燥静水圧コアに一段・二時刻レベルの
   iterative semi-implicit 経路を追加する。
2. `surface pressure -- horizontal divergence -- temperature -- hydrostatic pressure force` が作る
   高速波だけを reference-linear Jacobian に含め、現在の完全な非線形 RHS は変えない。
3. external/Lamb mode と、要求 step で陽的 wave CFL を超える最初の内部重力波 mode を陰的にする。
4. `run.time_step_s = 1800` を要求できるようにし、陽に残る移流、鉛直輸送、拡散、surface reservoir
   の制約に応じて安全に短縮する。
5. flat、terrain-following、baroclinic、forced climate の順に検証し、30 分を単なる「落ちない値」では
   なく、位相誤差・保存量・solver cost を伴う登録済み能力にする。

30 分は必達の運用目標であって hard ceiling ではない。1800/900/450 s の精度 gate を通した後、
2400、3600、5400、7200 s を順に探索し、安定性、位相・振幅、balance、climate envelope、retry、cost の
いずれかが最初に外れる点を case ごとの実用上限として記録する。単純な乾燥力学が30分超を許しても、
physics を結合した preset は独立により短い上限を持ち得る。

```text
Phase 9 explicit baseline
        |
        v
Milestone 0: 契約と基準測定
  ADR 0016 / explicit cost・CFL・wave spectrum
        |
        v
Milestone 1: 挙動を変えない分離
  RHS components / fast reference operator / CFL diagnostics
        |
        v
Milestone 2: solver を独立に検証
  matrix-free GMRES / preconditioner / shallow-water wave test bed
        |
        v
Milestone 3: dry external mode
  iterative Crank--Nicolson / 2-D Helmholtz / conservative back substitution
        |
        v
Milestone 4: multi-level gravity modes
  vertical structure matrix / modal transform / automatic mode selection
        |
        v
Milestone 5: nonlinear・terrain・physics integration
  retry / positivity / DCMIP / mountain wave / Held--Suarez
        |
        v
Milestone 6: 30-minute and performance gates
  1800/900/450 s convergence / 30-day comparison / 1200-day pilot
```

設計原則を六つ置く。

1. **半離散方程式を変えない。** `F(U)=L_ref(U-U_ref)+N(U)` とし、
   `N(U)=F(U)-L_ref(U-U_ref)` を同じ演算子で作る。reference 項を追加したまま nonlinear pressure
   force を二重に数えない。
2. **保存は後処理で直さない。** `ps` と `M*theta` の implicit correction は共有 edge flux の scatter
   から構成し、global mean を最後に合わせる mass fixer や clipping を使わない。
3. **explicit 経路を oracle として残す。** 既存 config に新 key を要求せず、既存 SSP-RK3 preset、
   checkpoint、Frame、診断を byte-exact に保つ。
4. **solver failure は accepted state にしない。** 線形残差未達、非線形反復未収束、positivity 違反、
   stage CFL 超過は step を半減して最初から再試行し、最小 step 未満なら診断付きで停止する。
5. **一段法を使う。** AB2/CN のような過去 tendency を必要とする方式は checkpoint layout を変えるため
   採用しない。restart は現在の瞬時 state だけから同じ次 step を再現する。
6. **安定性と速さを別々に gate する。** 1800 s で完走しても、一モデル日あたりの wall time が
   explicit より改善しなければ Phase 10 の性能完了とはしない。

## 2. 現状と解くべき stiffness

### 2.1 現在の時間積分

各 SSP-RK3 stage は次をすべて再評価する。

```text
state
  -> hybrid pressure / thermodynamics / hydrostatic diagnosis
  -> horizontal reconstructed fluxes
  -> column ps_dot and coordinate-relative vertical mass flux
  -> pressure-gradient and Coriolis sources
  -> diffusion and physics
  -> explicit stage update
```

水平 CFL は

```text
dt * sum_faces((abs(u_n) + sqrt(gamma*Rd*T)) * L_f) / A_cell <= cfl
```

であり、Earth-like 初期温度で fast speed は約 326 m/s である。1800 s は現在の初期上限に対して
wave Courant 約 3.4 に相当するため、SSP-RK3 のまま設定値だけを増やすことはできない。

一方 `dry_hydrostatic_flux` の jump dissipation は既に

```text
a_adv  = max(abs(u_n,L), abs(u_n,R))
a_fast = max(abs(u_n,L) + sqrt(gamma*Rd*T_L),
             abs(u_n,R) + sqrt(gamma*Rd*T_R))
```

へ分かれている。semi-implicit 経路では `a_adv` を explicit horizontal CFL に使い、`a_fast` は
implicit wave Courant の診断にだけ使う。explicit 経路は ADR 0011 の `a_fast` CFL を維持する。

### 2.2 30 分化後も残る制約

Phase 10 は次を陰的にしない。

| 制約 | Phase 10 の扱い |
|---|---|
| horizontal material advection | explicit。per-cell `a_adv` CFL を満たす |
| coordinate-relative vertical transport | explicit。既存 `vertical.cfl` を満たす |
| Laplacian/biharmonic diffusion | explicit。既存 diffusion limit を満たす |
| surface energy reservoir | explicit。既存 `surface.cfl` を満たす |
| Held--Suarez / planetary forcing | nonlinear iterationで再評価するが fast Jacobianには含めない |
| pressure and temperature bounds | accepted step 前に現在と同じ invariantを検査 |

したがって「1800 s」は無条件の固定 step ではなく最大要求値である。完了ゲートでは flat linear wave が
1800 s をそのまま受理することに加え、Held--Suarez で median accepted step が 1200 s 以上であることを
要求する。強い jet や極端な惑星設定で step が短縮されることは正しい挙動とする。

### 2.3 Milestone 0 で固定する基準測定

実装前に同一 compiler/build で次を保存する。

- `phase7_held_suarez` の day 0、1、10、30 における fast-wave、advective、vertical、diffusion、surface
  の各候補 step
- 200 s explicit run の seconds/model-day、cell-level-updates/s、peak RSS、RHS calls/step
- `phase5_linear_wave` の surface-pressure amplitude、phase、dry mass、energy drift
- small `N=2,4`, `K=2,4,8` reference column の線形化 Jacobian spectrum
- external mode と先頭 internal modes の phase speed、および 1800 s での mode 別 Courant number

この測定を `docs/validation/phase-10.md` の baseline 節へ先に記録し、後から target に合わせて基準を
動かさない。

## 3. 採用する semi-implicit 法

### 3.1 予報変数と reference state

checkpoint と prognostic layout は変えず、

```text
U = (ps, M*u, M*theta, M*q, optional Ts)
```

を用いる。`L_ref` が作用するのは `ps`, `M*u`, `M*theta` だけであり、passive tracer と surface
temperature に対する成分は 0 とする。

reference state は time-independent、水平一様な resting hydrostatic column とする。

- `ps_ref` は config に固定された正の値を使う
- `T_ref` は正の scalar とし、全層の reference temperature を定める
- `M_ref`, `theta_ref`, Exner、hydrostatic coefficients は `A/B` と惑星定数から driver construction 時に作る
- terrain は `L_ref` に含めず、full nonlinear residual と既存 pressure reference で扱う
- `L_ref(U_ref)=0` を unit test で要求する

架空惑星では Earth 固定 300 K を黙って使わない。reference temperature と pressure は canonical config と
fingerprint に入れる。

### 3.2 reference-linear fast operator

`L_ref` は現在の discrete path の tangent-linear 版として次を結合する。

1. `delta(M*u)` から centered reference mass flux を作る
2. unique edge scatter で level ごとの horizontal mass convergenceを得る
3. 既存 `B` recurrence の線形形で `delta ps_dot` と interface mass fluxを得る
4. reference `theta` の horizontal/vertical transportから `delta(M*theta)_dot` を得る
5. `delta ps`, `delta(M*theta)` から `delta T`, `delta Phi`, `delta alpha` を診断する
6. `delta grad(Phi) + alpha_ref*delta grad(p)` から pressure accelerationを得る

fast operator 自身は limiter、Rusanov jump、Coriolis、diffusion、physics、tracer を含めない。これらは
`N(U)` に残る。linear mass と thermal flux は同じ edge orientationを使い、global `ps` mass と
`M*theta` integral の変化が roundoff で 0 になることを property test にする。

### 3.3 iterative Crank--Nicolson

implicit weight `alpha` を使い、`alpha=0.5` を centered Crank--Nicolson とする。時刻 `n+1` の反復値
`U^(l)` に対し full residual を

```text
R(U^(l)) = U^(l) - U^n
           - dt * ((1-alpha)*F(U^n) + alpha*F(U^(l)))
```

とし、reference-linear quasi-Newton step

```text
[I - alpha*dt*L_ref] deltaU = -R(U^(l))
U^(l+1) = U^(l) + deltaU
```

を行う。初期値は `U^(0)=U^n` とする。

> **P10.14 で改訂。** 当初この節は「既定 2 回、残差未達なら上限 4 回、scaled nonlinear residual を
> 判定する」と定めていた。実装後の測定でこの残差判定が 1800 s を一度も受理せず、explicit より
> 2.4 倍遅い原因であることが判明した。Bénard (2003) の ICI 分類、Thuburn et al. (2014)、ENDGame の
> いずれも固定反復数を使い残差を判定しない。改訂後の契約は
> `semi_implicit.nonlinear_iterations` による固定反復で、残差は診断と発散検出にのみ使う。反復数は
> 反復誤差が時間打ち切り誤差を下回る最小値として実測で 3 に決めた。詳細は ADR 0016 と
> [Phase 10 検証報告](validation/phase-10.md) の P10.14 節を参照。

`alpha=0.5` で temporal convergence gate を通す。高周波 gravity-wave damping が必要な production
preset では `0.5 < alpha <= 0.55` を許すが、off-centering が時間精度と波の振幅を変えることを metadata
と validation に記録する。既定値の最終決定は P10.12 の flat/terrain 比較後に一度だけ行う。

### 3.4 vertical mode decomposition

reference-linear momentum、temperature、surface-pressure equations を消去して、level divergence の
3-D Helmholtz 系を作る。水平一様 reference により鉛直構造行列は全 cell で共通なので、一度だけ固有値
分解し、各 mode `m` の 2-D 問題へ変換する。

```text
(I - alpha^2 * dt^2 * c_m^2 * D*G) x_m = b_m
```

`c_m` は mode phase speed、`D*G` は現在の cell-centred flux divergence と pressure-gradient path に
対応する水平演算子である。最初は external mode だけを通して end-to-end を検証し、その後 internal
modes を追加する。

automatic selection は

```text
dt * c_m * sum_faces(L_f) / A_cell > semi_implicit.wave_cfl_threshold
```

となる mode を陰的対象とする。上限 mode 数を config で固定し、上限を超える必要が生じた run は
黙って残りを陽にせず拒否する。選択 mode 数、各 `c_m`、最大 implicit wave Courant を metadata に残す。

### 3.5 Helmholtz solver

現在の least-squares gradient と finite-volume divergence は厳密な随伴対ではないため、正定値対称性を
仮定する CG を最初から使用しない。

- outer solver: restarted matrix-free GMRES
- first preconditioner: area-scaled diagonal/Jacobi
- production preconditioner: symmetric finite-volume Helmholtz の multilevelまたは粗格子 correction
- workspace: driver-owned、warm-up 後 allocation 0
- stopping: scaled relative residualと absolute residualの両方
- failure: iteration上限到達時は step半減。未収束解を accept しない

`N=12` では 1 level が 864 cells であり、数個の modal 2-D solve は小さい。Phase 10 では MPI
multigrid を先回りせず、正しい operator、iteration count、wall timeを先に登録する。

### 3.6 back substitution と保存

modal solutionから divergence、temperature、`ln ps`、horizontal momentum correctionを戻し、最後に
prognostic `ps`, `M*u`, `M*theta` へ変換する。

- `ps` increment は解いた edge mass-flux correction の divergenceとして適用する
- `M*theta` は同じ mass flux と reference theta transportから作る
- momentum は cell中心 tangent planeへ射影する
- tracerは nonlinear finite-volume fluxだけで更新する
- solver residualを隠す global mean correction、temperature clipping、negative tracer clippingを行わない

線形 solve toleranceが保存 residualへ現れる場合は、toleranceを厳しくするか conservative flux formを
修正する。後処理で帳尻を合わせない。

## 4. 設定、互換性、診断

### 4.1 conditional config

最終名は ADR 0016 で固定するが、想定する設定は次である。

```ini
dry_hydrostatic.time_integrator = semi_implicit
dry_hydrostatic.advective_cfl = 0.45

semi_implicit.reference_surface_pressure_pa = 100000
semi_implicit.reference_temperature_k = 300
semi_implicit.reference_update = fixed
semi_implicit.implicit_weight = 0.5
semi_implicit.wave_cfl_threshold = 0.45
semi_implicit.maximum_implicit_modes = 5
semi_implicit.nonlinear_iterations = 3
semi_implicit.linear_relative_tolerance = 1e-4
semi_implicit.linear_maximum_iterations = 40
semi_implicit.minimum_time_step_s = 60

run.time_step_s = 1800
```

互換性規則。

- key が無い既存 config は SSP-RK3 のまま
- explicit mode では `dry_hydrostatic.cfl` が ADR 0011 の fast-wave CFL を意味し続ける
- semi-implicit modeだけが `dry_hydrostatic.advective_cfl` と `semi_implicit.*` を要求する
- writer は explicit legacy config に新 key を挿入せず、canonical text と fingerprintを変えない
- semi-implicit key は transport、shallow-water、vertical-column experimentでは拒否する
- tolerance、reference state、implicit weight、mode cap は結果へ影響するため fingerprint対象とする

### 4.2 step diagnostics

accepted stepごとに内部集計し、configured diagnostic intervalで次を出す。

```text
requested_dt_s
accepted_dt_s
advective_cfl
implicit_wave_cfl
vertical_cfl
selected_implicit_modes
linear_iterations_total / maximum
linear_relative_residual_max
nonlinear_iterations
nonlinear_relative_residual
retry_count
wall_seconds_rhs / linear_solve / total
```

`maximum_cfl` だけを出して implicit wave Courant と explicit stability Courant を混同しない。

### 4.3 preset

既存 Phase 5--8 preset は historical explicit baseline として変更しない。次を追加する。

```text
configs/phase10_linear_wave_semi_implicit.cfg
configs/phase10_dcmip_rest_semi_implicit.cfg
configs/phase10_linear_mountain_wave_semi_implicit.cfg
configs/phase10_held_suarez_semi_implicit.cfg
```

最後の preset は `N=12`, `K=20`, requested `dt=1800 s` とし、比較対象は
`configs/phase7_held_suarez.cfg` の explicit 200 s run とする。

## 5. 検証契約

### 5.1 algebra・solver unit gate

- zero RHS、identity、既知固有値を持つ小行列で GMRES の解と残差を検査
- iteration上限、NaN、zero norm、shape mismatchを明示的に拒否
- warm-up 後の solve が heap allocation 0
- vertical eigenvectorsの往復誤差と mode ordering
- `D*G` の numerical spectrumを small gridで測定し、solver選択と整合すること
- `L_ref(U_ref)=0`
- `F(U) = L_ref(U-U_ref) + N(U)` が componentごとに roundoffで一致
- fast operatorの global mass と `M*theta` tendencyが roundoffで 0

### 5.2 shallow-water test bed

production shallow-water driverの integratorは変更せず、同じ cubed-sphere geometry上の線形 shallow-water
系を solver test fixtureとして使う。

- gravity-wave Courant 0.5, 1, 2, 4, 8 で no-growth
- `alpha=0.5`、dt 半減列で2次時間収束
- low-wavenumber wave の振幅・位相誤差を explicit small-dt referenceと比較
- grid-scale divergent modeが増幅しない
- constant depth/rest と global massが roundoffで不変

この test bed が通る前に3-D dry driverへ solverを接続しない。

### 5.3 dry flat-wave gate

`phase5_linear_wave` を延長した flat isothermal caseで次を登録する。

- explicit 50 s referenceに対する 1800/900/450 s の surface-pressure amplitude・phase・L2 error
- `alpha=0.5` の observed temporal order
- 1800 s で最大 implicit wave Courantが3以上
- dry mass、`M*theta`、tracer mass driftが既存 explicit gateと同じ order
- pressure、temperature、layer thickness、tangencyが全 accepted stepで有効
- external-only と multi-mode の差を記録し、必要 internal mode数を事前基準から選ぶ
- 上記 gate 通過後に2400/3600/5400/7200 sを順次試し、最初の不合格理由と最後の合格 stepを記録する

### 5.4 balance、terrain、baroclinic gate

flat caseだけで accepted にしない。

1. isothermal rest: RHSとstateが roundoffで不変
2. DCMIP 2-0-0 rest: 既存 reference-balanced pressure forceを壊さず、implicit correctionが0
3. linear mountain wave: 30/15/7.5分列で波の振幅・位相・鉛直構造が収束
4. UMJS14/JW06 steady: 初期 balance residualを増幅しない
5. UMJS14/JW06 baroclinic: explicit small-step referenceに対する day 1--10 の norm/extrema envelope

terrain caseでだけ nonlinear/linear iterationが増える場合は、reference operatorへ terrainを安易に
混ぜず、pressure-gradient residualとpreconditionerを分けて診断する。realistic orographyで不安定なら
Phase 10 は flat-only success として完了扱いにしない。

### 5.5 forced climate と長時間 gate

- `phase10_held_suarez_semi_implicit` の30日 runを explicit 200 sと比較
- 最初の30日で median accepted step >= 1200 s
- 1800 sが受理された stepを明示し、短縮理由別の回数を報告
- zonal temperature/wind、Hadley circulation、mass/energy budgetを比較
- 同じ stateからの uninterrupted/restart が同じ step列とstateを再現
- 30日 gate後に1本の1200日 runを完走し、solver iteration、retry、climate statisticsの時間変化を記録

1200日結果の科学的 conformance判定は Phase 7 の責務を置き換えない。Phase 10 が要求するのは、
semi-implicit経路が長時間安定し、explicit referenceから説明不能に分岐せず、計算コストを実際に下げる
ことである。

### 5.6 performance gate

同じ `N=12`, `K=20`、同じ compiler、同じ診断間隔で model-dayあたりを比較する。

| 指標 | 登録目標 |
|---|---:|
| 30日 wall time speedup vs explicit 200 s | 3x 以上 |
| GMRES linear iteration | 最大40未満、p95 20以下 |
| nonlinear iteration | 最大4、通常2 |
| accepted median dt | 1200 s以上 |
| peak RSS | explicit の2倍以下 |
| warm-up後 solver allocation | 0 |

目標を外した場合は multigrid、mode cap、solver toleranceのどれが支配的かを profileし、数値許容値を
黙って緩めない。

## 6. スコープ境界

### 6.1 含めるもの

- dry-hydrostatic RHS component分離と reference-linear gravity-wave operator
- iterative Crank--Nicolson、vertical mode decomposition、2-D Helmholtz solve、back substitution
- matrix-free GMRESと少なくとも1つの実用 preconditioner
- explicit/implicit/vertical/surface CFLの分離診断
- conditional config、metadata、retry、restart、allocation/performance gate
- shallow-water linear test fixture
- flat、terrain、baroclinic、Held--Suarezの段階的 validation
- Phase 10専用 presetと `docs/validation/phase-10.md`

### 6.2 含めないもの

- semi-Lagrangian advection、departure-point探索、mass fixer
- production shallow-water solverの時間積分変更
- 非静水圧方程式、vertical acoustic solve、split-explicit subcycling
- vertical transport、diffusion、surface reservoir、physicsの全面的な implicit 化
- moist physics、灰色放射、境界層、乾燥対流調節
- MPI分散 multigrid、GPU、accelerator、並列 I/O
- 既存 pressure-gradient離散の全面交換、C-grid prognostic stateへの移行
- Phase 7の複数 seed climate conformance matrixの代行

これらは Phase 11 以後とする。Phase 10 中に advective CFL が 30 分化の主制約と実測された場合、
semi-Lagrangian化を同じ Phaseへ追加せず、測定を後続 ADR の入力として残す。

## 7. 変更境界

新規に想定するファイル。最終名は ADR 0016 で固定する。

```text
docs/adr/0016-semi-implicit-gravity-wave-integration.md
docs/validation/phase-10.md
include/myplanetsim/numerics/gmres.hpp
src/numerics/gmres.cpp
include/myplanetsim/numerics/helmholtz_solver.hpp
src/numerics/helmholtz_solver.cpp
include/myplanetsim/dynamics/dry_hydrostatic_fast_modes.hpp
src/dynamics/dry_hydrostatic_fast_modes.cpp
include/myplanetsim/dynamics/dry_hydrostatic_semi_implicit.hpp
src/dynamics/dry_hydrostatic_semi_implicit.cpp
tools/analyze_dry_wave_modes.cpp
tools/benchmark_semi_implicit.cpp
tests/test_gmres.cpp
tests/test_semi_implicit_shallow_water.cpp
tests/test_dry_hydrostatic_fast_modes.cpp
tests/test_dry_hydrostatic_semi_implicit.cpp
tests/test_phase10_balance.cpp
tests/test_phase10_long_step.cpp
configs/phase10_linear_wave_semi_implicit.cfg
configs/phase10_dcmip_rest_semi_implicit.cfg
configs/phase10_linear_mountain_wave_semi_implicit.cfg
configs/phase10_held_suarez_semi_implicit.cfg
```

既存変更は `experiment_config`、`dry_hydrostatic_driver`、`dry_hydrostatic_workspace`、
`dry_hydrostatic_coupling`、`dry_hydrostatic_sources`、run metadata/CSV、CMake/test登録、README、
`docs/plan.md` に留める。state/checkpoint/Frame layoutは変更しない。

## 8. コミット列

```text
1.  P10.01 docs: fix the semi-implicit gravity-wave contract
      ADR 0016、baseline測定、config、保存、retry、精度・性能閾値を固定

2.  P10.02 dynamics: separate dry RHS components and Courant diagnostics
      total RHSを変えずに pressure/continuity/thermal fast blockを識別、explicit byte-exact
3.  P10.03 numerics: add allocation-free restarted GMRES
      小行列、失敗経路、scaled norm、workspace test
4.  P10.04 numerics: validate Helmholtz solves on linear shallow water
      cubed-sphere fixture、large wave Courant、時間収束、mass/no-growth gate

5.  P10.05 config: add the conditional semi-implicit configuration
      legacy canonical text/fingerprint不変、invalid pairing、metadata
6.  P10.06 dynamics: build the hydrostatic reference column and vertical modes
      structure matrix、eigenpairs、mode ordering、automatic selection
7.  P10.07 dynamics: implement the conservative reference-linear fast operator
      ps/Mu/Mtheta tangent operator、add-subtract identity、global conservation

8.  P10.08 dynamics: integrate the external mode with iterative Crank--Nicolson
      2-D Helmholtz、back substitution、flat linear-wave 1800 s
9.  P10.09 dynamics: extend the solve to selected internal gravity modes
      modal transform、mode cap、multi-level phase/amplitude gate
10. P10.10 driver: add convergence, invariant, CFL, and retry control
      no unconverged accept、dt半減、minimum dt、restart determinism
11. P10.11 diagnostics: report solver work and separated Courant numbers
      CSV/metadata、allocation-free workspace、benchmark harness

12. P10.12 test: gate flat, balanced, terrain, and baroclinic dynamics
      DCMIP rest、mountain wave、UMJS14/JW06、GCC/Clang、sanitizers
13. P10.13 config: register Phase 10 long-step presets
      1800/900/450 s comparisonと30分超の上限探索。既存 Phase 5--8 presetは変更しない
14. P10.14 perf: meet the model-day speed and memory envelope
      profile、preconditioner改善、30日 explicit/semi-implicit比較
15. P10.15 validation: run the 1200-day pilot and close Phase 10
      長時間安定性、solver/retry統計、accuracy/cost、known gapsを記録
```

```text
entry -> P10.01 -> P10.02 -> P10.03 -> P10.04
                       |                    |
                       +-> P10.05 -> P10.06 +-> P10.07 -> P10.08 -> P10.09
                                                                  |
                                                                  v
                 P10.10 -> P10.11 -> P10.12 -> P10.13 -> P10.14 -> P10.15
```

依存の理由。

- P10.02 は numerical method変更前の最後の byte-exact 分割であり、以後の比較 oracleになる
- P10.03/04 は3-D couplingから独立して solver自体を falsifyする
- P10.06 の mode orderingが確定するまで automatic mode selectionを driverへ接続しない
- P10.07 の add-subtract identityと保存 gateが通るまで大 stepを許可しない
- P10.08 は external modeだけで最小 end-to-end pathを作り、P10.09 はその上へ internal modesを追加する
- P10.12 の terrain gateが通るまで climate presetを登録しない
- P10.14 で30日比較が accuracy/costの両方を通るまで1200日 runを開始しない

各コミットは build、対象 unit test、既存 Phase 0--9 regressionを通す。数値を変更するコミットに
baseline更新を混ぜず、旧 explicit結果と新 semi-implicit結果を別名で保持する。

## 9. 完了チェックリスト

P10.15 時点の状態。未達は1項目のみ。

### Method and compatibility

- [x] ADR 0016 が fast/slow split、reference state、implicit weight、mode選択、solver failureを固定している。
- [x] 既存 config、SSP-RK3 state、checkpoint、Frame、canonical text、fingerprintが byte-exactである。
- [x] `F=L_ref+(F-L_ref)` identityと `L_ref(U_ref)=0` が登録 testを持つ。
- [x] restartが追加履歴なしで uninterrupted runと一致する。

### Solver and conservation

- [x] GMRESが既知問題を tolerance内で解き、未収束・NaNを拒否する。
- [x] shallow-water fixtureが wave Courant 8まで安定し、centered設定で時間収束する。
- [x] vertical mode transformが往復し、必要 modeが再現可能に選ばれる。
- [x] implicit correctionが dry mass、`M*theta`、tracer massを後処理なしで保存する。
- [x] warm-up後の solver/RHS allocationが0である。

### Dynamics

- [x] flat linear waveが1800 sを受理し、wave Courant 3以上で振幅・位相 gateを満たす。
- [x] 2400/3600/5400/7200 s探索でcase別の最大合格stepと最初の不合格理由が記録される。
      登録 mode cap では 2400 s、cap を外すと 7200 s まで安定し、10800 s で explicit CFL に当たる。
- [x] isothermal restとDCMIP terrain restが既存の静止性を失わない。
- [x] linear mountain waveとUMJS14/JW06 baroclinic caseがstep細分化で収束する。
- [x] pressure、temperature、layer thickness、tracer、tangency違反をclipせず拒否する。
- [x] explicit CFL、implicit wave Courant、solver residual、retry理由が区別して記録される。

### Long step and performance

- [x] Held--Suarez 30日 runのmedian accepted stepが1200 s以上である。実測 1800 s、retry 0。
- [x] 同じモデル期間でexplicit 200 sより3倍以上速い。実測 5.78 倍。
- [x] GMRES p95 iterationが20以下、peak RSSがexplicitの2倍以下である。実測 1 反復、1.03 倍。
- [ ] 1本の1200日 semi-implicit runが完走し、iteration/retryの経時悪化がない。**未実施**。
      手順と合否条件は検証報告の P10.15 節にある。
- [x] GCC/Clang、ASan/UBSan、format-check、全 Phase 0--9 regressionが通る。
      Clang 開発スイート 85/85、GCC 15 と ASan/UBSan で phase10 ゲート 7/7。
- [x] `docs/validation/phase-10.md` に再現コマンド、比較値、失敗例、known gapsがある。

## 10. Phase 10 が既存の判断に与える影響

| 対象 | 影響 |
|---|---|
| ADR 0006 | 「semi-implicitには後続ADRが必要」を ADR 0016 で満たす。explicit SSP-RK3は残す |
| ADR 0011 | explicitのfast-wave CFL定義は維持。semi-implicitではadvective CFLとwave Courantを別診断にする |
| ADR 0014 | advective jump dissipationとLamb speedの分離をそのまま利用する |
| ADR 0015 | energy budget mitigationを維持し、時間積分差によるbudget変化を別 attributionする |
| Phase 5--8 preset | historical explicit baselineとして変更しない |
| checkpoint / Frame | prognostic layoutを変えないため変更なし |
| docs/plan.md | 従来の「Phase 10 物理拡張」を Phase 11 以後へ送り、Phase 10を本計画に固定する |

## 11. 参照

- D. Majewski et al. (2002),
  [The Operational Global Icosahedral--Hexagonal Gridpoint Model GME](https://doi.org/10.1175/1520-0493%282002%29130%3C0319%3ATOGIHG%3E2.0.CO%3B2):
  quasi-uniform全球格子、Eulerian dry dynamics、hybrid pressure、鉛直mode分解、external + 4 internal
  modesのsemi-implicit処理。
- J. Thuburn, C. J. Cotter, T. Dubos (2014),
  [A mimetic, semi-implicit, forward-in-time, finite volume shallow water model](https://doi.org/10.5194/gmd-7-909-2014):
  cubed-sphereを含む有限体積格子でのCrank--Nicolson fast-wave solve、Helmholtz反復、large wave Courant検証。
- S. Sandbach, J. Thuburn, D. Vassilev, M. G. Duda (2015),
  [A Semi-Implicit Version of the MPAS-Atmosphere Dynamical Core](https://doi.org/10.1175/MWR-D-15-0059.1):
  quasi-Newton Crank--Nicolson、geometric multigrid、緩いinner solve、terrain時のpressure-gradient注意点。
- ECMWF (2023),
  [IFS Documentation CY48R1 Part III](https://www.ecmwf.int/en/elibrary/81369-ifs-documentation-cy48r1-part-iii-dynamics-and-numerical-procedures):
  hydrostatic hybrid-coordinate equationsのwind divergence、temperature、surface pressureを結合する
  reference-linear semi-implicit定式化。
- ECMWF,
  [OpenIFS dynamical core](https://confluence.ecmwf.int/display/OIFS/3.1%2BOpenIFS%3A%2BDynamical%2BCore):
  hydrostatic semi-implicit/semi-Lagrangian構成とspectral Helmholtz solve。
- Met Office,
  [ENDGame](https://www.metoffice.gov.uk/research/foundation/dynamics/endgame):
  iterative semi-implicit couplingとHelmholtz問題を用いる実運用core。

文献の30--60分 stepはsemi-Lagrangian移流や異なる格子・演算子も含むため、その値をMyPlanetSimの保証には
使わない。Phase 10 の30分能力は上記の固有mode、Courant、accuracy、terrain、cost gateで独立に示す。

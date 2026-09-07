# Phase 12 — 乾燥対流調節・境界層の実装計画

状態: P12.01--P12.09 を実装し、独立 column・全球行列・1200日 pilot まで実施（2026-09-07）。
本書は計画のまま残す。ここに書いた設定名・係数・受け入れ閾値は計画時の提案値であり、
実際の測定値・成立範囲・未達項目は [Phase 12 検証報告](validation/phase-12.md) を参照する。

## 1. 目的と現状

Phase 11 の灰色放射に、乾燥不安定層の対流調節と、地表からの熱・運動量交換を鉛直に伝える
境界層を加える。最初に独立 column で保存・安定性を検証し、その後に全球 driver へ結合する。
湿潤対流、凝結、雲、予報 TKE、EDMF、完全な Monin–Obukhov 反復解法は後続とする。

確認した現状:

- [Phase 11 検証報告](validation/phase-11.md)の1200日 pilot は、200日以降の
  面積・層間等重みの対流不安定割合が平均 `0.09749`。日次最小 `d(theta)/dz` の最小値は
  `-0.010273 K/m`。対流調節を必要とする具体的な基準実験がある。
- `src/physics/gray_radiation_coupling.cpp` は放射に加え、
  `H=lambda*(Ts-Tbottom)` と `sigma>0.7` の Rayleigh drag を所有する。
  境界層を単純に足すと地表交換・摩擦を二重計上する。
- `src/dynamics/dry_hydrostatic_diffusion.cpp` は各 level の水平 diffusion。
  鉛直乱流混合の代用にはならない。
- state は `M*theta`、Cartesian 接ベクトルの水平運動量、受動 tracer、`Ts` を持つ。
  full-level Exner と静水圧 geopotential を再利用できる。
- driver は SSP-RK3 の stage と semi-implicit の反復内で RHS を評価する。
  瞬時の対流調節をこの RHS に副作用として入れることはできない。
- Phase 11 は局所放射保存を検証したが、全系の力学/時間積分残差の厳密分離は未完了。
  [ADR 0017](adr/0017-gray-radiation-column-and-energy-attribution.md)の有限上端・離散熱源
  attribution を引き継ぎ、対流導入だけでこの残差が解消したとは判断しない。

knowledge MCP では `MyPlanetSim` および `惑星|対流調節|Phase 11` を検索したが該当記録はなかった。
プロジェクトの現状判断は上記リポジトリ内のコード・ADR・検証報告による。

## 2. 既存モデルの調査と採用範囲

| 参照対象 | 確認できた実装・設計 | Phase 12 への反映 |
|---|---|---|
| climlab `ConvectiveAdjustment` | 指定 lapse rate を超える柱を瞬時・保存的に調節する。Akmaev 法を呼び、地表を含めるかを state で選択し、時間刻みに依存しない increment を返す | 大気のみの独立 adjustment kernel と、有限時間 tendency とは別の API |
| Isca surface flux / physics driver | 地表の熱・運動量交換を bulk 法で計算し、粗度・安定度を使う。対流 scheme と鉛直 diffusion を別に選択する | 放射、対流、地表交換・境界層の責務分離。地表の温度差には断熱補正を含める |
| GFDL Frierson–Held–Zurita-Gotor の理想化 GCM | bulk Richardson 数で境界層深さを決め、地表層につながり層頂でゼロになる diffusivity を使う。運動量・乾燥静的エネルギー・水蒸気を混合する | 安定度依存の地表交換と高さ依存 K を持つ簡易モデルの参照 |
| Isca mixed layer | 鉛直 diffusion と地表 reservoir を陰的に結合する | 薄い最下層・小熱容量でも顕熱を同一 flux で交換する column solve |

出典は [climlab の公式コード](https://climlab.readthedocs.io/en/latest/_modules/climlab/convection/convadj.html)、
[Isca 地表交換](https://execlim.github.io/Isca/modules/surface_flux.html)、
[Isca physics driver](https://execlim.github.io/Isca/latest/html/modules/idealised_moist_phys.html)、
[GFDL 著者公開論文のモデル記述](https://www.gfdl.noaa.gov/bibliography/related_files/dmwf0601.pdf)、
[Isca mixed layer](https://execlim.github.io/Isca/modules/mixedlayer.html)。2026-09-07 閲覧。

以下はこれらを参考にした MyPlanetSim の設計提案であり、特定モデルの完全移植ではない。
GFDL 論文の実験は独立した対流 scheme なしでも動作する構成で、境界層混合自体が深い対流を担う。
両方式の役割の重なりは本モデルでも比較する必要がある。初版の K profile は下記の簡略形とし、
同論文の地表層接続式・乾燥静的エネルギー拡散そのものとの一致は主張しない。
参照 fixture を作る段階で upstream の版/commit、入力、生成スクリプトを固定する。
既存コードを転記せず独立実装し、CI/runtime に Fortran や climlab の依存を追加しない。

## 3. スコープと設定の分離

推奨する初版は次の構成。

1. **乾燥対流調節:** 大気の温度のみを、中立な温位 profile へ保存的に調節する。
2. **境界層:** bulk 地表顕熱・応力、診断的な境界層深さ、熱・水平運動量・既存受動 tracer の
   鉛直拡散を扱う。地表との tracer flux はゼロ。
3. **結合:** 力学＋放射、境界層、対流調節の順で一つの試行 step を作り、全体を受理/棄却する。

`physics.kind=gray_radiation` を維持し、直交する `convection.kind` と `boundary_layer.kind` を追加する。
未指定はともに `none` とし、旧 config の canonical serialization・fingerprint・結果を維持する。
初版の本番組合せは gray 経路だけに限定する。独立 column test は放射なしでも呼べる。
Held–Suarez 等との組合せは silent ignore にせず、設定検証で未対応と伝える。

| 組合せ | 顕熱・下層摩擦の所有者 |
|---|---|
| 対流なし・境界層なし | 従来 gray adapter |
| 対流のみ | 従来 gray adapter。対流は大気の温度再配分のみ |
| 境界層あり、対流は任意 | 新境界層。gray adapter は放射・地表放射収支・内部熱のみ |

境界層ありで旧 `surface.air_exchange_coefficient_w_m2_k` の非ゼロを指定したら拒否する。
gray 経路内の既存 Rayleigh drag はこの場合に除外し、Held–Suarez の benchmark drag は変更しない。
新境界層では地表温度を対流混合ブロックに含めない。地表の熱容量・交換時間は境界層が扱う。

## 4. 乾燥対流調節

### 4.1 混合アルゴリズム

上から下へ `k=0..K-1` とし、固定した `ps`、層質量 `M=Delta p/g`、core の
`Pi=derived.exner_full` に対して `T=Pi*theta` を使う。乾燥安定条件は
`theta[k] >= theta[k+1]`。地球の `9.8 K/km` を固定値で埋め込まない。

不安定な隣接ブロックを併合し、その全層を次の一定温位へ更新する。

```text
w[k] = cp*M[k]*Pi[k]
theta_block = sum_block(w[k]*theta[k]) / sum_block(w[k])
T_new[k] = Pi[k]*theta_block
Delta(M*theta)[k] = M[k]*(theta_block-theta_old[k])
```

重み付き隣接ブロック併合法（stack を使う PAVA 型）なら、併合後に下流だけでなく直前ブロックとの
安定性も検査でき、`O(K)` で柱全体の安定性を得られる。二層だけの一回 sweep では終了しない。
安定柱・中立柱は恒等写像、再適用は冪等、不連続に離れた不安定領域は必要なときだけ併合する。
調節の境界判定には丸め誤差規模の温位 tolerance を使い、物理的な超断熱勾配を許す tuning にしない。

保存対象は固定質量での `E_h=sum(cp*M*T)`。単純な質量平均温位ではこの保存則にならない。
`M`、`ps`、水平運動量、tracer、`Ts` は直接変更しない。この限定的な対流調節は対流による
運動量・物質輸送を含まず、それらの完全な再現を主張しない。

### 4.2 全エネルギーとの関係

`E_h` 保存と既存の `E_dry=sum(M*(cv*T+Phi+|u|^2/2))` 保存を区別する。
調節前後で静水圧を再診断し、`Delta E_h` と `Delta E_dry` を両方記録する。
固定質量で `Phi` は温位に線形なので、ADR 0017 の熱源方向式を有限 increment に適用し、
実際のエネルギー差と一致することを検証する。境界仕事・上端・quadrature の寄与を導出して記録する。

`Delta E_dry` を全地球への一様熱源で打ち消さない。P12.01 でこの attribution を確定してから
調節 kernel を実装する。局所 enthalpy 保存を満たしても既存診断の差が残る場合は、それを数値結果として
公開し、全系の気候エネルギー閉鎖の完了条件を別に残す。

## 5. 境界層と地表交換

### 5.1 地表交換と K の初版

高さは海抜ではなく `z=(Phi-Phi_surface)/g`。地表 Exner `Pi_s=(ps/p0)^(R/cp)` と
`theta_s=Ts/Pi_s` を用いる。最下層の高さを `z_a`、状態を添字 `a` で表す。
静止地表に対する風を使い、顕熱 H は上向き正、`t_s` は大気に加わる地表応力とする。

```text
U_eff = sqrt(|u_a|^2 + U_g^2)
Ri_s = g*z_a*(theta_a-theta_s)/(theta_ref*U_eff^2)
f(Ri) = 1                     (Ri <= 0)
      = (1-Ri/Ri_crit)^2       (0 < Ri < Ri_crit)
      = 0                     (Ri >= Ri_crit)
C_D = kappa_vK^2 / log(z_a/z0m)^2 * f(Ri_s)
C_H = kappa_vK^2 / (log(z_a/z0m)*log(z_a/z0h)) * f(Ri_s)
H = rho_a*cp*C_H*U_eff*Pi_s*(theta_s-theta_a)
t_s = -rho_a*C_D*U_eff*u_a
```

`theta_ref` は地表・最下層温位の正の平均。`Pi_s` は温位差を地表圧力での温度差へ換算する
本計画の規約であり、Isca の式をそのまま写したものではない。
中立断熱柱と整合した地表温度なら H はゼロになる。粗度は `0<z0m,z0h<z_a` を検証する。
陸・海それぞれの粗度から exchange conductance を計算し、既存 land fraction で flux を面積平均する。
地表温度は既存どおり cell ごとに一つで、陸海二つの reservoir は追加しない。

`U_g` は未解像の風速に対する明示的パラメータ。初期試験値は `1 m/s` とし、`0/0.5/1` の感度を測る。
`U_g=0` かつ完全無風なら地表 flux をゼロと定義し、Ri の除算をしない。`U_g>0` でも
`u_a=0` の運動量 source はゼロ。非ゼロの gustiness が静止加熱柱の交換を可能にするという仮定を明記する。

境界層深さ h は最下層から上へ、最下層温位を基準とする bulk Ri が `Ri_crit` を初めて横切る
高さを補間して診断する。分母は各 level の `|u|^2+U_g^2`、高さは地表から測る。
交差なしは model top とし、その回数を出力する。完全無風の Ri 分類は浮力差の符号で定義する。
最下層より下の浅い安定境界層は未解像として記録し、架空の複数混合層を作らない。

初版の提案 K profile は `K_m=kappa_vK*u_star*z*(1-z/h)^2`（`0<z<h`）、
その外側でゼロ、`u_star=sqrt(C_D)*U_eff`、`K_h=K_m/Pr_t`、`K_q=K_h`。
粗度の異なる tile では各 K を計算して面積平均する。初期値は `kappa_vK=0.4`、`Pr_t=1`、
`Ri_crit=1`。これらは地球以外への普遍性を保証する値ではなく、`Ri_crit=0.25/0.5/1` も比較する。
背景 K は初版ゼロ、層頂での countergradient flux・entrainment mass flux は含めない。

### 5.2 共有界面 flux と離散化

内部界面 j の上側 U・下側 L に対し `dz=z_U-z_L>0` とする。
界面密度は column の圧力・温度から正値に再構成し、界面 Exner は core の half-level 値を使う。
熱に対しては乾燥中立柱の flux を厳密にゼロにする次の閉じ方を提案する。

```text
F_H[j] = -rho_j*cp*Pi_half[j]*K_h[j]*(theta_U-theta_L)/dz
F_m[j] = -rho_j*K_m[j]*(u_U-u_L)/dz
F_q[j] = -rho_j*K_q[j]*(q_U-q_L)/dz

cp*Pi[k]*d(M*theta)[k]/dt = F_H[k+1]-F_H[k] + D[k]
d(M*u)[k]/dt             = F_m[k+1]-F_m[k]
d(M*q)[k]/dt             = F_q[k+1]-F_q[k]
C_s*dTs/dt              = -H                  (境界層 substep)
```

上端は全 flux ゼロ。地表は `F_H[K]=H`、`F_m[K]=t_s`、`F_q[K]=0`。
同一の H を大気・地表に逆符号で使い、内部界面も一度だけ計算する。
これにより `sum(cp*M*Delta T)+C_s*Delta Ts=sum(D*dt)`、tracer 質量保存、
大気運動量変化＝地表応力 impulse が離散的に対応する。

この熱 flux divergence は `M*dtheta/dt=div(rho*K*grad(theta))` を直接解く方式とは異なる。
中立性と Phase 11 の熱源変換・enthalpy 収支をそろえるための選択であり、P12.01 で連続式との
関係と静水圧仕事を確認する。温位保存と enthalpy 保存を同時に保証すると記載しない。

### 5.3 陰解法と散逸熱

係数・幾何を substep 始点で評価して固定し、backward Euler の三重対角系で鉛直拡散を解く。
熱は大気 K 層と `Ts` の K+1 未知数を一緒に解く。運動量・tracer は K 未知数。
正の対角容量と非負の conductance を使い、M-matrix 性・tracer の最大値原理を検証する。
Cartesian 三成分は同じ係数で解き、接ベクトル性と柱の運動量収支を確認する。

静止地表の摩擦と内部 shear は運動エネルギーを減らす。初版はこの減少を大気の熱へ戻し、
物理的な shear/drag 散逸と backward Euler の数値散逸を別々に出力する。
固定係数では界面の `dt*G*|u_U_new-u_L_new|^2`、地表の
`dt*beta*|u_bottom_new|^2`、各層の `M*|Delta u|^2/2` の和が KE 減少に一致する。
界面分は隣接二層に半分ずつ、地表分は最下層、時間離散分は該当層へ配分し、D として熱 solve に加える。
これは局所の過程収支であり、既存力学の energy defect を補正する操作ではない。

係数固定の初版は一次精度。大きな dt での係数変化誤差は column の dt 細分化で測る。
必要なら同じ離散式を用いた局所 subcycling を追加し、無条件安定を精度保証と取り違えない。
Picard 反復・二次 IMEX は初版の必須条件にしない。

## 6. 時間積分・受理・再開

初版は次の一次 splitting を明示的に採用する。

```text
保存済み初期 state
  -> 既存 SSP-RK3 または semi-implicit による力学＋放射の候補
  -> 境界層・地表交換の陰的 step（有効時のみ）
  -> 乾燥対流調節（有効時のみ）
  -> 静水圧再診断・不変量検査
  -> state と全過程 budget をまとめて受理
```

- adjustment を RK stage / Newton 反復の途中に入れず、有限 increment として扱う。
- 初期出力・restart の読み込み時には黙って調節しない。初期の不安定度を記録し、最初の step から扱う。
- 境界層や調節が失敗したら、力学・放射を含む同じ初期 state から dt を半減して retry する。
  棄却分の熱量、impulse、調節回数を accepted budget に足さない。計算負荷 counter は別に残す。
- semi-implicit の既存方程式残差は「力学＋放射 substep」の残差として保存する。
  調節済み state をその方程式へ代入して solver 失敗と判定しない。
- 調節・混合後は RHS cache を必ず無効化する。時刻が同じでも state は変わる。
  最終 state の診断 flux と、各過程が実際に使った積算 flux を区別する。
- 放射 substep の budget は既存 RK/centered の重み、境界層は陰解法の flux 積分、
  対流は前後差で求める。放射だけの surface residual に境界層による `Delta Ts` を混ぜない。
- 既存 no-mixing 経路の三次/二次精度は維持する。新しい結合経路は一次で検証し、
  単に半 step を前後へ置いて二次 Strang と呼ばない。

予報 TKE・予報 h を追加しないので state/checkpoint payload は既存 surface layout を使える。
scheme・係数は fingerprint に含め、restart 後も同じ step 分割で再現する。
積算診断は Phase 11 の run segment 結合規約を拡張する。

## 7. 変更対象と出力

新規候補:

- `physics/dry_convective_adjustment.hpp/.cpp`: grid・I/O 非依存の柱調節。
- `physics/boundary_layer.hpp/.cpp`: bulk flux、Ri/h/K、陰的拡散・地表結合。
- `physics/dry_mixing_coupling.hpp/.cpp`: 3D adapter と過程別 increment/budget。

header は `include/myplanetsim/`、実装は `src/` に置く。workspace を再利用し、warm column/step の
heap allocation を抑える。現在の `gray_radiation` kernel の放射 flux 計算は再利用し、
coupling 側に顕熱・drag の所有者切替を設ける。
変更点は config/metadata、dry driver/workspace/diagnostics、CLI CSV dispatch、checkpoint 検証、
tests と CMake。Web Frame の意味や live UI は本 Phase で拡張しない。

主要設定案:

```text
physics.kind = gray_radiation
convection.kind = dry_adjustment
boundary_layer.kind = bulk_k_profile
boundary_layer.integrator = backward_euler
boundary_layer.critical_richardson = 1
boundary_layer.turbulent_prandtl = 1
boundary_layer.gustiness_m_s = 1
surface.air_exchange_coefficient_w_m2_k = 0
```

このほか陸海別 `roughness_momentum_m` / `roughness_heat_m` を追加する。
未知 enum、非有限値、非正の粗度・Pr・Ri 閾値、負の gustiness、不正な scheme 組合せを拒否する。
数値 tolerance・solver iteration limit は物理係数と分ける。

出力案:

- `convection_diagnostics.csv`: 調節前後の最小温位勾配、不安定割合、調節柱/層の割合、
  最大温度 increment、`Delta E_h`、`Delta E_dry`、過程 attribution 残差。
- `boundary_layer_diagnostics.csv`: H、地表応力、h、K の範囲、浅層未解像・model top 到達の割合、
  大気/地表の熱量、impulse、tracer 収支、物理/数値散逸と返還熱。
- `mixing_column.csv`: 選択柱の full-level T/theta/u/q と half-level flux/K、h、調節ブロック。
  level 種別を明示し、調節前後を対にする。
- retry 理由、受理 dt、kernel 呼出し/秒数、列 solve 数。瞬時 `W/m2`、
  受理区間 `J`・`N s`、面積平均と全球積算を区別する。

## 8. 実装順序と受け入れ条件

| 段階 | 作業 | 受け入れ条件 |
|---|---|---|
| P12.01 | ADR 0018、保存量・熱 flux・分割順・散逸熱の離散式、旧 baseline 固定 | enthalpy と E_dry の関係、境界仕事、gray の所有者切替を実装前に記録 |
| P12.02 | config と独立対流 kernel | 二層解析、不均一質量/Exner、多層併合、安定性、保存、冪等性 |
| P12.03 | 対流のみの driver 結合・過程 budget | 両 integrator の試行末尾で動作、retry/cache/restart、旧経路回帰 |
| P12.04 | 固定 K の鉛直拡散と地表 reservoir の陰解法 | diffusion 解析解、最大値原理、熱・運動量・tracer 収支、散逸熱の恒等式 |
| P12.05 | bulk flux、Ri/h/K と fractional surface | 中立・安定・不安定・無風極限、粗度検証、層頂 K=0、異なる惑星定数 |
| P12.06 | 境界層＋対流＋放射の driver 結合、CSV | 二重交換ゼロ、全過程 rollback、独立した残差、積算診断と restart 一致 |
| P12.07 | column 実験、地表付近を細分化した鉛直格子 | 加熱/夜間冷却/RCE、鉛直・時間感度、浅い境界層の解像範囲 |
| P12.08 | 全球行列・30日比較・性能 | 4通りの物理組合せ、陸海地形・同期回転、時間/格子/係数感度 |
| P12.09 | 登録条件の1200日 pilot、検証報告 | 長時間安定、局所収支、不安定度の前後、残差、solver 負荷の報告 |

主な tests は `test_dry_convective_adjustment.cpp`、`test_boundary_layer.cpp`、
`test_dry_mixing_core.cpp`、CLI/config/checkpoint 回帰。実装時に CMake へ登録する。
P12.03 で対流のみを先に評価できるため、境界層の効果を切り分けた状態で開発を進められる。

## 9. 検証行列と完了ゲート

### 9.1 独立 column

以下の数値は事前登録する提案閾値。変更時は旧結果と理由を検証報告へ残す。

1. **対流:** 二層解析、全柱不安定、離れた不安定領域、不均一 hybrid 質量、薄層、
   安定/中立柱、固定 seed のランダム柱。相対 enthalpy 誤差 `<=1e-12`、
   調節後 `theta[k]-theta[k+1]>=-1e-10 K`、二回目の最大変化 `<=1e-10 K`。
   単純な反復併合の独立比較器・条件をそろえた climlab fixture と照合する。
2. **熱 attribution:** 実際の静水圧再診断によるエネルギー差を解析式と照合する。
   差が十分大きいケースの相対誤差 `<=1e-6` と、丸め誤差用の絶対床を登録する。
3. **拡散:** 一定密度・一定 K・閉境界で cosine mode の減衰、線形 shear、受動 tracer。
   空間3解像度で次数 `>=1.7`、backward Euler は時間3刻みで `>=0.8`。
   非負 tracer と最大値原理、閉柱の運動量・tracer 保存を `1e-12` の相対誤差で検証。
4. **地表交換:** 一層＋地表の一定係数二熱容量系の解析的緩和、断熱中立で H=0、
   安定抑制、完全無風、純陸/純海と fraction 混合。対流を無効化して境界層単独を測る。
5. **散逸:** 内部応力相殺、地表 impulse、KE 単調減少、返還熱の一致。
   熱・地表・KE の局所収支残差は `max(1 J/m2, 交換量の絶対和)` で規格化し `<=1e-10`。
6. **過程結合:** 灰色放射＋地表加熱からの乾燥放射対流柱、夜間の安定境界層、薄い最下層、
   小熱容量。固定係数以外は dt/2、dt/4 で誤差減少を検証する。
   対流開始/停止の不連続点を滑らかな次数試験へ混ぜない。

独立 column の放射対流平衡は、一定外部入力に対し TOA と蓄熱が一致し、温度 drift と
調節後不安定度が十分小さいことを条件とする。日周期 forcing は周期平均で評価する。
対流 kernel が profile を中立にすることだけでは、放射対流平衡の成立とは判定しない。

### 9.2 鉛直格子

既存 uniform sigma の K=10 は接地境界層を解像する保証がない。既存 A/B 入力を使い、
地表付近を細かくした K=20/40/80 の column 用格子と K=20/40 の全球用格子を作る。
Earth-like 基準では最下層中心 30--100 m、1 km 以下に少なくとも5層を試験の開始目標とし、
実際の高さは温度・重力・地形ごとに診断して報告する。最浅の夜間境界層は未解像になり得る。
地球以外では scale height との比と、診断 h 内の実際の層数で評価する。
stretched preset の K を CLI/UI が uniform に再生成しない既存制約を維持する。

### 9.3 全球・長時間

同一初期状態・同一格子で「放射のみ」「放射＋対流」「放射＋境界層」「全部」の4通りを比較する。
これにより、境界層と対流の双方が不安定層へ及ぼす効果と、旧 drag を置換した効果を識別する。
調節前の不安定度も必ず測り、調節後ゼロという構造上当然の結果だけを改善の根拠にしない。

| preset 案 | 主な確認事項 |
|---|---|
| `phase12_gray_convection` | Phase 11 と同じ格子・旧地表交換で対流だけの差 |
| `phase12_gray_boundary_layer` | 新境界層だけの混合・摩擦・温度構造 |
| `phase12_gray_dry_mixing` | 全結合、地表付近を細分化した平坦基準 |
| `phase12_gray_land_ocean` | fractional 粗度/熱容量、地形、seam/corner |
| `phase12_gray_tidally_locked` | 昼側混合・夜側安定、gustiness 感度 |

短い CI は N=2/4・K=4/8 の接続試験。物理評価は N=6/12・K=20/40、上端1000/500/250 Pa。
1日比較は `dt=1800/900/450 s` と小刻み基準を使い、線形/非線形 solver 誤差と splitting 誤差を
別々に調べる。1800秒対450秒で温度の質量・面積重み RMS `<=0.5 K`、平均 TOA 差 `<=1 W/m2`、
平均 h の差 `<=10%` を仮の登録基準とする。平均 h がゼロのケースは相対差の対象外とする。
未達なら刻みを下げる。Phase 11 の30分能力を自動的に引き継がない。

30日は瞬時場の一致と、実際の時刻/flux 積算から計算する平均・変動幅を分ける。
係数感度は Ri・粗度・gustiness、数値感度は鉛直格子・dt・上端を一度に混ぜずに変更する。
登録条件で1200日 pilot を行い、200日以降の T/Ts/風/h/TOA/地表蓄熱・過程別収支を記録する。
乾燥質量 drift `<=1e-10`、非有限値・負の質量/温度/tracer ゼロ、column 保存ゲートを必須とする。
最初/最後100日の accepted dt、retry、solve 数、秒/モデル日も比較する。

### 9.4 完了条件

- [x] P12.01--P12.09 の実装・独立 column・結合・restart ゲートを通過した。
- [x] 対流が enthalpy を保存して安定性を満たし、E_dry との差を説明・計測できる。
- [x] 境界層の共有 flux、tracer 保存、地表 impulse、散逸熱が局所収支に対応する。
- [x] 旧顕熱・Rayleigh drag の二重計上がなく、無効時に旧結果・fingerprint が維持される。
- [x] 地表付近の解像範囲、結合の一次精度、30分刻み可否、係数感度が明記されている。
- [x] 全球行列・30日比較・1200日 pilot と全回帰、GCC/Clang、ASan/UBSan、format-check を実施した。
- [ ] `docs/validation/phase-12.md` に再現コマンド、失敗例と未解決事項を残した。
      版を固定した外部 climlab fixture との照合は未実施で、独立比較器は自前の反復併合に留まる。

物理過程の局所検証を完了することと、全球気候のエネルギー閉鎖・観測への conformance は別のゲート。
Phase 11 の未完了項目を閉じた扱いにはしない。全系残差の原因修正、二次の物理結合、非局所乱流、
湿潤対流の順序は本 Phase の測定結果をもとに次の計画で決める。

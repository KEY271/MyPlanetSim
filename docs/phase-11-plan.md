# Phase 11 実装計画 — 乾燥大気の灰色放射

**状態:** 実装中（2026-09-07）。Phase 11 を放射導入に割り当てる。方式・ゲートは
[ADR 0017](adr/0017-gray-radiation-column-and-energy-attribution.md) で固定し、段階ごとに実装・検証する。

## 1. 到達点

恒星からの短波、大気・地表からの長波、大気による吸収と下向き長波を計算し、各層の加熱率と
地表温度を、同一の界面放射 flux から更新する。短波・長波それぞれ一帯域の灰色近似
（two-band / semi-gray）を採用する。各水平 cell に独立な平面平行 1D column を置き、
水平放射輸送は扱わない。

Phase 8 の surface、orbit、fractional land と Phase 10 の semi-implicit driver を再利用する。
最初に column solver を解析解で検証し、その後に乾燥 3D core へ結合する。完了の意味は
「検証された放射輸送と熱収支を持つ乾燥モデル」であり、Earth climate や湿潤惑星の再現ではない。

| 範囲 | Phase 11 の扱い |
|---|---|
| 短波 | 恒星直達光の吸収、地表反射、反射光の大気による再吸収 |
| 長波 | 灰色二流近似による上下方向の吸収・熱放射、地表放射・反射 |
| 地表 | 既存の単一温度・面積熱容量、顕熱交換、内部熱 flux |
| 大気組成 | 指定する吸収係数で理想化。CO2 濃度からの自動換算はしない |
| 雲・水蒸気・散乱 | 対象外。Rayleigh 散乱、雲反射、潜熱、スペクトル線は後続 |
| 境界層・乾燥対流調節 | 対象外。既存の顕熱交換と下層 drag を利用する |
| 時間積分 | SSP-RK3 と既存 semi-implicit の両方。放射は fast-wave Jacobian に入れない |
| 出力 | C++ の診断 CSV と column profile。Web protocol 拡張は後続 |

放射だけでは対流不安定な温度勾配を解消できない。`d(theta)/dz < 0` の頻度・強度を記録し、
その状態を放射対流平衡と呼ばない。乾燥対流調節を必要とする気候実験は、別途その過程を検証して
から行う。数値安定性を保つための温度 clipping や隠れた Newtonian cooling は導入しない。

## 2. 現状と接続先

- `src/physics/surface_energy_balance.cpp` は地表の短波吸収、宇宙への直接長波放出、最下層との
  顕熱交換を持つ。大気の放射吸収・放出、下向き長波はまだない。
- `src/core/orbit_parameters.cpp` の恒星方向・距離依存 flux・昼夜判定を利用できる。
- `DryHydrostaticState` は既に `surface_temperature_k` を持つ。放射 flux は診断量なので新しい
  予報変数は不要。
- `src/dynamics/dry_hydrostatic_driver.cpp` は physics を RHS へ合成し、SSP-RK3 の各 stage と
  semi-implicit の各反復で評価する。surface budget と time-step limit もここへ接続済み。
- [ADR 0015](adr/0015-held-suarez-energy-budget-mitigation.md) にある乾燥力学のエネルギー残差は
  未解決。放射の flux 保存、大気への熱源変換、全系の数値残差をそれぞれ検証する必要がある。
- [Phase 10 検証報告](validation/phase-10.md) の 1800 s は Held--Suarez での結果であり、
  新しい放射・surface 結合で再登録する。水平依存 reference operator の改良は本 Phase に含めない。

新しい選択肢を `physics.kind = gray_radiation` とする。`none`、`held_suarez`、
`planetary_newtonian`、`surface_energy_balance` は意味・出力・既定動作を維持する。
新経路は放射 + 既存型の顕熱交換 + 下層 drag を所有し、旧 surface 放射を重ねて呼ばない。
Newtonian relaxation と灰色放射の同時指定は受け付けない。

## 3. 数値モデル

### 3.1 座標・光学的厚さ

half level を上から下へ `j=0..K`、層を `k=0..K-1` とする。`p_j=A_j+B_j*ps`、
`M_k=(p_{k+1}-p_k)/g` を各 RHS の現在の state から求める。光学的厚さは鉛直方向の値とし、
斜め入射・二流近似の係数をその定義へ混ぜない。

最初の opacity は温度・組成に依存しない指定値とする。

```text
Delta tau_SW,k = kappa_SW * (p_{k+1}-p_k) / g

kappa_LW(p) = kappa_LW,ref * (p/p_ref)^(n-1)
Delta tau_LW,k = kappa_LW,ref * p_ref / (g*n)
                 * [(p_{k+1}/p_ref)^n - (p_k/p_ref)^n]
```

`kappa` の単位は `m2 kg-1`、`p_ref` は放射用の固定基準圧力 Pa、`n>=1`。
`n=1` は一定質量吸収係数、`n=2` は圧力とともに増す簡易 opacity となる。
係数が一定でも `ps` と `g` によって柱の光学的厚さが変わる。地形上で `p/ps` に正規化して
全柱の厚さを人工的に同じにしない。実装では薄層の差分の桁落ちにも配慮する。

モデル上端を放射境界とし、その上の大気を透明として省略する。上端圧力を下げた感度試験を行い、
TOA と呼ぶ flux が有限 `p_top` に依存しない範囲を登録する。

### 3.2 短波

`mu0=max(0, r_hat dot s_hat)` とし、上端の直達下向き flux は `S(t)*mu0`。
昼側は次を上から下へ sweep する。

```text
SW_down[0] = S(t)*mu0
SW_down[k+1] = SW_down[k] * exp(-Delta tau_SW,k/mu0)
SW_up[K] = alpha_surface * SW_down[K]
SW_up[k] = SW_up[k+1] * exp(-D_SW*Delta tau_SW,k)
```

地表は拡散反射を仮定し、上向き反射光は二流の実効経路長 `D_SW` で減衰させる。
初期案は `D_SW=1.66`。この閉じ方は本計画のモデル選択であり、厳密な角度積分ではない。
夜側・`mu0=0` は全短波 flux をゼロにし、除算しない。小さな正の `mu0` を任意の floor へ
置換しない。夜側、terminator、透明層の極限を試験する。短波の大気熱放射はゼロ。

まず `kappa_SW=0` の透明短波から完成させ、その後に吸収を有効化する。大気散乱がないため、
`alpha_surface` は地表 albedo であり、観測された惑星 Bond albedo をそのまま代入する量ではない。

### 3.3 長波

層内温度一定・局所熱平衡を仮定し、`B_k=sigma_SB*T_k^4`、
`t_k=exp(-D_LW*Delta tau_LW,k)` と置く。初期案は `D_LW=1.66`。

```text
LW_down[0] = 0
LW_down[k+1] = t_k*LW_down[k] + (1-t_k)*B_k

LW_up[K] = epsilon_surface*sigma_SB*T_surface^4
           + (1-epsilon_surface)*LW_down[K]
LW_up[k] = t_k*LW_up[k+1] + (1-t_k)*B_k
```

地表は長波に不透明とし、放射率と吸収率を等しくする。`epsilon_surface<1` で吸収されない
下向き長波は上向きに反射させる。`1-t_k` は `-expm1(-D_LW*Delta tau_LW,k)` で計算し、
透明極限で精度を失わないようにする。厚い層では透過率のゼロへの underflow を正常な極限として扱う。

一往復の sweep で解けるので計算量は `O(K)`/column。層内一定の近似精度は鉛直細分化で測定し、
不十分なら後続の高次 source 積分として切り出す。

### 3.4 同一 flux から大気・地表を更新

正味 flux は上向き正とする。

```text
F[j] = LW_up[j] - LW_down[j] + SW_up[j] - SW_down[j]
Q_rad[k] = F[k+1] - F[k]                         [W m-2]
dT[k]/dt = Q_rad[k] / (cp*M_k)                   [K s-1]
d(M*theta)[k]/dt = Q_rad[k] / (cp*Pi_full[k])

H = lambda_air*(T_surface - T_air,bottom)         [W m-2, 上向き正]
C_surface*dT_surface/dt = -F[K] + Q_internal - H
d(M*theta)[bottom]/dt += H / (cp*Pi_full[bottom])
```

`Pi_full` は core が温度診断に使用する `derived.exner_full` を利用する。
`pow(p_full/p0,kappa)` で別の Exner を再構成しない。放射・顕熱は乾燥質量、トレーサ質量、
地表気圧に直接 source を加えない。

柱ごとに `sum(Q_rad)-F[K]=-F[0]` が telescoping で成立する。
顕熱も大気・地表で相殺する。地表の長波放出と TOA の OLR は異なる量として出力する。

この式による `cp*M*dT` の収支と、実装の `sum M*(cv*T+Phi+|u|^2/2)` の変化は
自動的に同一とは限らない。実装開始時に、静水圧再診断・有限上端の境界仕事・鉛直 quadrature を
含めて関係を導出する。既存の Held--Suarez の離散熱エネルギー attribution を参考に、固定質量で
熱源方向の全エネルギー微分を独立検算し、その差を明示する。全系の残差を放射 flux に戻す
global energy fixer は使わない。

## 4. 時間積分と安定性

初版はすべての RHS 呼出しで放射を再評価し、放射呼出し間隔の間引きや flux の時間補間を入れない。
SSP-RK3 は stage の時刻、semi-implicit は初期時刻と候補時刻の orbit と温度を使用する。
driver の `rhs_is_time_independent` 判定にも新 physics を追加し、最初の反復で日射を誤って使い回さない。

放射・顕熱を含む `K+1` 個の温度の source Jacobian を考え、熱容量で割った各行の絶対値和から
保守的な rate bound `L=max_i sum_j |d(Tdot_i)/dT_j|` を作る。
固定 opacity の flux sweep を利用した `O(K)` の bound 評価を設計し、小さな柱の有限差分
Jacobian と固有値で検証する。新経路では `dt_physics <= radiation.cfl/L` を既存 CFL の最小値へ加える。
初期案は `radiation.cfl=0.5`、全 rate がゼロなら制約は無限大。

特に地表だけの `C/(4*epsilon*sigma_SB*T_s^3)` では、薄い大気最下層との顕熱交換や大気の放射冷却を
制限できない。大気・地表の双方を含む bound を使い、state bounds、stage 制約、候補 state を検査する。
違反時は step を半減して最初から retry し、設定済み最小刻み未満で理由付き停止とする。

放射は fast-wave Jacobian の対象外なので、semi-implicit の固定反復数 3 回を無条件には流用しない。
`dt=1800/900/450 s` と反復数 `3/4/5` を比較して、反復誤差が時間離散化誤差を下回る組を登録する。
source-only の SSP-RK3 は三次、十分な反復の centered semi-implicit は二次を目標とする。
時間依存の日射では滑らかな区間と昼夜切替を分けて次数を評価する。

accepted step の budget は実際の時間積分の重みを使う。SSP-RK3 は `1/6,1/6,2/3`、
semi-implicit は既存の `1-w,w`。固定反復後の方程式残差に由来する差も記録する。
棄却した attempt は積算熱量・気候統計に入れず、retry の理由別累積数だけを残す。

## 5. API・設定・出力

新規候補は `include/myplanetsim/physics/gray_radiation.hpp` と
`src/physics/gray_radiation.cpp`。column kernel の入力は half-level 圧力、温度、地表温度、
`g`、`S`、`mu0`、opacity と地表係数。出力は SW/LW の上下界面 flux と層別 flux convergence。
力学 state・格子・I/O に依存させず、workspace を再利用する。

3D adapter が各 cell の column を呼び、温位 tendency、surface tendency、drag と診断を合成する。
既存 surface の顕熱・drag の共通化は、旧経路の結果が変わらない範囲に限定する。

設定キー案:

```text
physics.kind = gray_radiation
radiation.shortwave_absorption_m2_kg = 0
radiation.longwave_absorption_ref_m2_kg = 0.0003
radiation.reference_pressure_pa = 100000
radiation.longwave_pressure_exponent = 1
radiation.longwave_diffusivity_factor = 1.66
radiation.shortwave_diffuse_factor = 1.66
radiation.cfl = 0.5
```

これに既存 `orbit.*`、`surface.*`、`planet.*`、hybrid grid、time integrator 設定を組み合わせる。
上の数値は実装時の試験開始値であり、観測に合わせた opacity ではない。`g=9.81`、
`ps-p_top≈1e5 Pa` なら長波の鉛直光学的厚さは約 3.06 となる。
係数の非負・有限性、指数、正の基準圧力、albedo/emissivity の `[0,1]`、熱容量の正値を検査する。
新キーは gray 経路の canonical config と fingerprint に入れ、旧 config の canonical 表現を変えない。

- `radiation_diagnostics.csv`: TOA 入射 SW、反射 SW、OLR、正味 TOA flux、surface の上下 SW/LW、
  大気 SW/LW 吸収、surface 蓄熱、顕熱、内部熱、界面収支残差。平均 `W m-2` と積算 `J` を区別する。
- `radiation_column.csv`: 選択 column の時刻、cell ID、half-level 圧力・光学的厚さ・界面 flux、
  full-level 温度・加熱率。half/full level を曖昧に混ぜない schema を定義する。
- 全系診断: 放射入力、離散大気熱源寄与、drag work、diffusion、力学残差、時間積分残差を分ける。
  TOA の物理的な非平衡と数値的な収支誤差も分ける。
- 性能・安定性: 放射秒数、呼出し数、physics 制限刻み、実際の刻み、制約別・retry 理由別累積数。
  sampled CSV からだけ retry 回数を推測しない。

flux は現在の state と時刻から再計算できるため checkpoint payload は既存
`dry_hydrostatic_surface_v1` を再利用する。config 検証、surface 初期化、checkpoint dispatch、
CLI の CSV dispatch にある physics kind の分岐をすべて点検する。restart では orbit の原点を
引き継ぎ、途中時刻を新たな日射原点にしない。積算診断・統計は保存または run segment の
明示的結合で連続 run と比較する。FrameV2 の意味は変更しない。

## 6. 実装順序と受け入れ条件

| 段階 | 作業・主な変更対象 | 受け入れ条件 |
|---|---|---|
| P11.01 | ADR 0017、符号・opacity・離散エネルギー関係、旧 baseline 固定 | 完了: 境界条件、source 変換、有限上端、比較条件を実装前に記録 |
| P11.02 | config/metadata、gray column API、opacity | 不正設定拒否、薄層・`ps/g` 感度・全柱光学厚さ、旧 canonical identity |
| P11.03 | 長波 sweep、透明短波、地表境界 | 等温柱の解析 flux、透明・厚い極限、`epsilon<1` の反射と保存 |
| P11.04 | 短波吸収・地表反射の再吸収 | Beer--Lambert、昼夜・terminator、柱ごとの入射＝反射＋吸収 |
| P11.05 | 独立 column の熱源・surface 結合、安定刻み | 一層平衡、過渡冷却、顕熱相殺、source-only 時間収束、retry |
| P11.06 | dry driver/workspace、SSP-RK3 と semi-implicit | 同時刻 RHS、質量非干渉、積分重み、反復数・刻み感度、旧経路回帰 |
| P11.07 | CLI/CSV、checkpoint、統計、理由別 counter | 再開時の温度・flux・orbit・積算値一致、棄却分の二重計上なし |
| P11.08 | 全球 preset と比較・性能測定 | 透明大気回帰、地形・陸海・同期回転、水平/鉛直/時間解像度感度 |
| P11.09 | 長時間 pilot、検証報告、README 更新 | 登録した数値ゲートの結果・再現コマンド・未解決範囲を記録 |

各段階に必要な tests を追加する。主な対象は `tests/test_gray_radiation.cpp`、
`tests/test_radiative_column.cpp`、`tests/test_radiative_dry_core.cpp` と config/CSV/checkpoint 回帰。
実装時に CMake の登録も行う。

## 7. 検証行列と完了ゲート

以下は提案する事前登録値。基準を変更した場合は理由と変更前の結果を検証報告に残す。

### 7.1 column の厳密な検証

1. **輸送のみ:** 等温柱の `B+(F_in-B)*exp(-D*tau)`、直接短波の指数減衰、ゼロ opacity、
   厚い層、長波の非黒体地表。解析 flux 誤差は `1e-10 W m-2 + 1e-12*F_scale` 以下。
2. **界面保存:** `|sum(Q_rad)-F[K]+F[0]|/max(1 W m-2, sum(abs(boundary fluxes))) <= 1e-12`。
   SW/LW 別、地表反射込み、ランダムな妥当柱でも確認する。
3. **一層温室平衡:** 透明 SW、黒体地表、顕熱・内部熱ゼロ、一定吸収短波 `Q`。
   大気の長波吸収率を `0<a<=1` とすると `T_a^4=T_s^4/2`、
   `sigma_SB*T_s^4=Q/(1-a/2)`。平衡解の RHS と、摂動からの緩和を確認する。
   `a=0` の大気温度は平衡条件で拘束されないので別ケースにする。
4. **過渡・時間精度:** 解析的な地表単独冷却、顕熱交換のみ、滑らかな時変日射。
   3 刻み以上で SSP-RK3 は観測次数 `>=2.7`、centered semi-implicit は `>=1.7`。
5. **鉛直精度:** 固定した滑らかな `T(p)` を `K=10/20/40/80` で比較する。高精度積分または
   独立した細分化 column と TOA/surface flux・加熱率を照合し、誤差減少と観測次数を報告する。
6. **source 変換:** 使用する Exner と温度変化を照合する。全エネルギーの熱源方向微分は中心差分の
   perturbation ladder と比較し、上端・quadrature 寄与を含めて相対 `1e-6` を目標とする。

輸送の独立比較器として climlab の GreyGas も利用できる。使う場合は版、入力、生成手順を固定した
小さな fixture にし、同じ近似の比較であることを明記する。CI/runtime 依存には追加しない。

### 7.2 全球結合

| preset 案 | 目的・比較 |
|---|---|
| `phase11_transparent_surface` | 両 opacity ゼロで Phase 8 の surface tendency と短時間解へ一致。新しい安全刻みによる差は同じ刻みで比較 |
| `phase11_gray_uniform` | 平坦・一様地表、透明 SW、灰色 LW。放射結合・質量・解像度の基準 |
| `phase11_gray_shortwave` | 同じ条件で SW opacity だけ変更し、大気短波加熱と地表吸収の再配分を検証 |
| `phase11_gray_land_ocean` | fractional 熱容量、地形上の `ps`、海岸・seam/corner を検証 |
| `phase11_gray_tidally_locked` | 固定恒星方向、夜側 SW ゼロ、昼夜の蓄熱・輸送・放射を検証 |

短い CI は `N=2/4,K=4/8`。bounded comparison は `N=6/12,K=10/20`、追加の鉛直感度は
`K=40`。smooth な短時間解で空間・時間を別々に細分化する。
semi-implicit は 1800/900/450 s と十分小さい explicit の双方へ比較し、
初期 1 モデル日の 1800 s 対 450 s の目安を温度 RMS 差 `<=0.5 K`、
全球平均 TOA flux 差 `<=1 W m-2` とする。30 日は瞬時場の一致だけでなく平均・変動幅を比較する。
この範囲を外れたら刻みを下げ、1800 s 対応とは登録しない。

登録した平坦基準の 1200 日 pilot を 1 本実行し、200 日以降の温度・風・TOA flux・地表蓄熱・
対流不安定度・力学残差を記録する。乾燥質量の相対 drift `<=1e-10`、非有限値・負の層質量ゼロ、
放射界面保存を必須とする。最初と最後の 100 日の solver work、retry、accepted step の分布も比較する。
TOA imbalance は蓄熱とともに報告し、1200 日の完走だけで統計平衡を宣言しない。

30 分能力は平坦基準の median accepted step `>=1200 s` と上記精度で登録する。
薄い大気・小熱容量で physics CFL が厳しい場合は、短縮を正しい結果として別に記録する。
性能は同じ radiative case の explicit/semi-implicit を比較し、1 モデル日あたりの時間、
放射の割合、メモリ、初回以降の allocation を測定する。Phase 10 の 5.78 倍を流用しない。

### 7.3 完了条件

- [ ] P11.01--P11.09 と column の解析・保存・時間収束ゲートを通過した。
- [ ] 大気・地表は同一界面 flux を使用し、source 変換誤差と既存力学の energy defect を区別できる。
- [ ] 指定した全球行列と 1200 日 pilot が終了し、刻み・精度・性能の適用範囲が明記されている。
- [ ] restart と積算診断が連続 run と一致し、旧 config/checkpoint/Frame に回帰がない。
- [ ] C++ 全回帰、GCC/Clang、ASan/UBSan、format-check が通る。
- [ ] `docs/validation/phase-11.md` に再現手順、定量結果、失敗例、既知の制約がある。

全系 energy defect が残る場合、放射 scheme の局所保存と気候のエネルギー閉鎖を混同しない。
climate conformance は未完了として残し、原因となる core の修正を後続課題にする。

## 8. 参考と後続の拡張

- [GFDL Idealized Moist Spectral Atmospheric Model](https://www.gfdl.noaa.gov/idealized-moist-spectral-atmospheric-model-quickstart/)
  は指定した光学的厚さを用いる灰色放射の GCM への利用例。ここでは乾燥の最小構成を提案しており、
  同モデルの湿潤物理や気候への conformance は主張しない。
- [climlab GreyGas の公式実装](https://climlab.readthedocs.io/en/latest/_modules/climlab/radiation/greygas.html)
  は上下 flux、層放射、地表反射、flux convergence の比較対象。本計画の pressure opacity と
  直接光/拡散光の角度近似は、別途ここで指定した設計である。

次の物理拡張は、まず保存的な乾燥対流調節・境界層を独立 column で検証し、その後に水蒸気・凝結・
雲と多帯域/相関 k 放射を検討する。放射の呼出し間引きや陰的 source solver は、Phase 11 の
性能・stiffness の実測で必要性を確認してから設計する。

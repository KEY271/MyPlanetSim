# Phase 13 — 蒸発・水蒸気・降水の実装計画

状態: P13.01--P13.08 を実装済み（2026-09-08）。全球比較・長期検証は未実施。
以下の全球実験用係数・受け入れ閾値は提案値であり、実測値ではない。

## 1. 目的と設計判断

Phase 12 の乾燥大気・灰色放射・境界層に、輸送される水蒸気、地表蒸発、潜熱、即時降水、
簡易湿潤対流、陸面 bucket を加える。ユーザー提案の希薄水近似、複数 tracer、
operator splitting、instantaneous saturation adjustment、Simple Betts–Miller / Frierson 型対流、
bucket model を採用する。

最初の到達点は水循環と潜熱による循環変化を過程収支つきで検証できる理想化 GCM とする。
雲量・雲水・雨水の予報、降水の落下時間・再蒸発、氷・雪・融解、植生・地下水・河川輸送、
水蒸気・雲による放射吸収の変化は後続 Phase に分ける。

| 項目 | 初版の選択 |
|---|---|
| 水の力学的寄与 | 乾燥空気の質量、`ps`、`R_d`、`cp`、`cv` を変えない。仮温度補正も使わない |
| 熱力学的寄与 | 一定の蒸発潜熱 `L_v` による地表冷却・大気加熱を含める |
| 大気の水相 | 予報は水蒸気のみ。液相への変換量は同じ物理 substep で地表へ渡す |
| 格子スケール凝結 | 温度上昇と飽和量の変化を同時に解く瞬時飽和調節 |
| 対流 | 有限時間の SBM 型。深い降水対流と非降水の浅い対流を実装する |
| 地表 | 既存の cell ごとに一つの `Ts` と熱容量を維持。陸の水だけ bucket で予報 |
| 海 | 蒸発に水量制限のない reservoir。水の受け払いを積算し、全系の水収支に含める |
| 放射 | Phase 11 の既定灰色吸収係数を使い、水蒸気量には依存させない |

仮温度を省くことは比熱・質量を固定する条件から自動的に決まるものではなく、追加の簡略化である。
EOS・静水圧・境界層 Ri・parcel 浮力で一貫して省き、潜熱による温度変化はすべてに反映する。
将来仮温度を導入する際はこれらをまとめて変更する。初版は地球型の希薄水蒸気の範囲を対象とし、
蒸気が大気の主要成分になる惑星への適用は主張しない。

## 2. 現状と参考モデル

確認したコード上の前提:

- `DryHydrostaticState` は `M*theta`、水平運動量、単一の `M*q`、任意の `Ts` を持つ。
  `Derived`、face reconstruction/flux、水平拡散、鉛直結合、RK 合成、semi-implicit 残差、
  workspace、診断、checkpoint が単一 tracer を前提とする。
- checkpoint は `dry_hydrostatic_cell_column_v1` / `dry_hydrostatic_surface_v1`。
  大気 payload の長さは `cells + 5*cells*levels` に固定されている。
- `dry_hydrostatic_diffusion.cpp` の tracer 拡散は `M*Laplacian(q)`。
  不均一な層質量での水の保存と非負性を新しい経路で検証・確保する必要がある。
- Phase 12 の境界層は tracer 鉛直拡散を陰的に解くが地表 tracer flux はゼロ。
  陸海別粗度と面積平均、地表との熱交換、運動エネルギー散逸の返還熱を再利用できる。
- 両時間積分器に「力学＋放射 → 境界層 → 乾燥対流 → 検査 → 受理」と全体 retry がある。
  この境界を湿潤物理に拡張し、凝結を RK stage や semi-implicit 反復内へ入れない。
- [Phase 12 検証報告](validation/phase-12.md)では地表細分化 K=20/40/80 を検証済み。
  浅い境界層では内部界面の K がすべてゼロになる場合があり、湿潤対流との役割を比較する。
  既存の全系エネルギー残差の厳密分離は未完了である。

knowledge MCP を `MyPlanetSim` で検索したが該当記録はなかった。
現状判断はリポジトリのコード・ADR・検証報告による。

| 参照 | 確認した内容と採用範囲 |
|---|---|
| Frierson–Held–Zurita-Gotor (2006) | 予報水蒸気を放射から切り離した灰色大気、slab surface、bulk 蒸発、凝結による潜熱というモデル階層を採用 |
| Isca SBM 公式説明 | CAPE による起動、参照温湿度への緩和、enthalpy 補正、深い／浅い対流の分岐を採用。開始値は `tau=7200 s`、`RH_ref=0.8` |
| Isca `qe_moist_convection.F90` | Frierson 型の shallower shallow convection を参照。仮温度・厳密な湿度変換を含む upstream と、初版の希薄近似の差を明記する |
| Isca large-scale condensation | 飽和閾値 1、液体のみ、落下中の再蒸発なしの構成を参考にする |
| Isca surface flux / mixed layer | bulk 水蒸気交換、bucket の乾燥度による蒸発制限、地表潜熱の陰的結合を参考にする |

一次資料（2026-09-07 閲覧）:

- [Frierson et al. (2006), A Gray-Radiation Aquaplanet Moist GCM, Part I](https://www.gfdl.noaa.gov/bibliography/related_files/dmwf0601.pdf)
- [Isca: Simple Betts–Miller](https://execlim.github.io/Isca/modules/convection_simple_betts_miller.html)
- [Isca: SBM ソース](https://github.com/ExeClim/Isca/blob/master/src/atmos_param/qe_moist_convection/qe_moist_convection.F90)
- [Isca: Large Scale Condensation](https://execlim.github.io/Isca/modules/lscale_cond.html)
- [Isca: Surface flux](https://execlim.github.io/Isca/modules/surface_flux.html)
- [Isca: Mixed layer](https://execlim.github.io/Isca/modules/mixedlayer.html)

以下の離散式は MyPlanetSim への設計提案であり、Isca と同一の実装・気候を保証しない。
P13.01 で比較対象の commit と熱力学設定を固定し、差異の表と独立 column fixture を用意する。
既存コードの転記ではなく独立実装とし、runtime/CI に Fortran や外部 GCM の依存を追加しない。

## 3. 複数 tracer と保存的輸送

### 3.1 状態と識別

`Ntr` 個の registry と連続配列を導入する。各 tracer は一意の名前、役割
`passive | water_vapor`、初期条件、非負性要件、水平拡散の選択を持つ。
`water_vapor` は湿潤時にちょうど一つ必要。物理 kernel は registry が解決した index を受け取り、
配列の 0 番を水と仮定しない。混合比の和を 1 に正規化しない。

保存配列の提案は tracer-major の `Q[(tracer*cells+cell)*levels+level]`。
既存単一 tracer の cell-column 配置を保ち、柱を `span` として取り出せる。
`q=Q/M` はすべて乾燥キャリア質量あたりの量。水についても本書の `q_v` は
`kg_water/kg_dry_air` と定義し、厳密な specific humidity と混同しない。

一般化の範囲は 3D core と、それが使用する鉛直輸送・境界層・時間積分・I/O。
Phase 1/2 の独立 benchmark や Phase 4 の公開単一 column state を一括改名せず、
必要な低水準輸送 kernel を成分ごとに再利用する。既存 `DryHydrostatic*` 名も本 Phase で全面改名しない。

### 3.2 輸送・拡散・非負性

- 各成分は同じ乾燥質量輸送と整合させ、定数 `q_i` の保持、共有 face の反対称性、
  水平＋鉛直の成分別保存を検証する。tracer によって `ps` や fast mode を変えない。
- 移流の再構成・limiter は成分別に適用する。ゼロ場の short circuit も成分別にする。
  semi-implicit の候補・残差・補正・合成で全成分を更新する。solver の受理判定には全成分を含める。
- 湿潤経路の tracer 水平拡散は共有界面の `F_q=-rho*K*grad(q)` の有限体積形とし、
  各隣接 cell の質量増分を逆符号で更新する。傾いた hybrid 面の係数・計量は core に合わせて定義する。
  保存的な二次拡散を初版とし、水蒸気には非負性を保証できない biharmonic 拡散を許可しない。
  速度・温位の既存拡散は独立して選択できるようにする。
- 非負性は輸送 flux の制限または移流 subcycling で守る。semi-implicit の centered 更新についても
  SSP-RK3 と同じ保証があるとは仮定せず、湿潤柱・鋭い乾湿境界で検証する。
  dt 半減だけでは負値が解消しない場合は、受理済み乾燥質量 flux を用いる保存的 tracer substep を実装する。
- `max(q,0)` による無記録の水の生成はしない。丸め誤差の負値を補正する場合も、同じ成分内で
  保存的に再配分し、補正量と回数を記録する。供給元のない全ゼロ場は厳密にゼロのままとする。
- 境界層は同じ `K_q=K_h` で全成分を混合する。地表 flux は水蒸気にのみ与える。
  SBM は温度・水蒸気のみを変え、受動 tracer・運動量の対流輸送は初版の対象外と明記する。

単一の旧 tracer 経路は既存設定の結果を維持する。新経路の成分順序を入れ替えても、名前で対応した結果は
一致し、受動 tracer の追加で力学・水蒸気の結果が変わらないことを受け入れ条件とする。

## 4. 共通湿潤熱力学と瞬時飽和調節

### 4.1 熱力学の定義

初版は一定潜熱・液相のみの Clausius–Clapeyron を共通 kernel に置く。
蒸発、凝結、相対湿度、LCL、parcel ascent は同じ関数と定数を使用する。

```text
e_s(T) = e0 * exp[(L_v/R_v)*(1/T0 - 1/T)]
q_s(T,p) = epsilon*e_s(T)/p       epsilon=R_d/R_v
dq_s/dT = q_s*L_v/(R_v*T^2)
RH = q_v/q_s(T,p)
```

提案値は `T0=273.16 K`、`e0=611.657 Pa`、`L_v=2.5e6 J/kg`、`R_v=461.5 J/(kg K)`。
`R_d` と `cp` は惑星設定から取る。`q_s` は希薄近似で、厳密な混合比の分母 `p-e_s` と混用しない。
氷点下でも理想化した液相の式を継続するため、雪・氷の物性は再現しない。
温度・圧力の有限性と正値を検査し、`q_v` と `e_s/p` の最大値を出力する。
開始時の希薄性上限案はそれぞれ `0.1`。超過を clipping で隠さず、入力不適合または適用範囲外として
報告する。数値 solver の不収束とは区別する。

### 4.2 格子スケール凝結

固定した `p`、`M=Delta p/g`、`Pi` で、`q<=q_s(T,p)` なら恒等写像。
過飽和なら凝結混合比 `c` を次の一変数方程式で解く。

```text
q_new = q_old - c
T_new = T_old + (L_v/cp)*c
q_old - c = q_s(T_old + (L_v/cp)*c, p),   0 <= c <= q_old

Delta(M*theta) = M*(L_v/cp)*c/Pi
rain_lsc_amount = sum_k M[k]*c[k]        [kg/m2/呼出し]
```

関数は単調なので bracket 付き Newton/bisection を使う。温度を固定した
`q=q_s(T_old,p)` の clipping では代用しない。`cp*T+L_v*q` の保存、飽和への到達、
非負性、再適用の冪等性を同時に検証する。
通常解の停止条件は RH 残差と質量・enthalpy 残差の両方で定義し、極低湿度では絶対誤差も併用する。

凝結水は全層から同じ cell の地表へ渡す。落下中の再蒸発・水の heat capacity・落下の
位置／運動エネルギーは含めない。凝結潜熱を地表へもう一度与えない。
SBM の降水量と格子スケール降水量は別々に積算する。

## 5. 蒸発・露・陸面 bucket

### 5.1 水 flux と reservoir

`E` は地表から大気へ上向き正の水 mass flux、`H` は上向き顕熱、`P` は下向き降水。
最下層値を `a`、陸海 tile を `t` とし、Phase 12 の `rho_a`、`U_eff`、安定度・粗度を再利用する。
初版は水蒸気交換係数 `C_E,t=C_H,t`、水蒸気粗度は熱粗度と同じとする。

```text
E_pot,t = rho_a*C_E,t*U_eff*(q_s(Ts,ps)-q_a)
beta(W) = min(1, W/(0.75*Wmax))

E_ocean = E_pot,ocean
E_land  = beta(W)*E_pot,land     (E_pot,land >= 0)
        = E_pot,land             (E_pot,land < 0; 露)
E_cell  = f_land*E_land + (1-f_land)*E_ocean
LE_cell = L_v*E_cell
```

露を乾燥度で抑制しない。地表に水がない場合でも大気から露が供給され、`LE<0` により地表を暖める。
最下層から抜く水量には非負性の制約を適用する。塩分による海面飽和量の低下は初版で省く。

`W` は **陸面積あたり** の液体水量 `[kg/m2_land]`。開始値案は `Wmax=150 kg/m2`、
`W_initial=0.5*Wmax`。水深換算ではそれぞれ約 150 mm、75 mm。
`f_land=0` の cell は bucket をゼロに保ち、`f_land` で割る操作を行わない。
同一 cell の降水強度は陸海で等しいと仮定する。

```text
dW/dt = P - E_land - runoff_land
0 <= W <= Wmax
```

降水は全量浸透させ、容量超過だけを流出とする。地下への排水時定数は追加しない。
有限 step では実際に交換した水量で W を更新する。まず蒸発・露で更新し、SBM/凝結の後に降水を加え、
その都度容量超過を流出へ渡す。同じ substep の末尾に降った水は次の substep から蒸発に利用する。
`dt*E_land<=W_start` を守る制約も surface solve に組み込み、事後 clipping だけで済ませない。

海は枯渇を計算しない。海への降水、海からの蒸発、陸からの流出を受ける診断的な累積海洋水量差
`Delta W_ocean` を保持する。河川の経路や遅延はなく、流出は全球海洋 reservoir へ即時移す。
全陸実験では流出を系外流出として別台帳へ記録し、存在しない海へ分配しない。

### 5.2 陰的境界層・地表結合

既存の顕熱・運動量・鉛直 tracer solve を拡張し、同じ受理済み水 flux を最下層、地表熱、bucket、
海洋台帳で使用する。

```text
sum_k Delta(M*q_v) = dt*E_cell              (境界層 substep)
C_s*Delta Ts = -dt*(H + L_v*E_cell)
sum_k cp*M*Delta T = dt*H + returned_dissipation_heat
```

輸送された水蒸気の潜熱は `L_v*M*q_v` に蓄える。蒸発時に大気へ `L_v*E` の顕熱を直接加えない。

幾何・Ri・K・bulk 係数は各局所 substep の始点で固定する。`q_s(Ts_new,ps)` と水量制限は
新しい `Ts`、`q_a`、`W` と結合させる。柱の三重対角 solve を再利用し、surface の少数未知数を
bracket/Picard と active-set（蒸発／露、bucket 空／部分充填／満水）で解く設計とする。
収束後に全過程を同じ flux で再評価する。水量だけを制限して潜熱 flux を据え置くことは許さない。

初版から物理 subcycling を用意し、係数固定と splitting の精度を管理する。
開始値案は `physics_substep<=300 s`、さらに `<=tau_conv/4`。局所 solve 失敗時は substep を細分化し、
下限でも失敗した場合は力学を含む試行全体を棄却する。陰解法であることを大刻みの精度保証にしない。

## 6. Simple Betts–Miller / Frierson 型対流

### 6.1 Parcel と参照場

最下層から parcel を持ち上げ、LCL までは乾燥断熱、その上では共通の希薄・一定潜熱モデルの
pseudoadiabat を積分する。最下層が過飽和なら parcel のみを enthalpy 保存で飽和調節してから持ち上げる。
これは参照場の診断であり、実際の大気から水を差し引いたり降水台帳に加えたりしない。

```text
dT_parcel/dln(p) = (R_d*T_parcel + L_v*q_s)
                  / (cp + L_v^2*q_s/(R_v*T_parcel^2))
```

LCL と浮力ゼロ点は圧力間で補間し、`ln(p)` 上で少なくとも二次の parcel 積分を使う。
浮力は `g*(T_parcel-T_env)/T_env` で計算し、CAPE と CIN を別々に積算する。
最初の LFC から最初の上側中立点を対流域とする。対流調節の支持域は最下層からその中立点まで。
中立点が model top の外なら領域を top で切り、top 到達を診断する。
初版の起動条件は丸め誤差を除く `CAPE>0`。CIN は診断のみで、経験的な CIN 起動閾値は追加しない。

暫定参照場は `Tref0=T_parcel`、`qref0=RH_ref*q_s(Tref0,p)`。
対流域外は恒等写像。基準は `RH_ref=0.8`、`tau_conv=7200 s` とする。
温湿度の生の緩和量をそのまま使わず、以下の柱保存条件を満たすよう補正する。

### 6.2 深い対流と浅い対流

対流域の有効質量 `mu[k]=M[k]*f[k]` を使う。内部層の `f=1`、境界層は補間による `0<=f<=1`。
降水に対応する水量と加熱量を次の符号で定義する。

```text
Dq = sum mu*(q-qref0)
DT = (cp/L_v)*sum mu*(Tref0-T)
```

- **深い対流:** `Dq>0` かつ `DT>0`。初版は共通緩和時定数と一様温度 offset による
  enthalpy 補正を採用する。

  ```text
  delta_Tref = L_v*(Dq-DT)/(cp*sum mu)
  Tref = Tref0 + delta_Tref
  qref = qref0
  ```

  `qref` は補正後の温度で再計算しない。指定 RH は暫定 parcel 参照場への条件であり、
  補正後の RH を厳密に固定する条件ではない。過飽和が残れば後段の飽和調節が処理する。
  upstream Isca の場合分けによる湿度緩和時定数の変更は初版で採用せず、この差を fixture に記録する。
- **浅い対流:** `Dq<=0` かつ `DT>0`。対流頂を下げ、水増減の柱積分がゼロになる位置を探す
  shallower 型を採用する。最上部の参加率を補間し、`sum mu*(qref-q)=0` とする。
  同じ有効質量で温度参照場の offset を求め、`sum mu*(Tref-T)=0` を満たす。
  成立する深さがなければ湿潤調節はしない。負の降水をゼロにするだけの処理では代用しない。
- **その他:** 湿潤対流なし。CAPE 不足、熱条件不成立、浅い解なしを別の理由として記録する。

有限 step の緩和は、固定参照場に対して同じ係数を温度・水蒸気に使う。

```text
a = dt/(tau_conv+dt)
Delta T[k] = a*f[k]*(Tref[k]-T[k])
Delta q[k] = a*f[k]*(qref[k]-q[k])
rain_conv_amount = -sum M*Delta q >= 0
sum M*(cp*Delta T + L_v*Delta q) = 0
```

`a` は backward Euler 型の有限 increment とし、`dt/tau>1` の外挿を避ける。
温度 floor、不適切な参照場、solver 不収束は検査し、必要なら substep を細分化する。
有限時間緩和なので一回で CAPE がゼロになること、冪等性、全層が目標 RH になることは要求しない。

### 6.3 乾燥対流との関係

`convection.kind=none|dry_adjustment|simple_betts_miller` とする。
SBM 選択時は既存の乾燥対流調節を前処理として一度実行し、乾燥超断熱だけを除く。
その後の状態から parcel を再診断して SBM を実行する。これにより水ゼロの極限でも
乾燥不安定層を放置しない。乾燥前処理は温度だけを変え、過程 budget を SBM と分離する。
driver で旧乾燥対流をさらにもう一度呼ばない。

LSC-only 比較は `convection.kind=none`、乾燥対流併用比較は `dry_adjustment` で明示する。
潜熱加熱後には新たな乾燥不安定が生じ得るため、受理状態の乾燥安定性を一律に要求しない。
乾燥前処理の前後、SBM 前後、最終飽和調節後の安定度を区別して報告する。

## 7. Splitting・受理・再開

湿潤時の標準順序:

```text
state n と台帳を保存
  -> 既存 integrator で力学＋灰色放射の候補を作る（全 tracer を輸送）
  -> 以下の物理 substep を dt の合計まで繰り返す
       静水圧・境界層係数を再診断
       境界層＋顕熱＋蒸発/露＋地表熱＋bucket を陰的に更新
       乾燥対流前処理（設定時）
       SBM（設定時）
       instantaneous saturation adjustment
       対流＋格子スケール降水を bucket/海へ移し、容量超過を流出へ
  -> 静水圧再診断・非負性・飽和・水/過程エネルギー収支を検査
  -> state と accepted budget をまとめて受理
```

対流後に飽和調節を置き、格子スケール凝結が先に不安定の水をすべて除く順序を標準にしない。
順序依存性は column で比較し、最終段には必ず飽和調節を残す。
凝結・降水の state 更新を dynamics RHS に入れず、純粋な輸送と物理 increment の境界を保つ。

全体の結合は一次精度。物理だけ細分化しても力学／放射との splitting 誤差は消えないので、
大気 dt と物理 substep 幅の両方を変えて検証する。
30分の力学 dt は Phase 12 の結果を自動的に引き継がず、潜熱加熱が強いケースで再評価する。

全体 retry では全 tracer、`Ts`、bucket、降水・蒸発・流出・海洋台帳を元へ戻す。
棄却した呼出しの雨を accepted 積算へ足さない。計算負荷 counter は分離する。
物理更新後は RHS cache を無効化し、semi-implicit の残差は力学＋放射候補について保存する。
初期化や restart 読み込みで黙って凝結させず、初期過飽和を診断して最初の受理 step で扱う。

新しい checkpoint layout は成分数・名前・順序・役割、全 Q、`Ts`、bucket、
継続に必要な累積水台帳を保存し、config fingerprint と照合する。
旧 layout は旧単一 tracer config のまま読み書きできるようにする。
旧 dry checkpoint へ勝手に水蒸気を追加する restart は拒否し、湿潤初期場の生成とは分ける。
湿潤の途中停止・再開は同じ step 分割で連続実行と一致することを要求する。

## 8. 保存量・診断・出力

### 8.1 水の台帳

格子面積を `A` とし、大気水量は `W_atm=sum A*M*q_v`、陸水量は
`W_land=sum A*f_land*W`。受理 step ごとに次を照合する。

```text
Delta W_atm = integral_A_time(E_cell - P_conv - P_lsc)
Delta W_land = integral_A_time[f_land*(P_conv + P_lsc - E_land - runoff_land)]
Delta W_ocean = integral_A_time[(1-f_land)*(P_conv + P_lsc - E_ocean)
                              + f_land*runoff_land]
Delta(W_atm + W_land + W_ocean) = 0
```

全陸時は最後の式の海洋台帳を累積系外流出に置き換える。
`W_ocean` は無限の海の絶対水量ではなく、初期からの受け払い差である。
任意の受動 tracer については外部 source がなければ成分ごとに全球質量保存を検証する。

### 8.2 エネルギー

物理 substep の一次検証量を `H_m=sum M*(cp*T+L_v*q_v)` とする。
凝結・SBM は H_m 保存、境界層・地表交換は `H_m+KE+C_s*Ts` の過程収支を検証する。
水の比熱と降水の位置／運動エネルギーを省いた液相のエネルギー基準を明記し、
bucket の水量変化で `C_s` を変えない。

既存の診断エネルギーには `L_v*W_atm` を追加した `E_dry+E_latent+E_surface` を併記する。
[ADR 0018](adr/0018-dry-convection-and-boundary-layer-coupling.md) の静水圧再診断の有限差分帰属を使い、
`Delta H_m` と診断全エネルギー差を区別する。水収支や局所潜熱保存が閉じても、
既存の力学・有限上端・時間積分の全系残差が解消したとは判断しない。

### 8.3 出力

- 成分 metadata と `tracer_diagnostics.csv`: 各 tracer の質量、min/max、輸送・拡散・物理収支。
- `moisture_diagnostics.csv`: 可降水量、RH、蒸発／露、対流／格子スケール降水、潜熱、
  bucket、流出、海洋受け払い、水収支残差、希薄性診断。
- `moist_convection_diagnostics.csv`: CAPE/CIN、LCL/対流頂、深い／浅い／非起動の割合、
  最大温湿度 increment、H_m 残差、model top 到達。
- 選択柱の `moist_column.csv` と surface snapshot: 温湿度・参照場・過程前後差・bucket・雨量。
  offline の雨量／可降水量／bucket 分布と column profile を描く小さな可視化 script も含める。
- solver 反復、物理 substep 数、retry 理由、受理 dt、過程別時間・allocation。

瞬時 flux `[kg/m2/s, W/m2]` と受理区間の水量／熱量 `[kg/m2, J/m2]`、
tile 面積あたりと cell 面積あたり、全球積分を区別する。
`mm/day` は水量密度による表示用換算とし、内部状態は SI 単位で扱う。
気候平均は実際の受理区間で重み付けする。既存 CSV schema は旧単一 tracer 設定で維持する。
live Web の field selector や Frame の多成分化は後続とし、既存 FrameV2 は旧 preset の意味を維持する。
湿潤 preset を現行 gateway の allowlist に追加しない。

## 9. 設定とファイル構成案

設定例（未実装）:

```text
physics.kind = gray_radiation
moisture.kind = dilute_water
tracers.names = passive,water_vapor
tracers.passive.role = passive
tracers.water_vapor.role = water_vapor
tracers.water_vapor.initial_relative_humidity = 0.6
moisture.condensation = saturation_adjustment
moisture.surface_exchange = bulk
moisture.maximum_physics_substep_s = 300
convection.kind = simple_betts_miller
convection.relaxation_time_s = 7200
convection.reference_relative_humidity = 0.8
boundary_layer.kind = bulk_k_profile
surface.air_exchange_coefficient_w_m2_k = 0
surface.hydrology.kind = bucket
surface.hydrology.capacity_kg_m2 = 150
surface.hydrology.initial_fraction = 0.5
```

既定は moisture/hydrology ともに `none`、tracer registry 未指定時は既存受動 tracer 一つ。
新オプションを省略した config の canonical serialization・fingerprint・軌道を維持する。
湿潤の標準実行は gray＋surface＋bulk 境界層に限定する。
輸送のみ／蒸発なしの凝結試験は明示的に選択できるようにし、独立 column kernel は放射なしでも呼べる。
Held–Suarez との混在、SBM に水蒸気がない設定、重複名／役割、不正な時定数・RH・容量、
非対応の水蒸気拡散、有限陸に bucket のない標準湿潤実験は設定検証で拒否する。

新規 kernel 候補:

- `physics/moist_thermodynamics.hpp/.cpp`
- `physics/saturation_adjustment.hpp/.cpp`
- `physics/simple_betts_miller.hpp/.cpp`
- `physics/surface_hydrology.hpp/.cpp`
- `physics/moist_physics_coupling.hpp/.cpp`
- `dynamics/tracer_registry.hpp` と成分アクセス helper

header は `include/myplanetsim/`、実装は `src/` に置く。
既存 `boundary_layer` を拡張し、湿潤専用の重複した鉛直拡散実装は作らない。
変更対象には config、driver/workspace、reconstruction/flux/coupling/diffusion、diagnostics、
checkpoint、CLI dispatch、tests/CMake を含む。column 数×成分数に応じて workspace を再利用する。

## 10. 実装順序と受け入れ条件

| 段階 | 作業 | 受け入れ条件 |
|---|---|---|
| P13.01 | ADR 0019、希薄近似・保存量・台帳・分割順、旧 baseline と upstream fixture の固定 | q の定義、液相エネルギー基準、乾燥前処理、SBM の差異、checkpoint 契約を実装前に確定 |
| P13.02 | registry、多成分 state/derived/workspace、輸送・両 integrator | 1/2/4 成分、定数保持、成分順序交換、ゼロ tracer、受動成分の独立性、旧結果の回帰 |
| P13.03 | 水に適した保存的拡散・非負輸送、metadata/checkpoint | 不均一 ps・地形・乾湿境界の保存と非負性、restart、未知 schema／成分不一致の拒否 |
| P13.04 | 共通熱力学・独立飽和調節 | 解析微分、乾燥／飽和／過飽和、非均一 M/Pi、enthalpy 保存、飽和残差、冪等性 |
| P13.05 | 蒸発・露・bucket と境界層の陰的結合 | 蒸発冷却、露による昇温、空／満 bucket、供給制約、陸海 fraction 極限、同一 flux による水・熱保存 |
| P13.06 | SBM の parcel・深い／浅い対流 kernel | LCL/断熱極限/CAPE、深浅分岐、非降水の水保存、H_m 保存、参照 fixture、dt/tau が大きい場合 |
| P13.07 | 両 driver へ湿潤 splitting と subcycling を結合 | retry 時の雨・W・台帳 rollback、cache 無効化、同一刻みで restart 一致、乾燥極限、旧経路回帰 |
| P13.08 | 診断 CSV、独立湿潤 column、offline plot、CI preset | schema・単位・面積重み、放射対流 column、水・局所熱収支、時間／鉛直収束 |
| P13.09 | 全球比較行列・係数感度・性能 | 下記の 1日／30日 matrix を完走し、許容刻みと適用範囲を登録 |
| P13.10 | 1200日 aquaplanet pilot と検証報告 | 水収支・非負性・再開・長期統計を報告し、未達の気候／エネルギー条件を明示 |

P13.02--03 の多成分輸送を先に完成させ、P13.04--06 は grid/I/O 非依存の column 単位で検証してから
全球結合へ進む。大きな潜熱加熱を輸送の不具合の調査と同時に導入しない。

### 10.1 数値検証の開始基準

- kernel の水・enthalpy 保存残差は `<=1e-12*scale + absolute_tolerance` を目標とする。
  scale は初期 reservoir と交換量の大きい方とし、絶対 tolerance は水 `1e-12 kg/m2`、熱 `1e-6 J/m2` を開始値にする。
  非線形境界層はまず `1e-10*scale` を許容し、反復 tolerance への収束を確認する。
- 飽和調節後は `q-q_s <= 1e-12 + 1e-10*q_s`、水蒸気は保存補正後に非負。
  調節のみの試験では dt に依存しない。非凝結層は変更しない。
- 30日の全球水収支残差は `max(初期大気＋陸水量, 積算絶対交換量)` に対して `1e-9` 以下を開始目標とする。
  大きい潜熱 flux の相殺だけで見かけの小残差を作らず、局所交換と全球台帳を両方確認する。
- 独立 column は dt を3段階以上細分化し、切替点のない滑らかな実験で結合の観測次数 `>=0.8`。
  飽和／bucket 空／対流起動の切替実験はイベント時刻・積算量の収束を別に測る。
- K=20/40/80、地表細分化 `S=5` で parcel/CAPE、P/E、境界層深さを比較する。
  初期湿度は `RH*q_s(T,p)` から各格子で再生成し、単純な level index コピーをしない。
- 保存・非負性・restart は必須 gate。外部モデルとの降水分布の完全一致は要求しない。
  fixture との違いは希薄式・仮温度・緩和補正の差と数値誤差に分けて説明する。

### 10.2 全球実験行列

開始 preset は N=6、地表細分化 K=20、灰色放射を Phase 12 と共通にする。
同一格子・放射・境界層係数で以下を比較する。

1. 水ゼロ・蒸発なしの乾燥回帰、および受動 tracer を複数にした輸送のみ。
2. aquaplanet: 蒸発＋格子スケール凝結のみ。
3. aquaplanet: 上記＋乾燥対流調節。
4. aquaplanet: 上記の対流を SBM に置換した標準湿潤ケース。
5. 陸海 fraction＋bucket: 湿潤／乾燥初期 bucket、地形あり／なし、全陸の閉じた水循環。
6. 同期回転: 昼側蒸発・夜側凝結と bucket 貯水を診断する。氷なしの適用範囲を明記する。

1日で異常・収支・所有者を確認し、30日で統計と timestep/係数感度を比較する。
大気 dt は `450/900/1800 s`、物理 substep は `75/150/300 s`。
短期の explicit 基準とも比較し、3反復／5反復など既存 semi-implicit の反復誤差を分離する。
開始時の精度目標案は、1800秒対450秒で質量・面積重み温度 RMS `<=0.5 K`、
全球可降水量と平均 P/E の差が `<=5%`、TOA 差 `<=1 W/m2`。
降水がほぼゼロなら相対誤差を使わず絶対値と onset 時刻を報告する。
30日の乱流場は瞬時場一致で判定せず、時間平均と区間変動を比較する。

物理感度は `RH_ref=0.7/0.8/0.9`、`tau=1/2/4 h`、`Wmax=50/150/300 kg/m2` を一つずつ変更する。
全組合せの直積を必須にせず、基準からの差を調べる。
物理 subcycling・成分数1/2/4の実行時間、最大メモリ、warm-step allocation を記録する。

gate 通過後に1200日 aquaplanet pilot を行い、spin-up を除いた時間帯の P/E、水蒸気、温度、風、
TOA、対流／格子降水比、台帳 drift、極端な局所加熱、retry・substep 数の変化を報告する。
陸海・同期回転はまず30日までとし、1200日での検証を行っていない構成を完了扱いにしない。
平均 P と E の近さだけで平衡とは判断せず、大気／陸 reservoir の傾向と block 平均の変化も確認する。

完了時は `docs/validation/phase-13.md` に実測値・実行コマンド・fingerprint・未達条件を残す。
Phase 13 の完了は、定義した希薄・液相のみの水循環の実装と数値検証を意味する。
実在地球の降水気候、雲・氷の再現、既存全系エネルギー残差の解消とは分けて報告する。

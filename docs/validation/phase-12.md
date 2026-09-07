# Phase 12 検証記録

状態: P12.01--P12.09 の実装・独立 column・全球比較・1200日 pilot を実施。
Phase 11 から引き継いだ全系の力学/時間積分残差の厳密分離は未完了で、climate conformance は宣言しない。
本書の数値は下記の再現コマンドで得た測定値であり、目標値ではない。

## P12.01--P12.03 — 規約と乾燥対流調節

[ADR 0018](../adr/0018-dry-convection-and-boundary-layer-coupling.md) が保存量、共有 flux、
分割順、散逸熱の帰属、旧 gray adapter との所有者切替を実装前に固定した。
`convection.kind` / `boundary_layer.kind` は `physics.kind` と直交し、既定は `none`。
未指定の旧 config は canonical serialization・fingerprint・trajectory を保つ。

対流調節は stack を使う PAVA 型の隣接ブロック併合で、重みは `w=cp*M*Pi`。
`unit.dry_convective_adjustment` が二層解析解、不均一質量・Exner、離れた不安定領域、
上位ブロックの再検査、冪等性、`E_h` 保存、静水圧再診断による `E_dry` 帰属を検証する。

一度も併合しなかった層は平均を計算し直さずに元の温位をそのまま複製する。
`w*theta/w` は一般に元の値へ厳密には戻らないため、この複製をしないと安定柱が
丸め誤差規模で「調節された」と報告され、恒等写像と冪等性の要件を満たさない。
発見の経緯は P12.08 の preset 実行で、`adjusted_block_count=0` のまま
`adjusted_column_count=168`、最大温度増分 `5.6e-14 K` が記録されたこと。
回帰は `unit.dry_convective_adjustment` の不規則な質量・Exner のケースが押さえる。

`unit.dry_mixing_core` は SSP-RK3 と semi-implicit の両方で、試行ステップ末尾での適用、
retry 時の完全な rollback、RHS cache の無効化、restart 分割の一致を検証する。

## P12.04--P12.06 — 境界層と driver 結合

`unit.boundary_layer` は固定 K の解析解（空間2次、backward Euler の時間1次）、
最大値原理、熱・運動量・tracer の収支、散逸熱の恒等式、bulk 交換の中立・安定・不安定・
無風極限、粗度検証、層頂での `K=0`、陸海 fraction の面積平均を検証する。

境界層有効時は gray adapter の `H=lambda*(Ts-Tbottom)` と `sigma>0.7` の Rayleigh drag を
除外する。`surface.air_exchange_coefficient_w_m2_k` の非ゼロ指定は設定検証で拒否する。
`unit.dry_mixing_core` は放射 RHS の顕熱 power と drag work がともに厳密にゼロで、
物理 tendency の運動量成分が全てゼロであることを確認する。

`phase12.dry_mixing_cli` は SSP-RK3 と semi-implicit の両方で
`convection_diagnostics.csv` / `boundary_layer_diagnostics.csv` / `mixing_column.csv` の
schema と有限性を確認し、途中停止と再開後の最終 checkpoint が byte 一致することを確認する。

## P12.07 — 地表付近を細分化した鉛直格子と column 実験

`surface_refined_sigma_coefficients(top, K, S)` は sigma 層厚が地表へ向かって等比で薄くなる
ramp を作る。`S` は最上層と最下層の sigma 厚の比で、`S=1` は既存 uniform ramp と厳密に一致する。
これは ADR 0007 の uniform 族ではないため、CLI/UI の K 再生成は従来どおり拒否される。

`S=5`、model top 1000 Pa、`ps=1e5 Pa`、温位 300 K 一定（乾燥断熱）の柱での診断値:

| K | 最下層中心高度 [m] | 最下層 `dp` [Pa] | 1 km 以下の層数 | 最上層中心高度 [m] |
|---:|---:|---:|---:|---:|
| 20 | 87.19 | 1970.2 | 5 | 18912 |
| 40 | 43.66 | 990.6 | 9 | 19974 |
| 80 | 21.85 | 496.6 | 18 | 20821 |

計画の開始目標（最下層中心 30--100 m、1 km 以下に5層以上）を K=20/40 が満たす。
K=80 はそれより細かい。高度は温度・重力・地形に依存するため、実行ごとに診断して報告する。

`phase12.dry_mixing_column` は driver を通さず、放射 → 陰的境界層 → 対流調節の同じ順序を
1本の柱で再現して次を測る。

- 加熱（`Ts=320 K`、4時間）: 顕熱は上向き、地表 Ri は負、調節後の不安定割合はゼロ、
  地表に接する中立層の厚さが単調に増える。
- 夜間冷却（`Ts=285 K`）: 顕熱は下向き、地表 Ri は正、安定度係数は1未満、
  診断 h は加熱側より小さく、調節が働く層は皆無。
- 放射対流柱（日周期日射、1日、`dt=120 s`）: エンタルピー＋運動エネルギー＋地表蓄熱の変化が
  TOA 正味入射の時間積分と `1e-9 * 積算 |flux|` 以内で一致する。
- 時間精度: `dt=3600/8, /16, /32 s` の `3600/256 s` 基準に対する誤差比の対数が
  いずれも `>=0.8`。結合は一次であり、Strang 二次とは呼ばない。
- 鉛直感度: K=20/40/80 で診断 h が収束し、K=40 の K=80 に対する差は `h` の10%未満。

弱風で `U_g=1 m/s` の場合、診断 h は最下層中心（K=20 で 87.19 m）より上、
最初の内部界面（174.21 m）より下に入り得る。全球 `gray_dry_mixing` の初期状態がこれに当たり、
面積平均 h は 93.1 m で `max K_m = max K_h = 0` だった。K profile が内部界面を1つも覆わず、
混合は地表 flux と対流調節だけが担う。計画が予告した「未解像の浅い境界層」である。
同じ設定を1日回すと面積平均 h は 276.2 m、`max K_h` は `0.53 m2 s-1` まで成長し、K profile が効き始める。
h が最下層より浅い場合の `shallow_unresolved_area_fraction` とは別の未解像状態なので、
両者を分けて読む必要がある。

## P12.08 — 全球 preset と1日行列

`configs/phase12_{gray_convection,gray_boundary_layer,gray_dry_mixing,gray_land_ocean,gray_tidally_locked}.cfg`
を追加した。`gray_convection` は Phase 11 と同じ N=4/K=4 uniform 格子と旧地表交換を保ち、
対流調節だけを加える。残りは N=6、`S=5` の K=20 で、旧交換と旧 drag を新境界層へ置き換える。
陸海・地形は Phase 8 と同じ非本番 Earth fixture（`data/earth/manifest.txt`）を使う。

CI の `phase12.dry_mixing_matrix` は N=2/4、K=4/8 の接続試験として、
4通りの物理組合せの所有者分離、乾燥質量保存、調節前後の不安定度、
陸海 fraction と地形の到達、同期回転の昼夜差、係数感度、時間刻み細分化を確認する。

```sh
cmake --preset release
cmake --build --preset release
python3 tools/compare_dry_mixing.py
```

24ケースすべて完走（AppleClang、Release、2026-09-07）。実行コマンド・上書き設定・数値・
失敗ログは `output/phase12-comparison/{results,comparisons}.json` に保存する。
係数を変えたケースは派生 `.cfg` を同じディレクトリに書き出すので入力を復元できる。
RMS は層質量と cell 面積で重み付けし、格子が同一のケース間だけで比較する。

TOA と h の平均は受理ステップの積算値を実際の積算秒数で割った値であり、
日次 sample の平均ではない。

### 1日の物理行列（N=6、K=20、`dt=450 s`、5反復）

| ケース | 平均温度 [K] | TOA 上向き [W m-2] | 平均 h [m] | 調節前 不安定割合（最大） |
|---|---:|---:|---:|---:|
| 放射のみ（旧交換・旧 drag） | 277.704 | 94.180 | — | 0 |
| 放射＋対流（旧交換・旧 drag） | 277.813 | 94.192 | — | 0.05263 |
| 放射＋境界層 | 277.633 | 94.164 | 264.13 | 0 |
| 全部 | 277.636 | 94.164 | 226.20 | 0.05263 |

調節後の不安定割合は全ケースで厳密にゼロ。調節前の値も併記しており、
調節後ゼロという構造上当然の結果だけを改善の根拠にしていない。
旧交換のままの2ケースは境界層診断がゼロ、新境界層の2ケースは旧顕熱・旧 drag の
積算がゼロで、二重計上がないことを積算量で確認した。

「全部」に対する質量・面積重み温度 RMS は、放射のみ `0.5884 K`、放射＋対流 `0.6265 K`、
放射＋境界層 `0.0611 K`。1日では境界層の効果が対流調節より大きい。

### 数値・係数感度（すべて「全部」構成、`matrix-all` との比較）

| 変更 | 温度 RMS [K] | TOA 差 [W m-2] | 平均 h 相対差 |
|---|---:|---:|---:|
| `dt=1800`・3反復 対 `dt=450`・5反復 | 0.001886 | 0.2600 | -1.16% |
| `dt=1800`・5反復 対 `dt=450`・5反復 | 0.001886 | 0.2600 | -1.16% |
| `dt=900`・5反復 対 `dt=450`・5反復 | 0.000636 | 0.0486 | -0.38% |
| `dt=450`・5反復 対 `dt=30` SSP-RK3 | 0.000603 | 0.0181 | -0.36% |
| uniform K=20（`S=1`） | 1.59322 | -0.1976 | -4.74% |
| `Ri_crit=0.25` | 2.06e-06 | 1.42e-07 | -6.76% |
| `Ri_crit=0.5` | 2.01e-06 | 1.31e-07 | -4.30% |
| `Pr_t=1.5` | 1.35e-06 | 1.08e-08 | -0.002% |
| `U_g=0` | 0.40833 | -0.0102 | -63.83% |
| `U_g=0.5` | 0.16391 | -0.0047 | -30.27% |
| 陸粗度 1/10（`z0m=0.01`, `z0h=0.001`） | 0.09139 | -0.0027 | -8.56% |
| 上端 500 Pa | 0.01396 | -0.0076 | +0.13% |
| 上端 250 Pa | 0.02093 | -0.0112 | +0.20% |

計画の仮登録基準（1800秒対450秒で温度 RMS `<=0.5 K`、TOA 差 `<=1 W m-2`、
平均 h の差 `<=10%`）を満たす。よって1200日 pilot に1800秒を登録する。
Phase 11 の30分能力を自動的に引き継いだのではなく、この構成で測り直した結果である。

`Ri_crit` と `Pr_t` は温度をほとんど動かさない一方、診断 h は数%動く。
この構成では h が最初の内部界面に届かない時間帯が長く、K profile が内部混合へ
寄与しないためで、係数が物理的に無関係であることの証明ではない。
逆に鉛直格子（uniform 対 surface-refined、RMS `1.59 K`）と gustiness の影響は係数感度より大きい。
uniform K=20 では調節前の不安定割合が最後までゼロで、対流調節が一度も働かない。
最下層が厚いほど地表加熱が希釈されるためであり、対流の要否は鉛直格子に依存する。

### 格子・地理

| ケース | 平均温度 [K] | TOA 上向き [W m-2] | 平均 h [m] | 秒/モデル日 |
|---|---:|---:|---:|---:|
| N=6, K=20 | 277.636 | 94.164 | 226.20 | 2.73 |
| N=6, K=40 | 277.649 | 93.960 | 277.02 | 4.81 |
| N=12, K=20 | 277.636 | 94.247 | 226.20 | 9.62 |
| N=12, K=40 | 277.649 | 94.043 | 277.02 | 19.81 |
| `gray_land_ocean` | 278.115 | 95.068 | 617.63 | 2.38 |
| `gray_tidally_locked` | 277.636 | 92.966 | 226.48 | 2.32 |

`gray_land_ocean` の平均 h が大きいのは、地形と陸粗度で地表風・不安定度が強まるため。
異なる格子間の差は解像度依存であり、RMS の収束次数とは解釈しない。

全24ケースで乾燥質量の相対 drift の絶対値は `<=1.99e-16`、tracer 質量変化は厳密にゼロ、
局所熱収支残差の規格化値は `<=5.73e-17`、運動エネルギー恒等式の残差は `<=4.03e-14`、
retry は0回、非有限値・非正の層質量はなかった。

## P12.09 — 30日比較と1200日 pilot

### 30日比較

```sh
build/release/tools/benchmark_dry_mixing configs/phase12_gray_pilot.cfg \
  30 6 20 1800 3 1000 5 all output/phase12-30day/dt1800.csv
build/release/tools/benchmark_dry_mixing configs/phase12_gray_pilot.cfg \
  30 6 20 450 5 1000 5 all output/phase12-30day/dt450.csv
```

30日後の質量・面積重み温度 RMS 差は `1.25253 K`。
1日の精度ゲートを長時間の瞬時場一致へ拡張しない。
両者とも accepted step は要求どおり（1440 / 5760）、retry 0、
乾燥質量相対 drift `1.99e-16`、tracer 質量変化ゼロ、調節後の不安定割合は常にゼロ。
調節前の最大不安定割合は `0.14838`（1800秒）と `0.10845`（450秒）、
平均 h の最大は `1827.81 m` と `1953.18 m`。

### 1200日 pilot

`configs/phase12_gray_pilot.cfg` は N=6、`S=5` の K=20、1800秒、3反復、
線形相対許容誤差 `1e-10`、境界層と対流調節を両方有効にし旧地表交換を無効にした構成。
通常の CLI でも実行できる。過程別の記録には次を使用する。

```sh
mkdir -p output/phase12-pilot
build/release/tools/benchmark_dry_mixing configs/phase12_gray_pilot.cfg \
  1200 6 20 1800 3 1000 5 all output/phase12-pilot/pilot.csv \
  > output/phase12-pilot/pilot.log 2> output/phase12-pilot/pilot-progress.log
python3 tools/summarize_dry_mixing_pilot.py output/phase12-pilot/pilot.csv
```

1200日を完走した。57,600 accepted step、median/最小/最大 accepted dt はすべて1800秒、
retry（CFL・invariant・solver）は全区間で0回。線形反復は accepted step あたり16.3--17.8。
peak RSS は 9.9 MiB、0.4713 秒/モデル日。境界層・対流の column solve はそれぞれ 12,441,600 回。

診断時の最大乾燥質量相対 drift は `6.5618e-15`、tracer 質量変化は全ステップで厳密にゼロ、
非有限値・非正の層質量による停止はなし。最終大気温度は 208.264--355.157 K。

| 指標 | 最初の100日 | 最後の100日 |
|---|---:|---:|
| accepted step 数 | 4800 | 4800 |
| median dt [s] | 1800 | 1800 |
| 線形反復/accepted step | 17.8327 | 17.2654 |
| solver / CFL / invariant retry | 0 / 0 / 0 | 0 / 0 / 0 |
| 時間平均大気温度 [K] | 273.972 | 303.089 |
| 時間平均地表温度 [K] | 309.074 | 344.877 |
| 最大風速の日次値の時間平均 [m/s] | 28.391 | 54.743 |
| 時間平均の面積平均 h [m] | 1903.21 | 2390.88 |
| 積算 TOA 上向き [W m-2] | -63.613 | 6.472 |
| 積算 地表蓄熱 [W m-2] | 116.924 | 41.808 |
| 積算 境界層顕熱 [W m-2] | 34.308 | 41.266 |
| 積算 返還散逸熱 [W m-2] | 0.7255 | 1.4528 |

200日以降（48,000 accepted step）の積算値は TOA 上向き `2.6715 W m-2`、
地表蓄熱 `44.138 W m-2`、境界層顕熱 `41.656 W m-2`、返還散逸熱 `1.1719 W m-2`。
時間平均の面積平均 h は `2161.59 m`（範囲 1728.17--2521.65 m）、
`max K_h` は `276.96 m2 s-1`、`model_top_area_fraction` は全区間で0、
`shallow_unresolved_area_fraction` の最大は `0.010669`。

局所の過程収支は全 accepted step で成立する。境界層の熱収支残差を交換量の絶対和で
規格化した最大値は `4.2887e-17`、運動エネルギー恒等式の残差の最大絶対値は `1.0267e-13 J`、
tracer 質量変化は0、運動量収支残差の最大値は `1157.65 N s`。
この運動量残差は cell ごとの残差の大きさを面積重みで足し上げた未規格化の量で、
符号の相殺がない。全球の地表応力 impulse と比べた比率は pilot の step CSV に
impulse 列がないため本 run からは出せない。1日 `dt=60 s` の
`boundary_layer_diagnostics.csv` では残差 `0.1123 N s` に対し impulse `2.68e6 N s`
（比 `4.2e-8`）だった。規格化した恒等式の検証は `unit.boundary_layer` と
`unit.dry_mixing_core` が担当する。
顕熱の積算 `2.16878e24 J` と地表側の `-2.16878e24 J` は符号反転で16桁一致し、
同一の H を大気と地表へ逆符号で使う要件を満たす。

対流調節後の不安定割合は全 accepted step で厳密にゼロ。
調節前の値は200日以降で最大 `0.14239`（面積・層間等重み）。
日次 sample の瞬時場（調節後）で `d(theta)/dz<0` となる界面割合は200日以降の時間平均で
`0.01330` だが、その最小 `d(theta)/dz` は時間平均 `-5.01e-16 K/m`、最小値 `-5.79e-16 K/m` で
丸め誤差規模である。Phase 11 の同区間は割合 `0.09749`、最小 `d(theta)/dz` の最小値
`-0.010273 K/m` だったので、物理的な超断熱勾配の大きさは13桁以上小さくなった。
割合そのものは丸め符号を数えるため、勾配の大きさと併せて読む必要がある。

`*.steps.csv` は全受理ステップと棄却負荷、`*.climate.csv` は約1日ごとの瞬時診断を記録する。
窓内の集計対象は終了時刻が窓内にある受理ステップで、flux 平均はその積算 J を
実際の積算秒数と球面積で割る。日次の TOA や地表 flux は日周期の alias を含み得るため、
熱収支には全ステップ積算値を使う。

`*.climate.csv` の `unattributed_energy_j` は、離散全エネルギー変化から熱源・drag・diffusion と
対流調節の厳密な `E_dry` 帰属を引いた残りで、1200日後は `2.41697e24 J`。
この列は境界層の寄与を引いていない。境界層は enthalpy 基準（`cp` 基準）の熱量しか
報告せず、`total_energy_j` の `cv*T+Phi+|u|^2/2` 基準への変換を実装していないためである。
同区間の境界層顕熱が `2.16878e24 J` であることから、この値のほぼ全量が
「未帰属」であって新たな欠損ではない。Phase 11 の残差 `2.52765e23 J` とは
定義が異なるので直接比較しない。境界層の `E_dry` 帰属と全系 residual の
力学/時間積分への分離は、いずれも本 Phase で閉じていない。

### 未解決範囲

Phase 11 から引き継いだ全系の離散エネルギー残差は本 Phase でも解消していない。
境界層と対流調節はそれぞれの局所収支を厳密に閉じるが、それは力学と時間積分に由来する
残差の原因を特定・修正したことを意味しない。温度 clipping、Newtonian cooling、
global energy fixer は追加していない。

湿潤対流、凝結、雲、予報 TKE、EDMF、完全な Monin--Obukhov 反復解法、
層頂の countergradient flux と entrainment mass flux は本 Phase の範囲外で、未実装のままである。
`U_g` は未解像の風速に対する明示的パラメータで、静止した加熱柱の交換を可能にする仮定を含む。
上の感度表のとおり結果への影響が大きいため、既定値 `1 m/s` は較正済みの物理定数ではない。

対流調節は温位の再配分だけを行い、対流による運動量・物質輸送を含まない。
結合は一次であり、より高次の物理結合は次の計画で扱う。
上端は有限で、上端感度は表のとおり残る。

## 全回帰

```sh
cmake --build --preset dev
ctest --preset dev --output-on-failure -j 4
cmake --build build/dev --target format-check
cmake -S . -B build/gcc-release -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 -DMPS_WARNINGS_AS_ERRORS=OFF
cmake --build build/gcc-release -j 4
ctest --test-dir build/gcc-release --output-on-failure -j 4
cmake --build --preset asan-ubsan
ctest --preset asan-ubsan --output-on-failure -j 4
```

- AppleClang Debug（`-Werror`）: 96/96。
- GCC 15.2 Release: 96/96。プロジェクト既定の `MPS_WARNINGS_AS_ERRORS=OFF` で実行。
- Clang ASan/UBSan: 96/96。format-check も通過。

`phase12_gate` は `phase12.dry_mixing_cli`、`phase12.dry_mixing_column`、
`phase12.dry_mixing_matrix` の3つ。`unit.dry_convective_adjustment`、
`unit.boundary_layer`、`unit.dry_mixing_core` は unit gate に含まれる。

生の snapshot と全ステップ記録は `output/phase12-*`、小さな比較・pilot 集計は
[phase-12-results.json](phase-12-results.json) に保持する。
GCC executable のパスは環境に合わせて変更する。

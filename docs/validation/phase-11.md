# Phase 11 検証記録

状態: P11.01--P11.09 の実装・全球比較・1200日 pilot を実施。全系の力学/時間積分残差の厳密分離は未完了。climate conformance は宣言しない。

## P11.07 — 診断と再開

`radiation_diagnostics.csv` は現在の状態の面積平均 flux（W m-2）と、前回出力後に
受理した全ステップの積算熱量（`interval_*_j`、J）を分ける。初期行の積算値はゼロ。
`segment_start_time_s` / `segment_start_step` は再開した区間の始点を示す。
再開前後の interval 列を加算すれば積算値を結合できる。初期行を重複して統計に数えない。
retry の理由別 count と放射呼出し数・秒数も出力区間内の合計であり、棄却 attempt の計算負荷を含む。
放射呼出し数・秒数は時間発展 RHS の adapter 評価を対象とし、診断用の再評価と最終 column 出力を含まない。
再開直後には RHS cache を再構築するため、連続実行と呼出し数は一致しない場合がある。

`radiation_stable_dt_s` は受理したステップの評価点における最小制限刻みで、初期行のみ
初期状態の制限刻み。熱源の温度微分がゼロの設定では `inf` が無制約を表す。
`radiation_column.csv` は最終時刻の cell 0 を出力する。`interface` 行は half-level の
圧力と flux、`full` 行は full-level 圧力・温度と層の光学厚さ・加熱を保持し、非該当欄は空欄。
地表温度は `surface_state.csv` と既存 `dry_hydrostatic_surface_v1` checkpoint に保存する。
気候統計と physics 診断は run segment 単位であり、checkpoint に accumulator は格納しない。

検証コマンド:

```sh
cmake --build --preset dev
ctest --preset dev -R 'phase11|radiati' --output-on-failure
ctest --preset dev --output-on-failure -j 4
cmake --build build/dev --target format-check
```

CLI テストは SSP-RK3 と centered semi-implicit の両方について、開始時刻 12345 s、
終了時刻 12385 s、10 s 刻みで最初のステップから再開する。最終 checkpoint と column CSV は
byte 一致、各受理区間の flux・積算熱量・retry counter は完全一致を確認する。
強い日射と小さい地表熱容量で CFL retry を起こす独立テストは、受理時間と入射短波の
積算を比較し、棄却した時間を積算へ含めないことを確認する。

停止時には診断間隔に達していなくても最後の区間を出力する。SSP-RK3 と semi-implicit の双方で、
100ステップごとの出力設定を1ステップで停止し、最後の積算値を一度だけ出力する回帰を追加した。

## P11.08 — 全球 preset と1日比較

`configs/phase11_{transparent_surface,gray_uniform,gray_shortwave,gray_land_ocean,gray_tidally_locked}.cfg`
を追加。陸海・地形は Phase 8 と同じ非本番の Earth fixture を用いる。CI は N=2/4、K=4/8。
透明な新旧 physics の120秒 trajectory、短波吸収の再配分、同期回転方向、乾燥質量と界面保存を確認した。
`unit.rhs_allocation` は灰色放射込みの warm RHS の heap allocation がゼロであることを検証する。

```sh
cmake --preset release
cmake --build --preset release
python3 tools/compare_radiation.py
```

20ケースすべて完走（AppleClang 21、Release、2026-09-07、`2b4ba9d` の計測ツール）。実行コマンド・数値・失敗ログは
`output/phase11-comparison/results.json`、温度 RMS 比較は `comparisons.json` に保存する。
RMS は各層を等重み、水平を面積重みとする。空間解像度の異なるケースは全球平均と範囲の感度を示し、
異なる格子間の RMS 収束次数とは解釈しない。

N=6、K=10、1800秒・3反復と450秒・5反復の1日後の差は温度 RMS `1.3840e-4 K`、
TOA `1.2876e-4 W m-2`。1800秒・5反復対30秒 SSP-RK3 は `1.5024e-4 K`、
`1.3683e-4 W m-2`。1800秒で3反復対5反復の差は `1.8275e-5 K`、4対5反復は
`2.1567e-6 K` で、時間刻みの差より小さい。全1800秒ケースの median accepted step は1800秒。
この節の比較は1日基準に限定し、30日・1200日の瞬時場一致を意味しない。長時間結果は後述。

| ケース | 平均温度 K | TOA 上向き W m-2 | 質量相対 drift | 秒/モデル日 |
|---|---:|---:|---:|---:|
| time-1800-i3 | 277.965097 | 78.257206 | 0 | 0.2073 |
| time-1800-i5 | 277.965097 | 78.257207 | 0 | 0.2583 |
| time-450-i5 | 277.965145 | 78.257335 | 0 | 1.0101 |
| explicit | 277.965148 | 78.257343 | -1.41e-13 | 6.4232 |
| space-n6-k20 | 277.852080 | 77.798385 | 0 | 1.8552 |
| space-n6-k40 | 277.738354 | 77.660768 | 0 | 3.7230 |
| space-n12-k10 | 277.965133 | 79.221618 | 0 | 3.9688 |
| space-n12-k20 | 277.852067 | 78.762660 | 0 | 7.5766 |
| gray_shortwave | 278.934909 | 31.983904 | 0 | 0.9913 |
| gray_land_ocean | 278.302146 | 80.144033 | 0 | 0.9941 |
| gray_tidally_locked | 277.965581 | 78.258120 | 0 | 0.9831 |
| top-500 | 277.975456 | 78.246807 | 1.98e-16 | 0.9803 |
| top-250 | 277.980605 | 78.241977 | 0 | 0.9840 |

同じ N=6/K=10 放射ケースの1800秒・3反復は30秒 explicit 比 30.99 倍、放射 RHS adapter の割合は 22.9%、peak RSS は 6.5 MiB。単発の短時間計測であり、長時間の速度保証ではない。
上端1000/500/250 Paの感度は表のとおり。有限上端に対して完全不変とは主張しない。

## P11.09 — column 補強と長時間 pilot

column の補強試験:

- `B(p)=200+0.001*(p-1000)` の独立解析積分を使い、K=10/20/40/80 の
  TOA/surface flux と層加熱の誤差減少次数がともに `>=1.9`。
- 一層温室平衡の ±5 K 摂動から解析平衡へ戻ることを `1e-4 K` で確認。
- 表面冷却の centered 積分と滑らかな時変日射を検証。centered は次数 `>=1.7`、
  SSP-RK3 は `>=2.7`。夜・terminator の輸送テストと、滑らかな区間の次数試験を区別する。
- 固定 seed 11 の100本の妥当なランダム柱で SW/LW 別の telescoping を検証。

### pilot の登録と失敗例

`configs/phase11_gray_pilot.cfg` は N=6、K=10、1800秒、3反復、線形相対許容誤差 `1e-10`。
通常の CLI でも実行できる。詳細な solver/対流不安定度の記録には次を使用する。

```sh
mkdir -p output/phase11-pilot
build/release/tools/benchmark_radiation configs/phase11_gray_pilot.cfg \
  1200 6 10 1800 3 1000 output/phase11-pilot/tight.csv \
  > output/phase11-pilot/tight.log 2> output/phase11-pilot/tight-progress.log
python3 tools/summarize_radiation_pilot.py output/phase11-pilot/tight.csv
```

最初は線形相対許容誤差 `1e-8` で試し、896.244日に停止した。
`dt=0.01373291015625 s` で初期非線形残差 `2.8822230617559398e-9` に対して最終残差
`3.516111323706257e-9` となり、増大判定を通過できず、次の半減が設定最小刻み0.01秒を下回った。
失敗ログは `output/phase11-pilot/pilot-progress.log`、受理ステップは `pilot.csv.steps.csv`。
保存ゲート・非線形増大判定・最小刻みは変更せず、線形許容誤差だけを厳しくした。
これで同じ1200日実験が完走したことは確認できるが、停止の原因が線形誤差だけだったとの証明ではない。

再実行は57,610 accepted steps、median 1800秒、最小112.5秒、solver retry 26回、
CFL/invariant retry 0回。診断時の最大乾燥質量相対 drift `6.5618e-15`、
界面収支残差の最大絶対値 `6.5392e-15 W m-2`。非有限値・非正の層質量による停止はなく、
最終大気温度は207.477--352.908 K。

| 指標 | 最初の100日 | 最後の100日 |
|---|---:|---:|
| accepted step 数 | 4800 | 4801 |
| median dt [s] | 1800 | 1800 |
| 最小 dt [s] | 1800 | 450（終端合わせ） |
| 線形反復/accepted step | 18.0373 | 16.6609 |
| solver retry | 0 | 0 |
| 時間平均大気温度 [K] | 262.392 | 287.873 |
| 時間平均地表温度 [K] | 308.410 | 342.207 |
| 最大風速の日次値の時間平均 [m/s] | 27.944 | 56.291 |

200日以降の面積・層間等重みの対流不安定割合は平均0.09749。各 sample の最小
`d(theta)/dz` の平均は `-0.007137 K/m`、最小値は `-0.010273 K/m`。
放射対流平衡とは呼ばず、乾燥対流調節・境界層を後続の物理検証として残す。

`*.steps.csv` は全受理ステップと棄却負荷を記録する。窓内の集計対象は終了時刻が窓内にある
受理ステップで、flux 平均はその積算 J を実際の積算秒数と球面積で割る。
`*.climate.csv` は約1日ごとの瞬時診断。温度等の平均は実際の sample 時刻で台形積分する。
日次の TOA や surface flux は日周期の alias を含み得るため、熱収支には全ステップ積算値を使う。

### 30日比較と未解決範囲

上のコマンドで日数を30、刻み/反復を1800/3と450/5へ置き換えて比較した。
30日後の温度 RMS 差は `1.61338 K`。日次平均温度の平均は263.216/263.238 K、
日次 TOA の平均は -40.641/-40.581 W m-2、最大風速の範囲は0.025--25.926/0.025--34.451 m/s。
1日の精度ゲートを長時間の瞬時場一致へ拡張しない。1日ゲートと1200日 median の範囲で
30分刻み能力を登録する。

最終の離散大気エネルギー変化から熱源・drag・diffusion の積算寄与を引いた残差は
`2.52765e23 J`。これは放射界面保存誤差ではない。全系 residual の力学と時間積分への
厳密な分離は未完了であり、surface の時間積分 residual と非線形 solver residual を別途記録する。
この残差、有限上端・鉛直 quadrature 感度、対流不安定度があるため climate conformance と
統計平衡を宣言しない。温度 clipping、Newtonian cooling、global energy fixer は追加していない。

全ステップ積算からの200日以降の TOA 上向き flux は 1.843440 W m-2、surface 蓄熱は 2.227520 W m-2。最後100日はそれぞれ 4.944279、-0.139995 W m-2。

### 全回帰

- AppleClang 21 Debug（`-Werror`）: 90/90。停止時出力修正の最終変更後も関連5テストを再検証。
- GCC 15.2 Release: 90/90。プロジェクト既定の `MPS_WARNINGS_AS_ERRORS=OFF` で実行。
  `ON` では既存の Frame/大気 state fixture 等の部分初期化警告でビルドが止まる。
  新しい省略可能 config メンバーの空初期値は明示したが、全既存 fixture の警告修正は本変更に含めない。
- Clang ASan/UBSan: 最終変更を含む90/90。format-check も通過。

```sh
cmake -S . -B build/gcc-release -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-15 -DMPS_WARNINGS_AS_ERRORS=OFF
cmake --build build/gcc-release -j 4
ctest --test-dir build/gcc-release --output-on-failure -j 4
cmake --build --preset asan-ubsan
ctest --preset asan-ubsan --output-on-failure -j 4
cmake --build build/dev --target format-check
```

生の snapshot と全ステップ記録は `output/phase11-*`、小さな比較・pilot 集計は
[phase-11-results.json](phase-11-results.json) に保持する。GCC executable のパスは環境に合わせて変更する。

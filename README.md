# MyPlanetSim

MyPlanetSim は、惑星定数・自転・大気組成・加熱条件・地形を変更できる、C++ 製の架空惑星向け全球大気シミュレーターを段階的に開発するプロジェクトです。水平格子には cubed-sphere、鉛直座標には hybrid sigma-pressure 座標を採用します。

最初から完全な GCM を作るのではなく、cubed-sphere 上の球面輸送、乾燥 3D 力学、地形、理想化物理の順に、解析解・標準テスト・保存則を通過した機能だけを積み上げます。初期検証に使った全球 shallow-water は Phase 15 で廃止しました。現在の対象は乾燥・湿潤の静水圧全球大気コア、理想化物理と、その保存済み結果の可視化です。

## ドキュメント

- [開発・検証計画](docs/plan.md): 数値方式の初期方針、開発段階、各段階の検証項目、完了条件、参考文献
- [Phase 0 実装計画](docs/phase-0-plan.md): 科学・ソフトウェア基盤をコミット単位に分けた作業順と受け入れ条件
- [Phase 1 実装計画](docs/phase-1-plan.md): cubed-sphere 幾何と球面輸送をコミット単位に分けた作業順と受け入れ条件
- [Phase 2 実装計画（履歴）](docs/phase-2-plan.md): 廃止済み全球 shallow-water の標準試験と水平離散化比較
- [Phase 3 実装計画（履歴）](docs/phase-3-plan.md): 廃止済み対話実行経路と初期 Web UI の設計
- [Phase 4 実装計画](docs/phase-4-plan.md): 鉛直 1D の C++ column core と、C++ 検証後に行う offline profile 可視化
- [Phase 5 実装計画](docs/phase-5-plan.md): 地形なし乾燥 3D 静水圧コアを C++ で先に完成し、最後に live visualizer を更新する設計
- [Phase 6 実装計画](docs/phase-6-plan.md): 固定地形、傾いた hybrid 面の pressure gradient、山岳 benchmark と最小の地形入力設計
- [Phase 7 実装計画](docs/phase-7-plan.md): Held--Suarez 型の理想化乾燥物理、source budget、長時間気候統計の限定実装
- [Phase 8 実装計画](docs/phase-8-plan.md): 惑星・軌道設定、Earth 地形、fractional land--ocean surface と理想化熱強制
- [Phase 9 実装計画](docs/phase-9-plan.md): 演算子・CFL・平衡ベンチマークの修正、静的 cache と性能ゲート、production gate の unblock
- [Phase 10 実装計画](docs/phase-10-plan.md): 高速な Lamb・重力波の semi-implicit 化、鉛直 mode 分解、30 分時間刻みの精度・性能ゲート
- [Phase 11 実装計画](docs/phase-11-plan.md): 短波・長波の灰色放射、大気・地表の熱収支、column から全球結合への実装・検証順序
- [Phase 12 実装計画](docs/phase-12-plan.md): 既存モデルを参考にした乾燥対流調節・境界層、陰的地表交換、過程別収支と検証順序
- [Phase 13 実装計画](docs/phase-13-plan.md): 複数 tracer、希薄水蒸気・蒸発・即時降水、SBM 型対流、陸面 bucket と水・潜熱収支の設計
- [Phase 14 実装計画](docs/phase-14-plan.md): RHS 評価削減、物理の頻度分離、OpenMP、flux 融合、配列・SIMD・tile 最適化の比較と採用基準
- [Phase 15 実装計画](docs/phase-15-plan.md): shallow-water・対話実行・Chromium テストの廃止、地形と期間平均場を閲覧する viewer への整理
- [科学・数値規約](docs/conventions.md): 単位、座標、添字、符号、誤差・保存量の共通定義
- [Phase 0 検証報告](docs/validation/phase-0.md): toolchain、テスト、時間収束、restart の検証結果
- [Phase 1 検証報告](docs/validation/phase-1.md): cubed-sphere 幾何、演算子、球面輸送の検証結果
- [Phase 2 検証報告（履歴）](docs/validation/phase-2.md): 廃止前の shallow-water 標準試験、保存量、方式比較
- [Phase 3 検証報告（履歴）](docs/validation/phase-3.md): 廃止前の native gateway、Web UI、control protocol
- [Phase 4 検証報告](docs/validation/phase-4.md): hybrid 鉛直 column、診断、restart の検証結果
- [Phase 5 検証報告](docs/validation/phase-5.md): 乾燥 3D 静水圧コアの検証結果と production gap
- [Phase 6 検証報告](docs/validation/phase-6.md): 固定地形と山岳 benchmark の検証結果と production gap
- [Phase 7 検証報告](docs/validation/phase-7.md): Held--Suarez forcing、online 統計、CI gate と保留中の production matrix
- [Phase 8 検証報告](docs/validation/phase-8.md): 惑星・軌道、fractional surface、結合熱収支の bounded gate と production gap
- [Phase 9 検証報告](docs/validation/phase-9.md): 正しさの修正、平衡 benchmark、性能 gate と既知の数値的制約
- [Phase 10 検証報告](docs/validation/phase-10.md): semi-implicit 実装の固定 baseline と段階的検証結果
- [Phase 11 検証報告](docs/validation/phase-11.md): 灰色放射の解析・再開検証、全球比較、1200日 pilot と適用範囲
- [Phase 12 検証報告](docs/validation/phase-12.md): 乾燥対流調節・境界層の column 検証、地表細分化格子、4通り比較と1200日 pilot
- [Phase 13 検証報告](docs/validation/phase-13.md): 希薄水蒸気、SBM 型対流、陸面 bucket と長期 pilot
- [Phase 14 検証報告](docs/validation/phase-14.md): RHS・物理・並列化候補の性能比較と採否
- [Phase 15 検証報告](docs/validation/phase-15.md): 期間平均 dataset、restart、viewer、削除内容とブラウザー不要 gate
- [cubed-sphere panel ADR](docs/adr/0001-cubed-sphere-panel-conventions.md): panel、edge、向き、flux 符号規約
- [shallow-water state ADR（廃止済み）](docs/adr/0002-shallow-water-state-and-staggering.md): Phase 2 の予報変数、staggering、flux/source 分割の履歴
- [水平離散化 ADR（廃止済み）](docs/adr/0003-shallow-water-horizontal-discretization.md): Phase 2 の Rusanov 法と compatible 候補の比較履歴
- [鉛直座標 ADR](docs/adr/0005-hybrid-vertical-coordinate-and-column-state.md): hybrid `A/B`、column state、鉛直 mass flux の規約
- [乾燥 3D 結合 ADR](docs/adr/0006-dry-hydrostatic-state-and-coupling.md): 3D state、水平・鉛直 flux、pressure gradient、実装順序の規約
- [固定地形 ADR](docs/adr/0008-fixed-orography-and-lower-boundary.md): 不変な surface geopotential、下部境界、地形 source と I/O 境界
- [惑星 forcing・fractional surface ADR](docs/adr/0010-planetary-forcing-and-fractional-surface.md): benchmark forcing、軌道、陸海 fraction、ゼロ深度 surface 熱収支
- [ベクトル演算子・時間刻み ADR](docs/adr/0011-vector-operators-and-time-step-normalization.md): 球面ベクトル Laplacian、cell 単位 CFL、地表 reservoir の安定性制約
- [静的 cache・性能ゲート ADR](docs/adr/0012-static-grid-cache-and-performance-gates.md): 時間不変な格子量の前計算、RHS workspace、正しさと分離した性能ゲート
- [平衡初期場 ADR](docs/adr/0013-balanced-benchmark-initial-states.md): `steady` の定義、UMJS14/JW06 の惑星整合、DCMIP 座標の修復
- [低散逸・reference balance ADR](docs/adr/0014-low-dissipation-and-reference-balanced-dry-core.md): Lamb CFL と移流散逸の分離、DCMIP terrain rest の離散保存
- [重力波 semi-implicit ADR](docs/adr/0016-semi-implicit-gravity-wave-integration.md): fast/slow 分割、反復法、solver failure、30分超の探索契約
- [期間平均 visual dataset ADR](docs/adr/0028-period-mean-visual-dataset.md): 受理 step の時間平均、保存形式、restart と読み取り専用 viewer

## ビルドと実行

CMake 3.25 以上、Ninja、C++20 対応コンパイラを使用します。通常の開発では
CMake preset を利用します。

```sh
cmake --preset dev
cmake --build --preset dev
./build/dev/my_planet_sim --config configs/phase0_ode.cfg
./build/dev/my_planet_sim --config configs/phase1_solid_body.cfg
ctest --preset dev
cmake --build build/dev --target format-check
```

`release` と `asan-ubsan` preset も同じ configure/build/test 手順で使用できます。

実行ファイルは `experiment.kind` に応じて、Phase 0 の製造 ODE、Phase 1 の
全球 cubed-sphere tracer 輸送、Phase 4 の独立した hybrid 鉛直 column、Phase 5 以降の
乾燥・湿潤 3D 静水圧コアと理想化物理を実行します。輸送ケースは診断を標準出力へ、
cell snapshot を `output.directory/tracer.csv` へ出力します。乾燥・湿潤コアでは既存の
診断 CSV、checkpoint と進捗・停止処理を維持しています。

```sh
./build/dev/my_planet_sim --config configs/phase4_isothermal.cfg
./build/dev/my_planet_sim --config configs/phase4_manufactured_transport.cfg
./build/dev/my_planet_sim --config configs/phase5_isothermal_rest.cfg
./build/dev/my_planet_sim --config configs/phase7_held_suarez_short.cfg
./build/dev/my_planet_sim --config configs/phase8_earth_geography.cfg
```

Held--Suarez run は `physics_diagnostics.csv` と `climate_statistics.csv` を
`output.directory` に出力します。`phase7_held_suarez_short.cfg` は CI 用であり、1200 日の
climate conformance を主張する preset ではありません。
standalone dry-hydrostatic run は10秒ごとに、step、モデル日、進捗率、経過時間、ETAを
標準エラーへ表示します。間隔は `--progress-interval-s SEC` で変更でき、`0` で無効化できます。
`Ctrl-C` はpartial diagnosticsと指定済みcheckpointを書いて `result.status = cancelled` で終了します。

Phase 8 の surface energy balance run は `surface_state.csv` と
`surface_diagnostics.csv` を追加で出力します。checked-in Earth surface product は
結合経路を検証する解析的 CI fixture であり、NOAA ETOPO/GSHHG から生成した production
Earth data ではありません。この境界と production matrix の保留項目は
[Phase 8 検証報告](docs/validation/phase-8.md)に記録しています。

## 開発状況

Phase 15 までの限定した実装・検証範囲は完了しています。球面輸送、乾燥・湿潤 3D
静水圧コア、固定地形、理想化物理、惑星・軌道・陸海 surface と期間平均 viewer を実装済みです。
初期の shallow-water と対話実行 gateway は履歴文書だけを残して廃止しました。現在の状態と
再現手順は [Phase 15 検証報告](docs/validation/phase-15.md)を参照してください。

Phase 10 では、乾燥静水圧コアの高速な Lamb・重力波だけを semi-implicit に扱い、陽に残る移流 CFL を
守りながら長い時間刻みを可能にしました。`N=12`, `K=20` の Held--Suarez は通常 1800 秒で進み、
explicit 200 秒に対して 5.78 倍の速度を示しました。反復停止則は残差許容値ではなく固定反復数で、
その回数は反復誤差が時間打ち切り誤差を下回る最小値として実測で決めています。1200 日 pilot は
20分57秒で完走し、day 412の一時的な強風では陽的移流CFLに従って安全に短縮した後、1800秒へ
復帰しました。長期的なsolver劣化と非有限値はなく、Phase 10は完了しています。設計と実装順は
[Phase 10 実装計画](docs/phase-10-plan.md)、実測値は
[Phase 10 検証報告](docs/validation/phase-10.md) に定めています。

本プロジェクトは研究・開発段階の実験的ソフトウェアです。Phase 5--7 の production 規模の
長時間の production climate matrix は未実施です。DCMIP 2-0-0 の terrain-following rest は
reference-state pressure force で離散保存され、移流散逸は Lamb CFL から分離されましたが、
現時点の結果を production climate validation として扱わないでください。

### 地形・期間平均 viewer

viewer はシミュレーターを起動せず、通常 CLI が保存した `VisualDatasetV1` をローカルで読みます。
小格子の例は次のコマンドで生成できます。

```sh
just dataset
# just を使わない場合
./build/dev/my_planet_sim --config configs/phase15_viewer_rest_n4.cfg --progress-interval-s 0
```

出力先は `output/phase15-viewer-rest-n4/viewer/` です。`manifest.json`、`terrain.bin` と
`means/period_*.bin` を含みます。集計は受理 step の終端値を実時間で重み付けし、
`statistics.start_time_s` と `statistics.period_s` が定める半開区間ごとに保存します。
途中終了した期間は coverage と `complete=false` を持つため、完全な期間平均と区別できます。

viewer の依存関係と開発 server は次のように準備します。

```sh
just setup
just dev
# または: cd web && npm ci && npm run dev --workspace @myplanetsim/ui
```

表示された URL を開き、dataset の `viewer` フォルダーを選択してください。地形、集計期間、
変数、model level と cell を切り替えられ、2D map・3D globe・鉛直 profile が選択を共有します。
計算開始、初期条件編集、瞬時 frame の再生は行いません。全自動 gate は Node と CTest で
ファイル契約・数値処理を検証し、Canvas/WebGL の実描画は手動 smoke の対象です。

```sh
just check
```

## 将来のファイル構成

```text
MyPlanetSim/
├── CMakeLists.txt          # C++ ビルドとテスト登録
├── README.md               # プロジェクト入口
├── docs/                   # 設計、数式、検証記録、ADR
├── include/myplanetsim/    # 公開 C++ API
├── src/                    # 実行入口、格子、数値演算、力学、物理、I/O
├── apps/                   # 実行プログラム
├── tests/                  # unit / convergence / regression テスト
├── configs/                # 再現可能な実験設定
├── web/                    # 読み取り専用 viewer と保存 dataset protocol
└── tools/                  # 可視化・比較などの補助ツール
```

構成は各開発段階の開始時に必要最小限だけ追加し、計画上の候補を先回りして空ディレクトリ化しません。

## ライセンス

このプロジェクトは [MIT License](LICENSE) の下で公開します。

灰色放射は `physics.kind = gray_radiation` で選択できます。短波吸収、地表反射、上下長波と
大気・地表の熱収支を同じ界面 flux から計算し、診断 CSV と地表温度を含む checkpoint を出力します。
`configs/phase11_gray_uniform.cfg` は1日の基準、`configs/phase11_gray_pilot.cfg` は検証済みの
N=6/K=10・1200日候補です。pilot は1800秒の median accepted step で完走しましたが、
全系エネルギー残差と対流不安定層が残るため、気候平衡の再現とは位置付けていません。
詳細な比較と再現コマンドは [Phase 11 検証報告](docs/validation/phase-11.md) にあります。

乾燥対流調節は `convection.kind = dry_adjustment`、bulk 境界層は
`boundary_layer.kind = bulk_k_profile` で選択できます。両者は `physics.kind` と直交し、
既定は `none` なので既存 config の結果と fingerprint は変わりません。境界層を有効にすると
灰色放射 adapter の旧顕熱交換と下層 Rayleigh drag は無効になり、地表熱・運動量・tracer の
鉛直交換を境界層が所有します。`configs/phase12_gray_dry_mixing.cfg` は地表付近を細分化した
K=20 の1日基準、`configs/phase12_gray_pilot.cfg` は1200日候補です。過程別収支は
`convection_diagnostics.csv` / `boundary_layer_diagnostics.csv` / `mixing_column.csv` に出力します。
測定値・感度・適用範囲は [Phase 12 検証報告](docs/validation/phase-12.md) にあります。

Phase 13 の bounded 湿潤結合は `configs/phase13_moist_ci.cfg` で実行できます。
`tracer_diagnostics.csv`、`moisture_diagnostics.csv`、
`moist_convection_diagnostics.csv`、`moist_column.csv`、`moist_surface_state.csv` を出力します。
柱・水循環・bucket の dependency-free SVG は次のコマンドで生成できます。

```sh
python3 tools/plot_moist_column.py output/phase13_moist_ci/moist_column.csv
```

P13.09 の全球行列・30日感度と P13.10 の1200日 aquaplanet pilot は完走済みです。
pilot は foreground segment runner で実行・再開できます。

```sh
python3 tools/run_phase13_pilot.py
# graceful stop 後の再開
python3 tools/run_phase13_pilot.py --resume
```

適用範囲、実測値、出力構成は [Phase 13 検証報告](docs/validation/phase-13.md) にあります。

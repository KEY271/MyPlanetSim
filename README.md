# MyPlanetSim

MyPlanetSim は、惑星定数・自転・大気組成・加熱条件・地形を変更できる、C++ 製の架空惑星向け全球大気シミュレーターを段階的に開発するプロジェクトです。水平格子には cubed-sphere、鉛直座標には hybrid sigma-pressure 座標を採用します。

最初から完全な GCM を作るのではなく、cubed-sphere 上の球面輸送、全球 shallow-water、乾燥 3D 力学、地形、理想化物理の順に、解析解・標準テスト・保存則を通過した機能だけを積み上げます。初期の到達目標は、乾燥・静水圧・全球大気力学コアと理想化強制です。湿潤過程、放射、化学、海洋結合はその後の拡張対象です。

## ドキュメント

- [開発・検証計画](docs/plan.md): 数値方式の初期方針、開発段階、各段階の検証項目、完了条件、参考文献
- [Phase 0 実装計画](docs/phase-0-plan.md): 科学・ソフトウェア基盤をコミット単位に分けた作業順と受け入れ条件
- [Phase 1 実装計画](docs/phase-1-plan.md): cubed-sphere 幾何と球面輸送をコミット単位に分けた作業順と受け入れ条件
- [Phase 2 実装計画](docs/phase-2-plan.md): 全球 shallow-water、標準試験、水平離散化比較の作業順と受け入れ条件
- [Phase 3 実装計画](docs/phase-3-plan.md): C++ 制御付き独立 Web UI、初期条件編集、3D/2D cubed-sphere 表示の設計とコミット列
- [科学・数値規約](docs/conventions.md): 単位、座標、添字、符号、誤差・保存量の共通定義
- [Phase 0 検証報告](docs/validation/phase-0.md): toolchain、テスト、時間収束、restart の検証結果
- [Phase 1 検証報告](docs/validation/phase-1.md): cubed-sphere 幾何、演算子、球面輸送の検証結果
- [Phase 2 検証報告](docs/validation/phase-2.md): shallow-water 標準試験、保存量、方式比較の検証結果
- [cubed-sphere panel ADR](docs/adr/0001-cubed-sphere-panel-conventions.md): panel、edge、向き、flux 符号規約
- [shallow-water state ADR](docs/adr/0002-shallow-water-state-and-staggering.md): 予報変数、staggering、flux/source 分割
- [水平離散化 ADR](docs/adr/0003-shallow-water-horizontal-discretization.md): 基準 Rusanov 法と compatible 候補の比較と採否

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
全球 cubed-sphere tracer 輸送、Phase 2 の全球 shallow-water を実行します。
輸送ケースは診断を標準出力へ、cell snapshot を `output.directory/tracer.csv`
へ出力します。shallow-water ケースは診断と保存量 drift、計算コストを標準
出力へ、cell snapshot を `output.directory/shallow_water.csv`、区間診断を
`output.directory/diagnostics.csv` へ出力します。

```sh
./build/dev/my_planet_sim --config configs/phase2_williamson2.cfg
```

Phase 0 の科学・ソフトウェア基盤、Phase 1 の cubed-sphere 幾何・球面受動
輸送、Phase 2 の全球 shallow-water はローカル検証を完了しています。次の
開発対象は Phase 3 の C++ 制御付き独立 Web 可視化 UI です。

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
├── web/                    # 独立 UI、local C++ control gateway、共有 protocol
└── tools/                  # 可視化・比較などの補助ツール
```

構成は各開発段階の開始時に必要最小限だけ追加し、計画上の候補を先回りして空ディレクトリ化しません。

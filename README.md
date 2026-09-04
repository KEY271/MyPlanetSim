# MyPlanetSim

MyPlanetSim は、惑星定数・自転・大気組成・加熱条件・地形を変更できる、C++ 製の架空惑星向け全球大気シミュレーターを段階的に開発するプロジェクトです。水平格子には cubed-sphere、鉛直座標には hybrid sigma-pressure 座標を採用します。

最初から完全な GCM を作るのではなく、cubed-sphere 上の球面輸送、全球 shallow-water、乾燥 3D 力学、地形、理想化物理の順に、解析解・標準テスト・保存則を通過した機能だけを積み上げます。初期の到達目標は、乾燥・静水圧・全球大気力学コアと理想化強制です。湿潤過程、放射、化学、海洋結合はその後の拡張対象です。

## ドキュメント

- [開発・検証計画](docs/plan.md): 数値方式の初期方針、開発段階、各段階の検証項目、完了条件、参考文献
- [Phase 0 実装計画](docs/phase-0-plan.md): 科学・ソフトウェア基盤をコミット単位に分けた作業順と受け入れ条件
- [Phase 1 実装計画](docs/phase-1-plan.md): cubed-sphere 幾何と球面輸送をコミット単位に分けた作業順と受け入れ条件
- [Phase 2 実装計画](docs/phase-2-plan.md): 全球 shallow-water、標準試験、水平離散化比較の作業順と受け入れ条件
- [Phase 3 実装計画](docs/phase-3-plan.md): C++ 制御付き独立 Web UI、初期条件編集、3D/2D cubed-sphere 表示の設計とコミット列
- [Phase 4 実装計画](docs/phase-4-plan.md): 鉛直 1D の C++ column core と、C++ 検証後に行う offline profile 可視化
- [Phase 5 実装計画](docs/phase-5-plan.md): 地形なし乾燥 3D 静水圧コアを C++ で先に完成し、最後に live visualizer を更新する設計
- [科学・数値規約](docs/conventions.md): 単位、座標、添字、符号、誤差・保存量の共通定義
- [Phase 0 検証報告](docs/validation/phase-0.md): toolchain、テスト、時間収束、restart の検証結果
- [Phase 1 検証報告](docs/validation/phase-1.md): cubed-sphere 幾何、演算子、球面輸送の検証結果
- [Phase 2 検証報告](docs/validation/phase-2.md): shallow-water 標準試験、保存量、方式比較の検証結果
- [cubed-sphere panel ADR](docs/adr/0001-cubed-sphere-panel-conventions.md): panel、edge、向き、flux 符号規約
- [shallow-water state ADR](docs/adr/0002-shallow-water-state-and-staggering.md): 予報変数、staggering、flux/source 分割
- [水平離散化 ADR](docs/adr/0003-shallow-water-horizontal-discretization.md): 基準 Rusanov 法と compatible 候補の比較と採否
- [鉛直座標 ADR](docs/adr/0005-hybrid-vertical-coordinate-and-column-state.md): hybrid `A/B`、column state、鉛直 mass flux の規約
- [乾燥 3D 結合 ADR](docs/adr/0006-dry-hydrostatic-state-and-coupling.md): 3D state、水平・鉛直 flux、pressure gradient、実装順序の規約

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
全球 cubed-sphere tracer 輸送、Phase 2 の全球 shallow-water、Phase 4 の独立した
hybrid 鉛直 column、Phase 5 の地形なし乾燥 3D 静水圧コアを実行します。
輸送ケースは診断を標準出力へ、cell snapshot を `output.directory/tracer.csv`
へ出力します。shallow-water ケースは診断と保存量 drift、計算コストを標準
出力へ、cell snapshot を `output.directory/shallow_water.csv`、区間診断を
`output.directory/diagnostics.csv` へ出力します。

```sh
./build/dev/my_planet_sim --config configs/phase2_williamson2.cfg
./build/dev/my_planet_sim --config configs/phase4_isothermal.cfg
./build/dev/my_planet_sim --config configs/phase4_manufactured_transport.cfg
./build/dev/my_planet_sim --config configs/phase5_isothermal_rest.cfg
```

Phase 0 の科学・ソフトウェア基盤、Phase 1 の cubed-sphere 幾何・球面受動
輸送、Phase 2 の全球 shallow-water、および Phase 3 の C++ 制御付き Web 可視化
基盤はローカル検証を完了しています。Phase 3 の検証結果と既知 gap は
[Phase 3 検証報告](docs/validation/phase-3.md)に記録しています。

### Web UI と native gateway

[`just`](https://github.com/casey/just) が利用できる場合、repository root から次の1コマンドで
native gateway と live UI を起動できます。必要な C++ build または `node_modules` がなければ
初回だけ自動セットアップし、接続用URLを標準出力へ表示します。終了は `Ctrl-C` です。

```sh
just dev
```

既定portを変更する場合は `just dev 9876 5174`、offline/mock UIだけなら `just ui`、
gatewayだけなら `just gateway` を使用します。全検証は `just check` で実行できます。

`just` を使わない場合、UI は C++ を必要としない mock/offline preview として起動できます。

```sh
cd web
npm ci
npm run build
npm run dev --workspace @myplanetsim/ui
```

native simulator を local gateway から制御する場合は、gateway 起動時に binary と
allowlist preset を固定します。起動時に表示される JSON の `port` と `token` を
API request に使います。

```sh
MPS_SIMULATOR_BINARY="$PWD/../build/dev/my_planet_sim" \
MPS_REST_PRESET="$PWD/../configs/phase3_rest_n4.cfg" \
MPS_DRY_PRESET="$PWD/../configs/phase5_visualizer_rest_n4.cfg" \
MPS_RUN_ROOT="$PWD/../.runs" \
npm run start --workspace @myplanetsim/gateway
```

`MPS_DRY_PRESET` は任意です。指定すると capabilities が Phase 5 の
dry hydrostatic preset を additive に返し、UI の preset selector に現れます。gateway は
起動時に `--describe-control` と preset descriptor の compatibility を検査し、preset ごとの
N 上限、edit 可否、累積 published byte 予算 (既定 256 MiB) を強制します。shallow-water
preset は Phase 3 の FrameV1 と Gaussian edit を維持し、dry hydrostatic preset は FrameV2 を
publish して edit を拒否します。

gateway API は loopback と Bearer token を要求します。run request を送信した後、
`POST /api/v1/runs/{id}/cancel` で停止し、`GET /api/v1/runs/{id}/bundle` で request、
control text、event log、diagnostics をまとめて取得できます。frame は
`GET /api/v1/runs/{id}/frames/{sequence}` で `frame.ready` event 後に取得します。
実際の small N=4 run、frame 0/後続 frame、cancel の再現は次で検証できます。

```sh
npm run test:live --workspace @myplanetsim/gateway
```

UI を実ブラウザで開いて gateway 接続・preset・model level・column profile まで確認する
gate も用意しています。Chromium が未取得なら skip します。

```sh
npx playwright install chromium
npm run test:browser --workspace @myplanetsim/ui
```

ブラウザ UI を live gateway へ接続するには、別terminalでUIを起動し、gateway起動時に
表示された `port` と `token` をURL fragmentへ指定します。fragmentはclient生成後にURLから
除去され、tokenはgateway以外へ送信されません。**Vite が最後に表示する token なしの
`http://127.0.0.1:5173/` を開くと offline demo にフォールバックし、shallow-water preset
しか出ません**（UI 上部に警告が出ます）。`just dev` は token 付き URL を自動で開きます。

```text
http://localhost:5173/#token=<token>&gateway=http%3A%2F%2F127.0.0.1%3A<port>
```

dry hydrostatic preset を選ぶと、UI は model level selector、surface pressure/pressure/
potential temperature/temperature/tracer/wind speed の field selector、選択 column の
対数 pressure profile、選択 cell の inspector を表示します。level slice は 2D map と 3D globe が
共有し、cell click は edit ではなく column selection になります。

鉛直解像度 `K` は N と同様に run ごとに指定できます（`K (vertical resolution)`）。
[ADR 0007](docs/adr/0007-interactive-vertical-resolution-override.md) により、preset が
uniform sigma の hybrid 座標を使っている場合に限り、preset の model top を保ったまま
`A_pa[k] = p_top*(1-k/K)`, `B[k] = k/K` を再生成します。stretched な座標を持つ preset は
黙って平坦化されず拒否されます。K を変えた run は fingerprint が変わる別 configuration です。
UI の `Model level` slider は表示する層を選ぶ view 設定で、`K` は run 設定です。

tokenを指定しない場合は安全なoffline/mock clientを使用します。offline/mock client は
shallow-water preset のみを提供し、dry hydrostatic は native gateway を必要とします。Chromium/Firefox/WebKitの
自動matrixは次の拡張対象です。

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

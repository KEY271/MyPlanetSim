# MyPlanetSim

MyPlanetSim は、惑星定数・自転・軌道・大気組成・放射・地表条件・地形を変更できる、C++20 製の研究用全球大気シミュレーターです。水平格子に cubed-sphere、鉛直座標に hybrid sigma-pressure 座標を用い、乾燥・希薄水蒸気の静水圧大気と、保存済み期間平均データを表示する Web viewer を実装しています。

現行機能、設定、数値方式、出力形式、制約は [実装仕様書](docs/specification.md) に集約しています。`configs/`、`tests/`、`tools/` に残る `phase*` という名前は既存の再現・回帰 fixture の識別子であり、進行中の開発計画ではありません。

## クイックスタート

必要なものは CMake 3.25 以上、Ninja、C++20 対応コンパイラです。viewer には Node.js と npm も必要です。

```sh
cmake --preset dev
cmake --build --preset dev
./build/dev/my_planet_sim --config configs/phase5_isothermal_rest.cfg
ctest --preset dev
```

期間平均 dataset の生成と viewer の起動:

```sh
just dataset
just setup
just dev
```

全チェック:

```sh
just check
```

実行結果は設定の `output.directory` 以下に生成され、Git 管理されません。

## ライセンス

[MIT License](LICENSE)

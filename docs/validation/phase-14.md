# Phase 14 検証記録

状態: P14.01 を実装・検証中。P14.02 以後は未着手。
速度の数値目標を満たすための diffusion、残差許容値、物理パラメータの調整は行わない。

## P14.01: 計測基盤

- RHS の initial / iteration / final を実呼出しで計数。discarded は破棄試行の部分集合であり総数に再加算しない。retry 前に共有する initial は一回だけ数える。線形 solve で先に失敗する試行は RHS を増やさない。
- 任意の RHS profiling を追加。診断/setup、reconstruction、flux/CFL、column coupling、source、diffusion、physics、result copy の排他的8領域。reconstruction 内の pack/limiter、線形 solve 内の modal/Helmholtz はまだ個別に分離していない。
- observer は full RHS を呼ばないので diagnostic full RHS は0。observer の物理診断は積分 RHS の回数に含めない。
- `benchmark_phase14.py` は同じコマンドの warm-up 後に通常5回、profile5回、allocation別1回を実行し、中央値/範囲、環境、binary/config SHA256 を保存する。
- moist benchmark の allocation 計数は `MPS_COUNT_ALLOCATIONS=1` 指定時だけ有効。通常runの0をallocationゼロの根拠にしてはならない。
- [精度 manifest](phase-14-manifest.json) を数値候補の試験前に登録。

```sh
cmake --build --preset release -j 4
ctest --preset release --output-on-failure -j 4
python3 tools/benchmark_phase14.py --output output/phase14/baseline.json -- \
  build/release/tools/benchmark_moist configs/phase13_moist_pilot.cfg \
  1 6 20 1800 3 20000 5 output/phase14/baseline.csv
```

Release の105テストが成功。計測有効/無効を含む継続積分のstate一致を確認。
Phase 13 の既存1200日CSVを修正版で再集計し、(200日,1200日] は
86,400,000秒、100日ブロック10個（最初200–300、最後1100–1200）となった。
積算fluxの区間が境界をまたぐ場合は推定配分せずエラーにする。

P14.01 の残作業: 発達後checkpointを含む比較行列、詳細physics/copy/solver計測、
高解像度・thread条件別baseline。初日の測定だけでSBM性能や方式の採否は判定しない。

### 発達後の比較 fixture

`benchmark_moist` の末尾に `INITIAL_CHECKPOINT SOURCE_CONFIG` を指定すると、
元configのfingerprintで検証した湿潤checkpointから指定日数だけ積分する。
物理設定・tracer registry・格子を維持し、時間積分と診断の設定だけ変更できる。
CLIのtop/refinementは元の座標と整合することを確認した上で、丸めを含めた元の
鉛直座標をそのまま用いる。これは比較専用のimportであり、production restartの契約は変えない。
蒸発・降水はimport時点の累積値を差し引き、水台帳も開始時の海洋・外部流出を含む。
共通runnerはcheckpointのSHA256も記録する。

```sh
python3 tools/benchmark_phase14.py --output output/phase14/developed-baseline.json -- \
  build/release/tools/benchmark_moist configs/phase13_moist_pilot.cfg \
  1 6 20 1800 3 20000 5 output/phase14/developed.csv \
  output/phase13_moist_pilot/pilot.chk configs/phase13_moist_pilot.cfg
```

1200日checkpointから1800秒進めるsmokeでRHS=5/step、retry=0、対流降水比0.9845。
`phase14.fixture` は実際に生成したcheckpointを使い、経過秒数・区間蒸発・RHS=4の
2反復import・同条件再実行のsnapshot byte一致・物理/鉛直座標不一致の拒否を検証する。
既存Release 105テストと新しいfixtureテストが成功。

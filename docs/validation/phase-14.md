# Phase 14 検証記録

状態: P14.01 を実装・検証中。P14.02 の固定反復比較を実装。
predictor/IMEX 比較、P14.03 以後は未着手。
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

P14.01 の残作業: 地形等を含む比較行列、詳細physics/copy/solver計測、
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

## P14.02: 固定2/3/5反復の最初の比較

[比較契約 ADR 0020](../adr/0020-phase14-fixed-iteration-comparison.md) と
[実測JSON](phase-14-ici-results.json) を登録。数値コード・preset・線形許容値は変更していない。

```sh
python3 tools/compare_phase14_ici.py --output output/phase14/ici-initial
python3 tools/compare_phase14_ici.py --output output/phase14/ici-developed \
  --checkpoint output/phase13_moist_pilot/pilot.chk \
  --source-config configs/phase13_moist_pilot.cfg
```

両ケースでN=6/K=20、450/900/1800秒×2/3/5反復、各1モデル日を測定した。
各点はwarm-up＋通常5回＋profile5回＋allocation別1回。
同一binary、Release、OpenMP OFF。電源・CPU affinityは固定していない。
JSONに中央値と範囲、binary/config/checkpoint SHA256、build設定を記録した。
18条件すべてでretry/physics retry=0、RHSは2/3/5反復に対して4/5/7回/step。
通常/profile/allocationで比較対象の物理診断とRHS回数が完全一致した。
Release buildと107テストが成功。比較器の回帰テストは参照質量による重み付け、
格子・行数不一致、非有限温度の拒否を含む。

以下の温度RMSは450秒・5反復基準。速度は同一dtの3反復に対する実測である。

| 開始状態 | 力学dt [s] | 2反復 [s/day] | 3反復 [s/day] | speedup | 2反復 温度RMS [K] |
|---|---:|---:|---:|---:|---:|
| 初期場 | 450 | 4.5662 | 5.0693 | 1.1102 | 0.000000793 |
| 初期場 | 900 | 3.0011 | 3.2441 | 1.0810 | 0.0002637 |
| 初期場 | 1800 | 2.5323 | 2.6570 | 1.0493 | 0.0008910 |
| 1200日checkpoint | 450 | 4.7425 | 5.2265 | 1.1021 | 0.005092 |
| 1200日checkpoint | 900 | 3.0990 | 3.3697 | 1.0874 | 0.07370 |
| 1200日checkpoint | 1800 | 2.6280 | 2.7343 | 1.0404 | 0.09404 |

現行schedulerは最大300秒を等分するため、450秒の実physics intervalは225秒、
900/1800秒では300秒になる。dtを変えた差にはこの差も含まれる。
同一dtの5反復に対する2反復の温度RMSは、発達後で順に
0.005092 / 0.004543 / 0.04122 K。反復の影響と時間刻みの影響を混同しない。
発達後1800秒・2反復のPW差は0.00923%、降水差0.2265%、TOA差0.1774 W/m2。
測定した1日gateは18条件すべて通過したが、0.09404 Kは登録値0.1 Kに近い。

**採否:** 2反復をこの湿潤1日比較の候補として残す。標準presetは3反復を維持する。
既存線形波テストでは2反復の反復誤差が時間離散誤差を上回るため、全ケースへの一括採用はしない。
1800秒での全体短縮は初期場4.70%、発達後3.89%に留まり、RHSの20%削減を
全体の20%削減として報告しない。目標達成のための追加調整はしていない。

P14.02 の残作業はpredictor/IMEX比較、terrain/baroclinic等の適用範囲、75秒physics基準、
CAPE/降水tail・規格化enthalpy残差・長期統計。P14.01の詳細計測・高解像度baselineと、
P14.03以後の実装も未完了。今回の測定だけでこれらを完了扱いにしない。

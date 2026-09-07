# Phase 11 検証記録

状態: P11.07 まで実装。全球比較と長時間 pilot は未完了。

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

未完了: P11.08 の全球 preset・空間時間感度・性能測定、P11.09 の 1200 日 pilot、
対流不安定度と長時間の solver 分布、GCC/Clang/ASan/UBSan の最終全回帰。
既存力学のエネルギー残差について climate conformance を宣言しない。

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

20ケースすべて完走（AppleClang 21、Release、2026-09-07）。実行コマンド・数値・失敗ログは
`output/phase11-comparison/results.json`、温度 RMS 比較は `comparisons.json` に保存する。
RMS は各層を等重み、水平を面積重みとする。空間解像度の異なるケースは全球平均と範囲の感度を示し、
異なる格子間の RMS 収束次数とは解釈しない。

N=6、K=10、1800秒・3反復と450秒・5反復の1日後の差は温度 RMS `1.3840e-4 K`、
TOA `1.2876e-4 W m-2`。1800秒・5反復対30秒 SSP-RK3 は `1.5024e-4 K`、
`1.3683e-4 W m-2`。1800秒で3反復対5反復の差は `1.8275e-5 K`、4対5反復は
`2.1567e-6 K` で、時間刻みの差より小さい。全1800秒ケースの median accepted step は1800秒。
この登録は1日基準に限定し、30日・1200日の能力はまだ含めない。

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

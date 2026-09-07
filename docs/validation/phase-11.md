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

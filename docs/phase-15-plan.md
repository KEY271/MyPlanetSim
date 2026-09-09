# Phase 15 — shallow-water・対話実行の廃止と地形・期間平均の可視化

状態: 実装中（2026-09-09）。P15.01 完了。
調査基準: `742f11f`。knowledge MCP を `MyPlanetSim` で検索したが該当記録はなかった。
以下のファイル名・設定名・データ形式は、既存と明記したものを除き実装予定である。

## 1. 目的とスコープ

初期検証に使った shallow-water、visualizer からの計算実行、Chromium を使うテストを廃止する。
visualizer は **地形と、事前計算した月平均などの期間平均場の閲覧** に集中する。
各セルの空間分布は保持するが、全時刻の全セル値を保存・転送・再生する機能は作らない。

ユーザーが必要とする操作は、地形または計算結果を開き、集計期間・変数・鉛直層を選び、
地図・地球儀・選択セルの値と鉛直プロファイルを確認することである。
平均は C++ の CLI 計算中に逐次集計し、Web は保存済みの集計値を表示する。
前案の瞬時 FrameV2 の通常 CLI 出力、全フレーム読込、時間スライダー、Play/Pause は採用しない。

Phase 15 の必須範囲:

- shallow-water 本体・専用数値処理・専用設定・テストを削除し、共有処理を分離する。
- Web の計算操作、初期条件編集、gateway、C++ control protocol、mock solver を削除する。
- 地形だけで開けるデータ出力と、セル・層ごとの期間平均出力を追加する。
- 2D 地図、3D globe、セル選択、鉛直プロファイルを読み取り専用 viewer に接続する。
- Chromium/Playwright と live/browser gate を削除し、数値・ファイル契約をブラウザー不要で検証する。

新しい物理、MPI/GPU、瞬時場アニメーション、等圧面補間、年をまたぐ月別平年値、
異常値・分散・極値・降水集計の新規表示、3D 地形の凹凸強調は後続候補とする。
既存の物理診断・収支 CSV、benchmark、CLI の進捗・停止・checkpoint/restart は維持する。
月平均は viewer の中心であり、既存の数値検証を月平均だけに置き換えるものではない。
Phase 14 の長期統計・性能の未完了ゲートは別途残し、この整理で達成済みとは扱わない。

## 2. 現状と設計上の注意

| 現在の実装 | Phase 15 での扱い |
|---|---|
| `web/apps/ui/src/main.tsx` は計算設定、編集、実行、表示をまとめて管理 | データ読込、集計期間、変数、層、セルの選択へ整理 |
| `Map2D`、`Globe3D`、`ColumnProfile`、投影・セル選択・色付けの部品がある | 表示処理を再利用し、編集マーカーと実行 protocol への依存を除去 |
| `FrameV2` は温度・圧力・風などの瞬時値を固定レイアウトで格納 | 配列順と数値処理を参考にし、期間平均の独立した契約へ置換 |
| FrameV2 出力は `apps/my_planet_sim.cpp` の machine control 経路にある | 通常 CLI に期間平均 writer を追加。瞬時 frame 出力は移植しない |
| `ClimateStatisticsAccumulator` は緯度帯・層ごとの面積平均と時間積算 | 既存の気候検証用に維持。セルごとの期間平均 accumulator を別に追加 |
| `SurfaceOrography` / `SurfaceBoundary` に計算用の地形・陸域率がある | 実際に補間・平滑化した計算格子上の値を一度だけ出力 |
| 乾燥 driver が `shallow_water_diffusion.cpp` の `stable_diffusion_time_step()` を使用 | 共通 numerics モジュールへ移してから shallow-water を削除 |
| `FrameV2` に N≤24、K≤30 の制限がある | 新形式は配列サイズ・bytes の上限で検証。対話実行の制限を引き継がない |

## 3. 集計する量と平均の意味

### 3.1 初版のフィールド

| データ | フィールド | 集計 |
|---|---|---|
| 静的地形 | 標高 m、地表ジオポテンシャル m²/s²、陸域率 1 | 平均せず保存 |
| 地表 | 地表気圧 Pa | 各水平セルの時間平均 |
| 大気 | 温度 K、圧力 Pa、温位 K | 各水平セル・model level の時間平均 |
| 風 | 東西風・南北風 m/s、風速 m/s | 成分と瞬時風速をそれぞれ時間平均 |
| 任意の既存場 | 地表温度 K、水蒸気・凝結水など登録済み tracer の混合比 | モデルが持つ場合に限り、その状態から時間平均 |

最低受け入れ対象は地形・温度・圧力・地表気圧・風とする。
任意場は名前・単位・地表/大気の区別をフィールド一覧から取得し、存在しない量をゼロで捏造しない。
tracer は registry の識別子を保存し、配列の先頭成分を無条件に水蒸気と解釈しない。

温度は各状態から診断した温度を平均する。平均した保存変数から温度を再計算しない。
平均風速 `mean(sqrt(u²+v²))` と平均風の大きさ `sqrt(mean(u)²+mean(v)²)` は区別する。
初版の鉛直座標は model level で、圧力も同じ層の平均圧力を表示する。
地形追随層は時刻によって圧力が変わるため、この表示を「等圧面の月平均」と呼ばない。
セルごとの時間平均に面積重みは不要。将来の全球・領域平均ではセル面積を別途重みに使う。

### 3.2 「月」と集計期間

架空惑星へ地球の暦を暗黙に適用しない。初版は SI 秒で指定する固定長の期間平均とする。
設定候補は `statistics.enabled`、`statistics.start_time_s`、`statistics.period_s`、
`statistics.fields`。既存 config は集計無効のまま読めるようにする。
Phase 15 の利用例では `period_s = 2592000`（30×86400秒）を明示し、UI に「30日平均」と表示する。
86400秒を表示上の1日とし、惑星の自転周期・太陽日・公転周期とは区別する。
地球の暦月、惑星年の12分割、季節区間による集計はこの Phase で追加しない。

- spin-up 除外時刻 `start_time_s` を明示する。除外時刻以前の値を平均に混ぜない。
- 窓は `start_time_s + j*period_s` を起点とする半開区間 `[start,end)` とする。
- 各連続区間を別々に保存する。複数年の同じ月をまとめる平年値ではない。
- 計算終了・停止で窓が途中なら、実際の集計範囲・秒数・coverage と `complete=false` を保存する。
  完全な30日平均とは区別して表示し、寄与のない空の窓は出力しない。
- 任意に切り出した期間を UI で再平均する機能は初版に含めない。

### 3.3 逐次集計の数値契約

初版は、物理更新まで完了した **受理済み step の終端値** をその step 区間の代表値とする
時間重み付き矩形則を採用する。`x_i` が `(t_(i-1),t_i]` の代表値、窓が W のとき、

```text
w_i = length([t_(i-1), t_i] ∩ W)
mean(x; W) = sum(w_i * x_i) / sum(w_i)
```

これは瞬時値の単純な標本平均ではなく、離散的な時間平均の近似である。
可変 dt・retry で step 数が変わっても実際の秒数で重みを付ける。
step が窓境界や spin-up 終端をまたぐ場合は重なり時間を各窓に分配する。
集計境界のためにモデルの時間刻みや物理 scheduler を変更しない。
粗い dt の境界近似は上記矩形則に従うことを metadata と検証記録に明記する。

- 初期状態・restart 直後の同時刻通知は積算しない。棄却 step、RK stage、ICI candidate も含めない。
- 既存 observer の `derived` は診断間隔でしか渡らないため、それだけに依存しない。
  受理 step ごとに必要な温度・圧力・風等を再利用または軽量に診断する経路を設ける。
  集計のための full RHS・放射・SBM の再評価や physics cache の変更は禁止する。
- 集計は `diagnostics.interval_steps` や Web 操作に依存させない。
- 現在の窓の積分値と秒数だけを保持し、窓が閉じたら平均を書き出して再利用する。
  メモリは O(セル数×層数×フィールド数)、出力は O(期間数×セル数×層数×フィールド数)。
  step 数に比例した state/frame の保持・出力は行わない。
- NaN/Inf・配列不一致は集計エラーとし、欠測を黙ってゼロ置換しない。

## 4. 保存形式と restart

### 4.1 地形と集計結果の独立したデータ契約

予定構成:

```text
output/viewer/
  manifest.json
  terrain.bin
  means/
    period_000000.bin
    period_000001.bin
```

`manifest.json` は dataset schema version、grid mapping、N/K、panel/セル/層の配列順、
惑星の半径・重力・時間単位、hybrid A/B、設定 fingerprint、地形 source fingerprint、
build 情報、フィールド識別子・単位・shape・格納位置、各期間の予定境界・実集計範囲・
重み秒数・coverage・完全性・ファイルサイズ、平均則を記録する。
binary は little-endian float64 を基本に、magic/version と shape/byte length を検証する。
計算途中の一時ファイルを完成済みデータとして manifest に載せず、原子的に公開する。
地形のみの場合は期間一覧が空でよい。地形にはシミュレーションに使用した値と source を記録する。
陸域率の定義がない地形では当該フィールドを省略する。

既存 `FrameV1/V2` と run/event protocol は新形式へ混在させず削除する。
期間平均を瞬時 frame に偽装しない。旧 binary の配列処理は必要な部分だけを再利用する。
旧 frame 互換 reader・一括変換器は Phase 15 の必須範囲に含めない。
既存の診断 CSV や最終 checkpoint だけから、失われたセルごとの月平均を復元することはできない。
過去の計算成果物は削除せず、新しい期間平均は集計を有効にした計算から取得する。

### 4.2 計算の中断と継続

モデルの checkpoint を維持し、集計有効時は checkpoint に対応する専用の集計 sidecar を保存する。
内容は窓番号・積分値・重み秒数・最終受理時刻/step・集計設定・出力済み期間の識別子とする。
checkpoint 内容の hash と照合し、別状態や別格子の sidecar は読まない。
同じペアからの restart では二重積算・窓の重複出力を防ぎ、連続実行と平均を一致させる。
途中窓を閲覧用に出した場合でも、再開時は sidecar の積分値から継続し、同じ期間 ID の結果を更新する。

sidecar のない旧 checkpoint からもモデルは再開可能とするが、集計は再開時刻から開始し、
窓の欠けた部分は coverage に反映する。過去の寄与を補って完全な平均とは表示しない。
sidecar が存在するのに不整合な場合は明示的なエラーとする。
既存の乾燥・湿潤 checkpoint layout と設定 fingerprint の互換性を回帰検証する。

## 5. 削除・移動するファイルと機能

| 対象 | 作業 |
|---|---|
| `src/dynamics/shallow_water_*`、対応ヘッダー、専用 diagnostics | 全削除。`stable_diffusion_time_step` を先に共通 `numerics/diffusion_stability.*` へ移動 |
| `src/numerics/compatible_operators.cpp`、`src/grid/dual_topology.cpp` と対応ヘッダー | 現在 shallow-water 専用の利用なので関連テストと削除。削除前に残存参照を再確認 |
| `spherical_operators.*` | shallow-water 専用渦位関数を削除。共通の球面・ベクトル演算子は維持 |
| `experiment_config.*` | `kShallowWater`、専用 enum/parameters、読み込み・検証・出力を削除。共通 limiter/reconstruction/diffusion 型は維持 |
| `apps/my_planet_sim.cpp` | shallow-water 実行・CSV・診断分岐と machine control を削除。通常 CLI に集計・地形出力を追加 |
| `configs/phase2_*.cfg`、`phase3_rest_n4.cfg`、`phase6_williamson5.cfg` | 削除。Phase 名だけでなく `experiment.kind` と参照を確認 |
| Williamson 5 の解析地形 | 専用 benchmark と一緒に削除。DCMIP/JW06/linear bell/CSV/Earth の地形は維持 |
| `test_shallow_water_*`、Williamson 2/6、Galewsky、`test_semi_implicit_shallow_water.cpp` | 専用テストと CMake 登録を削除。GMRES/Helmholtz と dry semi-implicit の検証は維持 |
| `test_orographic_responses.cpp`、`test_time_step_normalization.cpp`、config/地形/checkpoint テスト等 | shallow-water ケースだけ削除・fixture 置換。共通・乾燥大気の検証を残す |
| `src/control/`、`include/myplanetsim/control/`、`test_control_request.cpp` | 全削除。`--control-request`、`--event-stream`、`--describe-control` も廃止 |
| `src/io/frame.cpp`、対応ヘッダー・FrameV1/V2 テスト | 新しい dataset I/O に置換後、旧形式を削除 |
| `web/apps/gateway/` | process 起動・停止、capabilities、events、frames、bundle、認証 token、gateway/live テストを全削除 |
| `web/packages/protocol/src/client.ts`、`index.ts` の実行 protocol | SimulationClient、mock、要求・イベント・状態遷移を削除。パッケージを表示データ契約に絞る |
| `web/packages/protocol/src/visual.ts` | shallow-water CSV/派生量・旧 frame decoder を削除し、期間平均の field/slice/sample API に置換 |
| `web/apps/ui/src/editor.ts`、`controller.ts` と専用テスト | 編集・実行 controller を削除し、dataset 読込と閲覧状態の管理へ置換 |
| `main.tsx`、`Map2D.tsx`、`Globe3D.tsx`、CSS | Run/Cancel、preset、N/K 編集、dt/終了時刻/出力間隔、Gaussian 編集・リング、実行診断・bundle・request 表示を削除 |
| 時系列 UI | 全 frame の保持、追尾、再接続、時間スライダー、Play/Pause/Replay を削除。保存済み集計期間の選択へ置換 |
| `interactive_dry_presets.txt` | 削除。`phase5_visualizer_rest_n4.cfg` は専用用途を廃止し、必要な小規模 config は Phase 15 fixture へ整理 |
| `web/apps/ui/test/browser.test.js`、package.json、lockfile | Playwright/Chromium、`test:browser` と関連環境変数・依存を削除 |
| `.github/workflows/ci.yml`、`Justfile`、`vite.config.ts` | live/browser job、`just gateway`、gateway 同時起動、token 注入を削除。`just dev` は viewer 起動へ簡略化 |

ODE・球面輸送・独立鉛直 column の数値検証は今回の削除対象に含めない。
Web のテストが shallow-water のダミーデータに依存する場合は地形・平均場 fixture に置き換える。
`VisualDatasetV1` など表示汎用型と削除する瞬時 `FrameV1` を名前だけで同一視しない。
境界層・湿潤対流で使う `shallow` は別概念なので削除しない。

## 6. viewer の構成

- 起動時はデータを開く画面とする。偽の計算・自動 mock run は行わない。
- フォルダー選択を基本に、manifest と関連ファイルの複数ファイル選択でも読める入口を設ける。
  ローカル File の相対パスを索引化し、絶対パスや親ディレクトリー参照・欠損を拒否する。
  backend・認証 token・シミュレーター binary を不要とする。
- 地形だけで 2D/3D、標高・陸域率、格子境界、セル情報を表示できる。
- 平均場では集計期間・変数・model level を選ぶ。N/K はファイル由来の読み取り専用情報とする。
- セル選択を 2D/3D/鉛直プロファイルで共有し、位置・地形・平均値・単位・平均圧力を表示する。
- 地表量には層選択を適用しない。地形表示は集計期間の有無に依存しない。
- カラーバーと単位、範囲を表示し、期間比較用に色範囲を固定できるようにする。
- 期間の開始・終了、実集計時間、不完全期間の表示を設ける。単に「時刻 t の温度」と表示しない。
- 選択期間のファイルだけを遅延読込し、小さい上限付き cache を使う。
  期間切替や dataset 切替で古い非同期読込が画面を上書きしないようにする。
- 描画・投影・セル選択ロジックと、データの形式・読込・状態管理を分離する。

## 7. 実装順と受け入れ条件

| 順番 | 作業単位 | 主な変更先 | 受け入れ条件 |
|---|---|---|---|
| P15.01 | 依存整理と集計契約の固定（完了） | 共通 numerics、config 設計、dataset schema/ADR | dry core の拡散処理を分離。平均則・暦・配列順・互換方針を固定 |
| P15.02 | セル・層ごとの逐次期間平均（完了） | `diagnostics/period_mean.*`、driver の受理 step hook | 定数場・可変 dt・複数窓・spin-up・retry 除外の数値テストが通る |
| P15.03 | 通常 CLI の集計出力・地形出力・継続（完了） | `io/visual_dataset.*`、CLI、config、集計 sidecar | 地形だけで出力可能。乾燥/湿潤とも期間平均を出力し、restart 前後で一致 |
| P15.04 | TypeScript の保存データ読込（完了） | protocol の dataset/mean API、loader、fixtures | C++ 出力と値・shape・セル順・期間情報が一致し、破損を検出 |
| P15.05 | viewer を地形・期間平均へ移行（完了） | main、Map2D、Globe3D、ColumnProfile、閲覧状態 | 地形だけ/平均場の両方を開け、期間・層・セルの選択が連動 |
| P15.06 | 対話実行と旧 frame を除去 | gateway、control、editor、clients、frame I/O、設定 | viewer と通常 CLI の置換経路が成立し、旧制御・瞬時 frame 参照が消える |
| P15.07 | shallow-water を除去 | dynamics、diagnostics、config、専用設定・テスト、CMake | 対象本体と専用基盤を削除し、残るモデル・共通演算の回帰が通る |
| P15.08 | テスト・CI・起動・文書の最終整理 | CI、package/lock、Justfile、README、旧 ADR、検証報告 | Chromium 不要で全 gate が通り、新しい利用手順を再現できる |

実装はこの順に一項目ずつ進める。置換出力・reader の準備前に現行の表示経路だけを削除しない。
新しい API・CLI option・config は実装時に最終名を揃え、この計画と利用例も更新する。

## 8. 検証計画

### 数値・継続

- 定数場はどの窓・dt でも同じ平均。既知の時間関数で矩形則の期待値と時間細分化時の収束を検査する。
- 異なる長さの step、窓をまたぐ step、複数窓をまたぐ step、spin-up、ちょうど窓終端、
  端数期間、ゼロ寄与、重複通知・逆行時刻を検証する。
- 反対向きの風で平均風と平均風速が区別され、温度・圧力が保存変数からの誤った平均にならない。
- 同一入力の集計 ON/OFF でモデル state と既存収支診断を変えず、診断出力頻度でも平均が変わらない。
- 窓の途中・境界・retry 後の restart を連続実行と比較。sidecar 欠落/不一致と partial 更新も検証する。

### I/O・viewer

- 小格子 C++ CLI → manifest/terrain/mean binary → Node/Vitest decoder の結合テストを設ける。
  少なくとも2期間を短い試験用 period で生成し、温度・圧力・地形・期間情報を直接照合する。
- N/K、flatten order、tracer ID、欠損ファイル、version、bytes 上限、非有限値、
  地形と平均場の fingerprint 不整合を検証する。
- 地形のみ、地表量、大気量、期間・層・セル切替、不完全期間、色範囲固定を単体テストする。
- Chromium/Playwright の代替ブラウザー自動化は導入しない。
  Canvas/WebGL の実描画は手動 smoke で確認し、ブラウザー不要テストが描画を保証するとは扱わない。

### 回帰・性能・記録

- C++ build/CTest/format、既存 GCC/Clang/sanitizer、Web lint/typecheck/test/build を通す。
- 生きているコード・設定・依存・CI に shallow-water/control/gateway/Playwright が残らないことを確認する。
  廃止を説明する履歴文書の語句は許容する。
- 同一期間・dt・field 数で集計 ON/OFF の wall time、診断回数、peak memory、出力 bytes を比較する。
  RHS/physics 再評価を追加していないこと、窓数に比例して accumulator memory が増えないことを確認する。
- `docs/validation/phase-15.md` に削除一覧、schema、再現コマンド、数値結果、性能差、
  手動表示確認と残る制約を記録する。
- 旧 Phase 2/3 の検証報告は履歴として保存。ADR 0002/0003/0004/0007 の廃止・置換対象を明記し、
  現在も有効な格子・力学・物理の設計と区別する。

完了条件は、**CLI で地形と期間平均を作り、シミュレーターなしで viewer に読み込めること**、
および不要な実行・瞬時再生・shallow-water・Chromium 依存を除去した状態で既存大気モデルの検証が通ることである。

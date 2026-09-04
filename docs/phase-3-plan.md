# Phase 3 実装計画 — C++ 制御付き独立 Web 可視化 UI

## 1. 到達点

Phase 3 では、ブラウザから shallow-water の初期条件と安全な実行パラメータを編集し、
native C++ simulator を開始・停止して、計算中の field を 3D 球体と 2D 全球地図で
逐次確認できる local Web application を実装する。

Web UI は C++ source/header へ直接依存せず、C++ simulator も HTTP や UI library に
依存しない。両者の間に localhost 専用 control gateway と version 付き protocol を置く。
UI 自体は mock simulator と保存済み snapshot だけでも単独 build/test/preview でき、
live mode を使うときだけ gateway と native `my_planet_sim` executable を必要とする。

成果物は次を含む。

- React/TypeScript 製 UI と Node.js 製 localhost control gateway
- browser/gateway/C++ 間の version 付き run、event、frame contract
- preset 選択、実行時間、frame interval、初期条件 edit の UI
- click 位置へ Gaussian depth perturbation を加える initial-condition editor
- C++ による実際の shallow-water 時間積分と run/cancel lifecycle
- initial、intermediate、final state の binary frame 出力と逐次表示
- scalar field を cell 単位で描画する interactive 3D globe
- equirectangular projection を初期表示とする interactive 2D global map
- 3D/2D の双方で panel seam と全 cell edge を切り替える cubed-sphere grid overlay
- offline CSV import、診断、provenance、request/result export
- unit、component、protocol、process、browser integration、performance の検証

Phase 3 完了時の「クリックで波を起こす」は描画効果ではない。静止 shallow-water stateへ
局所 depth perturbation を加え、その初期状態を native C++ solver で積分して生じた
gravity wave を表示することを意味する。

## 2. 境界と前提

### 2.1 独立性の定義

同一 repository 内で次の境界を守る。

- `web/` は独自の package、lockfile、TypeScript 設定、test、build output を持つ。
- CMake から npm script を呼ばず、Web build から CMake を呼ばない。
- UI/gateway は C++ header や library へ link/import しない。
- C++ core/application は HTTP、React、Node.js、Web schema validator に依存しない。
- C++ と gateway の結合点は CLI、strict control request、NDJSON event、binary frame のみ。
- C++ CI、Web CI、live integration CI を分け、失敗の責務を判別できるようにする。
- offline UI は Node.js gateway/C++ binary がなくても fixture と local CSV で動作する。
- live UI は gateway 経由でのみ C++ を操作し、browserからprocessやfilesystemへ直接触れない。

従来案の「C++ sourceを変更しない」「実行時backendを持たない」は撤回する。初期条件を
科学計算へ反映するには、C++側に初期状態edit、frame observer、machine-readable CLIを
追加する必要がある。ただしこれらはWeb固有の描画処理ではなく、batch automationからも
利用できる汎用application/core interfaceとして実装する。

### 2.2 Phase 3 の実行範囲

Phase 3 は single-user、local workstation 向けとし、gateway は loopback interface のみに
bindする。remote deployment、共有cluster、multi-user queue、TLS/OAuth、container sandboxは
含めない。それらが必要なら、認証・authorization・resource isolation・quota・auditを伴う
別Phase/ADRで設計する。

### 2.3 含めないもの

- browser内へのshallow-water solverの再実装またはWebAssembly化
- arbitrary shell command、任意binary、任意filesystem pathを受け取るAPI
- source code編集、build、Git操作を行うWeb endpoint
- checkpoint内容を無検証で書き換える機能
- 複数user/hostにまたがるscheduler、remote execution
- 3D volume rendering、鉛直断面、複数layer。3D core完成後に追加する。
- Earth coastlineや地形をsimulation resultにない情報として暗黙に重ねること

## 3. 全体アーキテクチャ

```text
Browser UI (web/apps/ui)
  ├─ offline: CSV / saved frames / mock SimulationClient
  └─ live: HTTP RunRequestV1 + authenticated event stream
                         │
                         ▼
Local control gateway (web/apps/gateway, 127.0.0.1 only)
  ├─ validates protocol and allowlisted preset/overrides
  ├─ creates an isolated run directory
  ├─ writes ControlRequestV1 text
  ├─ spawns one configured my_planet_sim binary without a shell
  ├─ relays NDJSON lifecycle/frame events
  └─ serves only frame files below that run directory
                         │
                         ▼
Native C++ my_planet_sim
  ├─ loads the existing experiment config
  ├─ validates and applies ControlRequestV1
  ├─ constructs the benchmark initial state
  ├─ applies InitialConditionEditV1 in C++
  ├─ integrates with the existing shallow-water driver
  └─ atomically writes FrameV1 + NDJSON events
```

browserはgatewayとだけ通信し、gatewayはC++ executableをchild processとして管理する。
gatewayがsimulation equationsやcubed-sphere fieldを再計算しないこと、UI previewより
C++が返すframe 0をauthoritative initial stateとすることを不変条件にする。

### 3.1 Web workspace

```text
web/
├── package.json                    # npm workspaces と共通script
├── package-lock.json
├── apps/
│   ├── ui/                         # React/Vite browser application
│   └── gateway/                    # local HTTP/process controller
├── packages/
│   └── protocol/                   # RunRequest/Event/Frame metadata schema
├── test/                           # shared fixture/fake simulator
└── e2e/                            # Playwright offline/live scenario
```

`web/apps/ui` は `SimulationClient` interface にだけ依存し、HTTP clientとmock clientを
差し替えられるようにする。`web/apps/gateway` は Node.js built-in HTTP/stream/process APIを
基本とし、process起動ではargument arrayと `shell: false` を使用する。production local modeでは
gatewayがbuild済みUI assetと `/api/v1` を同じoriginから配信する。developmentではViteが
`/api/v1` をloopback gatewayへproxyし、UI単独開発時はmock clientを使う。

### 3.2 C++ 側の責務分離

```text
include/myplanetsim/control/        # control request と initial-condition edit model
include/myplanetsim/io/             # frame/event writer contract
src/control/                        # strict parse、validation、edit適用
src/io/                             # little-endian frame、NDJSON event出力
src/dynamics/                       # observerを受けるsolver driver
apps/my_planet_sim.cpp              # CLI/process boundary、signal/result handling
```

core libraryはnetworkやjob IDを知らない。gatewayのHTTP requestを直接渡さず、applicationが
versioned control requestへ変換してからcoreへ渡す。

## 4. Control protocol

### 4.1 Browser → gateway: `RunRequestV1`

概念schemaは次とする。

```text
RunRequestV1
  protocolVersion = 1
  presetId
  run
    endTimeSeconds
    maximumTimeStepSeconds
    frameIntervalSteps
  initialCondition
    edits[]
      id
      kind = gaussian_depth
      centerUnit = [x, y, z]
      amplitudeMeters
      sigmaRadians
      massPolicy = preserve_global | allow_change
```

`presetId` はgateway起動時に設定したallowlistから選ぶ。requestからconfig path、output path、
binary path、environment variable、command argumentを直接指定できない。run overrideも
capability responseで公開した数値範囲だけを受理する。

### 4.2 Gateway → C++: `ControlRequestV1`

C++側へ新しいJSON dependencyを入れず、既存configと同系統のstrict `key = value` textを使う。

```text
control.format_version = 1
control.run_id = <gatewayが発行したopaque ID>
control.end_time_s = 86400
control.maximum_time_step_s = 600
control.frame_interval_steps = 10
control.frame_directory = <gatewayが決めたrun-root内path>
initial_edits.count = 1
initial_edits.0.kind = gaussian_depth
initial_edits.0.center_x = ...
initial_edits.0.center_y = ...
initial_edits.0.center_z = ...
initial_edits.0.amplitude_m = ...
initial_edits.0.sigma_rad = ...
initial_edits.0.mass_policy = preserve_global
```

gatewayは browser JSON を検証してこのfileを生成する。C++も独立に全key、型、range、重複、
有限値、unit vector、run ID、edit数上限、output pathを再検証する。最終的な実行可否はC++が決める。

### 4.3 C++ → gateway: event stream

`--control-request PATH --event-stream ndjson` を指定したmachine modeでは、stdoutの各行を
一つのJSON objectとする。human-readable metadataはrun directoryのfileへ書き、stdoutへ
混在させない。

```text
run.accepted
run.started
frame.ready { sequence, timeSeconds, step, relativePath, byteLength }
diagnostics.sample
run.completed
run.cancelled
run.failed { code, message }
```

eventには `protocolVersion`, `runId`, sequence numberを持たせ、未知eventはversion規約に従う。
stderrはgatewayがサイズ上限付きlogとして保存し、browserにはsanitized summaryだけを返す。
gatewayは起動時に `my_planet_sim --describe-control` を一度実行し、C++が返すprotocol version、
supported edit、numeric limitとgateway実装が一致しなければlive endpointを開始しない。

### 4.4 Gateway → browser: local HTTP API

```text
GET  /api/v1/capabilities
POST /api/v1/runs
GET  /api/v1/runs/{id}
GET  /api/v1/runs/{id}/events
GET  /api/v1/runs/{id}/frames/{sequence}
POST /api/v1/runs/{id}/cancel
GET  /api/v1/runs/{id}/bundle
```

event endpointはauthenticated `fetch` のstreaming responseを用い、browserは増分NDJSONを読む。
WebSocketは双方向通信が必要になるまで導入しない。commandは明示的なPOST、進捗は一方向stream、
field payloadはbinary GETに分け、再接続とcacheを単純にする。

### 4.5 `FrameV1`

frame metadataはevent/manifestに置き、payloadは明示的なlittle-endian `Float64` binaryとする。

```text
FrameV1 payload order
  depth[6*N*N]
  momentum_x[6*N*N]
  momentum_y[6*N*N]
  momentum_z[6*N*N]
```

metadataはschema version、mapping、`N`、flatten order、field descriptor、time、step、byte length、
config fingerprintを持つ。C++はtemporary fileへ書いてclose後にatomic renameし、gatewayは
`frame.ready` 後のfileだけを公開する。既存CSVとcheckpointは後方互換のまま維持する。

## 5. Interactive initial conditions

### 5.1 edit lifecycle

UIのrun stateは次とする。

```text
draft -> submitting -> running -> cancelling -> cancelled
                           └──────> completed
                           └──────> failed
```

`draft` では preset と edit list を変更できる。Runを押した時点でrequestをimmutable snapshot化し、
実行中のcontrol変更はactive runへ反映しない。変更する場合はcancelして新しいrunを作る。

### 5.2 click による Gaussian depth perturbation

3D/2Dのclickを共通unit vector `o` へ変換する。cell center `x_i` との球面角距離を

```text
d_i = acos(clamp(dot(o, x_i), -1, 1))
delta_h_i = A * exp(-d_i^2 / (2*sigma^2))
```

とし、C++でbenchmark initial stateへ加える。defaultは `massPolicy=preserve_global` とし、
cell area `S_i` を使って

```text
mean_delta = sum(S_i * delta_h_i) / sum(S_i)
h_i <- h_i + delta_h_i - mean_delta
```

とする。`allow_change` を選んだ場合は平均を引かず、追加したglobal massをprovenanceへ記録する。
どちらの場合も全cellで `h > depth_floor_m`、finite、momentum tangent invariantを検査し、違反時は
runを開始しない。

depth perturbation時のvelocityは変更しない。静止caseならpressure gradientにより実際の
gravity waveが励起される。既存の流れを持つpresetでは移流・回転・波が重なるため、UIにpresetと
initial diagnosticsを明示する。

### 5.3 preview と authoritative state

UIはclick直後の応答性のため、同じGaussian式で半透明の `Pending edit` overlayを表示してよい。
これは入力位置/幅のpreviewだけで、field valueやdiagnosticsとして扱わない。Run開始後にC++が
出力するframe 0を受信したらpreviewを消し、以後はC++ frameだけを表示する。

したがってPhase 3の波の表示経路は次になる。

```text
click -> RunRequestV1 -> gateway -> C++ initial edit -> C++ integration
      -> frame 0, 1, 2, ... -> gateway -> 3D/2D renderer
```

解析的ring animationは採用しない。見えている波はC++ stateの `depth_m` またはそのinitial差分である。

### 5.4 edit拡張点

Phase 3のproduction editは `gaussian_depth` 一種類に限定し、protocolはtagged unionとしておく。
将来のtracer paint、tangent velocity impulse、vortex、temperature perturbationは、それぞれ
unit、保存量、positivity、投影規約をC++ testで固定してからkindを追加する。

## 6. 3D/2D visualizer

### 6.1 表示用dataset

UI内部の `VisualDatasetV1` はsourceがoffline CSVでもlive FrameV1でも同じ形に正規化する。

```text
VisualDatasetV1
  schemaVersion = 1
  sourceKind = offline_csv | live_run
  grid
    topology = cubed_sphere
    mapping = equiangular_gnomonic_v1
    cellsPerPanel = N
    flattenOrder = panel_major_then_j_then_i
  frame
    timeSeconds, step, configFingerprint
  fields[]
    id, label, unit, kind, provenance, Float64Array values
```

現行 `tracer.csv` と `shallow_water.csv` のoffline importも維持する。CSVにはtime/stepがないため
unknownと表示し、推測しない。shallow-water frameから `speed_m_s=norm(momentum)/depth` と
initial frameとの差分 `depth_anomaly_m` をderived fieldとして提供する。

### 6.2 cubed-sphere geometry

ADR 0001のpanel basisと

```text
x = normalize(n + tan(alpha) * e_alpha + tan(beta) * e_beta)
```

をTypeScript側でも実装する。C++ sourceはimportせず、C++が出力したgolden coordinate fixtureで
全panel、seam、corner、forward/inverse mappingを照合する。

### 6.3 cubed-sphere grid overlay

3D/2D共通の `Grid` controlを置く。

- `Off`: grid lineを描かない。
- `Panel seams`: 異なるpanelを接続するedgeだけを太線で描く。
- `All cells`: panel seamに加え、全cell境界を細線で描く。

全cellの4辺からcanonical unique edge listを一度だけ作る。閉球面で `12*N*N` unique edges、
panel seam segmentが `12*N` になることをtestする。defaultは `Panel seams` とする。

3Dではunique-edge line geometryをfield surfaceよりわずかに外側へ描き、2Dでは同じedge IDを
D3で投影・antimeridian clipしてoverlay canvasへ描く。高解像度を理由にユーザー指定を黙って
無効化しない。

### 6.4 3D globe

- Three.jsの単一 `BufferGeometry` とtyped color/cell-id bufferを使う。
- dragで回転、wheel/pinchでzoom、resetで初期姿勢へ戻す。
- scalar colorはcell内で一定とし、panel seamへ補間artifactを作らない。
- Raycaster pickingをunit vectorとcell IDへ変換し、edit/inspect双方で共有する。
- selected cell、pending edit、gridをfieldとは別overlayで描く。
- WebGL context loss、resize、resource disposalを処理する。

### 6.5 2D global map

- D3 equirectangular projection + Canvas 2Dを使い、cellごとのSVG nodeは作らない。
- dragでpan、wheel/pinchでzoom、resetで全球へ戻す。
- base field、grid、selection/editをseparate canvasへcacheする。
- inverse projectionからunit vectorへ戻し、3Dと同じcell lookup/edit originを使う。
- antimeridianを横切るcell/edgeをclipし、地図全幅を横断するpolygonを作らない。
- 3D/2Dでfield、range、selection、edit、current frameを同期する。

### 6.6 field/diagnostics control

- depth、depth anomaly、speed、momentum componentを選択できる。
- sequential/diverging scale、auto/manual range、legendを持つ。
- constant/signed/extreme rangeでfinite colorを生成する。
- inspectorはclip前のraw C++ value、unit、panel/index、lon/lat、time/stepを示す。
- current mass/energy/enstrophy/AAMとinitial driftをeventから表示する。
- animationは受信済みframe内だけでplay/pause/scrubし、存在しないtimeを補間して科学値として
  表示しない。視覚補間を後で追加する場合は明示的にlabelする。

## 7. Gateway lifecycle と安全性

### 7.1 process管理

- gateway起動時にsimulator binary、preset directory、run rootを固定する。
- browser requestからこれらのpathを変更できない。
- `child_process.spawn(binary, args, {shell: false})` 相当で起動し、command stringを組み立てない。
- Phase 3は同時run数1をdefaultとし、超過requestは明示的に拒否またはqueueする。
- timeout、最大N、最大end time、最大frame数、最大output bytesを設定する。
- cancelはprocessへgraceful terminationを要求する。C++はstep境界のcancellation predicateで停止し、
  final eventを出す。期限内に終了しない場合だけgatewayが強制終了し、status遷移を記録する。
- gateway終了時にorphan processを残さない。

### 7.2 local HTTP防御

- `127.0.0.1`/`::1`だけにbindし、`0.0.0.0`をdefaultにしない。
- 起動時にrandom 256-bit session tokenを作り、URL fragment経由でUIへ渡してmemoryだけに保持する。
- API requestはauthorization header、Origin、Hostを検証する。
- request body、header、log、frame、run数にsize/rate limitを設ける。
- run/frame pathはopaque IDからserver側で解決し、`..`、absolute path、symlink escapeを拒否する。
- error responseにlocal absolute path、environment、raw command lineを含めない。
- UI assetに外部CDNやtelemetryを使わず、local simulation dataを外部送信しない。

### 7.3 provenance と再現性

run directoryに次を保存する。

- original `RunRequestV1`
- resolved preset/configと `ControlRequestV1`
- model version、Git revision、dirty flag、compiler、config fingerprint
- initial edit listとmass policy、initial diagnostics
- ordered event log、stderr log、frame manifest、binary frames
- terminal statusとexit code

UIからrun bundleをdownloadできるようにする。既存checkpoint/restartとの意味を変えず、Phase 3の
requestから同じrunをCLIで再実行できるcommand/documentationを残す。

## 8. 検証方針

### 8.1 C++ gate

- control requestのunknown/duplicate/missing key、version、range、非有限値、pathを拒否する。
- Gaussian editが回転した任意originで同じshapeとなり、seam/cornerに不連続を作らない。
- `preserve_global` のarea-weighted mass changeがrounding scale、`allow_change`が解析的/数値積分値と
  一致する。
- edit後もdepth positive、state finite、momentum tangent invariantを満たす。
- frame 0がedit適用後・積分前のstateとbitwise一致する。
- observerがinitial、指定interval、finalだけを順序通り通知し、solver resultを変えない。
- cancellation predicateをstep境界で検査し、graceful cancel時のstate/frame/eventが整合する。
- FrameV1 round-tripがbitwise一致し、truncated/wrong-endian/wrong-size/versionを拒否する。
- machine modeのstdoutはvalid NDJSONだけで、human modeの既存出力/CSV/checkpointを壊さない。
- rest + Gaussian depth perturbationで非zero momentumと伝播するdepth anomalyが生じ、mass driftが
  既存gate内、全depthがpositiveとなる。

### 8.2 protocol/gateway gate

- UI/gateway schemaのvalid/invalid fixtureを共有し、protocol version mismatchを明示する。
- fake simulatorでspawn args、shell無効、event ordering、partial line、stderr上限、nonzero exitをtest。
- arbitrary preset/path/argument、traversal、symlink escape、oversized request、bad origin/tokenを拒否する。
- cancel、timeout、gateway shutdownでchild processとstatusが収束する。
- frame.ready前のfile、temporary file、別runのfileを配信しない。
- reconnect時にlast sequence以降のeventを取得できる。

### 8.3 geometry/renderer/browser gate

- CSV/frameの長さ、field、flatten order、derived valueを検証する。
- `N=1,2,4,16` でcell数が `6*N*N`、unique edge数が `12*N*N`、seam数が `12*N`となる。
- 全panelのforward/inverse mappingがinteriorで `1e-12 rad` 以内に一致する。
- 3D/2Dで同じ既知cellを選択し、同じedit originをRunRequestへ格納する。
- grid controlの3 modeが両viewへ反映される。
- antimeridian、pole、panel seam/cornerで描画・picking・editが破綻しない。
- field/range/selection/current frameが両viewで同期する。
- keyboard、focus、accessible name、reduced motion、narrow viewportを検証する。

### 8.4 live end-to-end gate

小さい `N=4` rest presetを使い、実binaryで次を自動化する。

1. capabilitiesを取得する。
2. 2Dまたは3DでGaussian depth editを一つ作る。
3. runをsubmitし、accepted/started/frame 0/completedを受信する。
4. frame 0の最大depth位置がclick origin近傍であることを確認する。
5. 後続frameでdepth anomalyとmomentumが変化することを確認する。
6. UIの3D/2Dが同じframe/fieldを表示することを確認する。
7. 別runをcancelし、orphan process/partial published frameがないことを確認する。

### 8.5 performance gate

- `N=48` と `N=96` でframe encode、disk write、HTTP transfer、decode、GPU uploadを分けて計測する。
- field payloadはJSON/base64へ変換しない。
- latest frame優先時も受信済みframe順序とdrop policyを明示する。
- camera/map操作中にcellごとのDOM/object、frameごとのgeometry再生成をしない。
- `All cells` gridもedgeごとのscene objectを作らずbatch描画する。
- CIは構造指標/timeoutをhard gate、load/frame timeはreference machineのmedian/p95をreportする。

## 9. コミット規約

各コミットは一つのprotocol/core/UI責務と対応testを含む。C++変更は既存Phase 0–2 test、
GCC/Clang、ASan/UBSan、formatを維持する。Web変更はlint、format、typecheck、unit、buildを通す。
protocol変更はC++、gateway、UI fixtureを同じcommitで更新するか、後方互換な順序で分割する。
generated frame、run directory、`web/dist`、`node_modules`、browser reportはcommitしない。

## 10. コミット列

### P3.01 — control architecture ADR

```text
docs: define the interactive simulator control boundary
```

変更:

- ADR 0004にUI/gateway/C++境界、local-only前提、protocol versioningを固定。
- `RunRequestV1`, `ControlRequestV1`, `EventV1`, `FrameV1`を記述。
- synthetic wave previewを撤回し、C++ frameをauthoritativeとする。

受け入れ条件: process、data、security、failure、reproducibilityのownerが一意に決まる。

### P3.02 — C++ control request schema

```text
control: add versioned run control requests
```

変更:

- control request model、strict parser、validationをcoreへ追加。
- end time、maximum dt、frame interval、gateway-owned outputを既存configへ適用。
- canonical writeとfingerprint対象を定義。

テスト: round-trip、unknown/duplicate/missing key、version/range/path、既存config不変。

### P3.03 — C++ initial-condition edits

```text
dynamics: apply validated shallow-water initial edits
```

変更:

- tagged `InitialConditionEditV1` と `gaussian_depth` を追加。
- benchmark state構築後、積分前にeditを適用。
- preserve/allow mass policy、positivity、tangency、edit diagnosticsを追加。

テスト: origin回転、seam/corner、mass、positivity、複数edit順序、invalid edit。

### P3.04 — solver frame observer

```text
dynamics: publish shallow-water states to an observer
```

変更:

- driverへnon-owning frame observerを追加。
- step境界で検査するnon-owning cancellation predicateを追加。
- initial、interval、finalの通知規約を固定。
- observerなしの既存runにallocation/I/O/結果差を持ち込まない。

テスト: schedule、final重複防止、observer exception、cancel境界、bitwise result regression。

### P3.05 — binary frames と machine CLI

```text
io: stream versioned shallow-water run events
```

変更:

- FrameV1 little-endian writer/reader test helperとatomic publishを追加。
- `--control-request`, `--event-stream ndjson` CLIを追加。
- `--describe-control` とsignalからcancellation predicateへの安全なbridgeを追加。
- metadata、diagnostics、lifecycle/frame eventをmachine modeで出力。
- human CLI、CSV、checkpointの既存挙動を維持。

テスト: binary round-trip/corruption、NDJSON parse、failure event、legacy CLI回帰。

### P3.06 — 独立 Web workspace

```text
build(web): scaffold the visualizer workspace
```

変更:

- `web/apps/ui`, `web/apps/gateway`, `web/packages/protocol` とlockfileを追加。
- React/TypeScript/Vite、Node gateway、Vitest/Playwrightのscriptを追加。
- repository外source import禁止、Web専用CI、ignore ruleを追加。

テスト: C++なしで `npm ci`, lint, typecheck, unit, UI production build。

### P3.07 — shared Web protocol と mock client

```text
feat(web): define the simulation client protocol
```

変更:

- RunRequest/Event/Frame metadata validatorとrun state machineを追加。
- HTTP `SimulationClient` interfaceとdeterministic mockを追加。
- version mismatchとtyped errorをUIへ伝える。

テスト: valid/invalid contract、state transition、event sequence、mock cancellation。

### P3.08 — visual dataset と offline import

```text
feat(web): normalize simulator fields for visualization
```

変更:

- immutable `VisualDatasetV1` とfield/derived field modelを追加。
- FrameV1 decoder、現行transport/shallow-water CSV adapterを追加。
- field/range/legend/inspectorの基本controlを追加。

テスト: binary/CSV、row shuffle、missing/duplicate/nonfinite、depth anomaly/speed。

### P3.09 — cubed-sphere geometry と grid

```text
feat(web): reconstruct the cubed-sphere mesh and grid
```

変更:

- panel mapping、cell geometry、inverse lookup、unique edge/seam分類を追加。
- `Off`, `Panel seams`, `All cells` grid stateを追加。
- C++ golden coordinate fixtureを追加。

テスト: mapping、seam/corner、cell/edge/seam count、hit test。

### P3.10 — interactive 3D globe

```text
feat(web): render live fields on a 3D globe
```

変更:

- single BufferGeometry、color/cell ID、batched grid lineを追加。
- rotate/zoom/reset、inspect/edit picking、selection/pending-edit overlayを追加。
- resize、context loss、resource disposalを追加。

テスト: buffer finite/size、known picking、grid mode、unmount/context handling。

### P3.11 — interactive 2D global map

```text
feat(web): render live fields on a 2D global map
```

変更:

- D3 projection + layered Canvas、pan/zoom/resetを追加。
- field、grid、selection、pending edit、inverse-projection pickingを追加。
- antimeridian/pole clippingを追加。

テスト: projection round-trip、clipping、known picking、grid mode、resize。

### P3.12 — secure local gateway と run lifecycle

```text
feat(gateway): manage native simulator runs
```

変更:

- loopback HTTP、session token、origin/host check、capabilities/run endpointを追加。
- startup時に `--describe-control` を照合し、incompatible binaryをfail fastする。
- preset allowlist、isolated run directory、request/control translationを追加。
- configured binaryをshellなしでspawnし、status/log/limitを管理。

テスト: fake binary、exact argv、bad request/auth/path、nonzero exit、limit/timeout。

### P3.13 — event/frame delivery と cancellation

```text
feat(gateway): stream simulator frames and cancellation
```

変更:

- stdout NDJSON parser、ordered authenticated event streamを追加。
- published FrameV1 GET、bundle export、reconnect sequenceを追加。
- cancel/grace timeout/shutdown cleanupを追加。

テスト: partial/malformed event、frame isolation、cancel race、orphan/temporary file防止。

### P3.14 — interactive initial-condition editor

```text
feat(web): edit shallow-water initial conditions
```

変更:

- preset/run parameter form、draft edit list、undo/clearを追加。
- 3D/2D clickからGaussian amplitude/sigma/mass policyを設定。
- pending previewとrequest summary/validationを追加。

テスト: click座標、edit ordering、unit/range、immutable request、keyboard操作。

### P3.15 — live C++ simulation controller

```text
feat(web): run and inspect native simulations
```

変更:

- submit/progress/cancel/retry、event reconnectをUIへ接続。
- frame 0でpending previewを置換し、frame timeline/play/pause/scrubを追加。
- 3D/2Dのcurrent frame/field/selection/gridを同期。

E2E: mock gatewayでcomplete/fail/cancel/reconnect、last-good frame保持。

### P3.16 — diagnostics と reproducible run bundles

```text
feat(web): expose run diagnostics and provenance
```

変更:

- mass/energy/enstrophy/AAM、drift、time/step、CFLを表示。
- resolved config、initial edits、build metadata、event logをinspection可能にする。
- run bundle downloadとoffline reopenを追加。

テスト: metadata欠落、diagnostic sequence、bundle round-trip、raw value保持。

### P3.17 — live integration、安全性、性能 gate

```text
test: gate interactive simulator control end to end
```

変更:

- `N=4` native C++ live E2Eとcancel scenarioをCIへ追加。
- Chromium/Firefox/WebKitのoffline/control flowを追加。
- security regression、`N=48/96` frame/render benchmarkを追加。
- accessibility、responsive、WebGL failure/error UXを仕上げる。

受け入れ条件: click editがC++ frame 0へ一致し、後続stateがC++で進化する証拠を自動確認できる。

### P3.18 — 使用手順と Phase 3 検証報告

```text
docs: record Phase 3 interactive visualizer validation
```

変更:

- offline UI、gateway起動、binary/preset指定、run/cancel/export手順を記載。
- `docs/validation/phase-3.md` にprotocol、scientific gate、browser、security、performanceを記録。
- synthetic previewではなくC++ solver resultであること、local-only制約、既知gapを記録。

受け入れ条件: clean checkoutから小さいclick-wave runとbundle再生を再現できる。

## 11. 依存関係

```text
P3.01
  ├─ P3.02 ─ P3.03 ─ P3.04 ─ P3.05 ───────────────┐
  └─ P3.06 ─ P3.07 ─ P3.08 ─ P3.09 ─┬─ P3.10 ─┐ │
                                      └─ P3.11 ─┤ │
                                                ▼ ▼
                             P3.12 ─ P3.13 ─ P3.14 ─ P3.15
                                                       │
                                                       ▼
                                             P3.16 ─ P3.17 ─ P3.18
```

C++ machine modeとmock clientは並行実装できる。P3.12以降で実processを接続し、P3.15までは
mock E2Eも維持する。rendererは `SimulationClient` と直接結合せず `VisualDatasetV1` だけを読む。

## 12. 完了チェックリスト

- [x] UI/gateway/C++のprotocolと責務がADRで固定されている。
- [x] `web/` はC++なしでlint/test/build/offline previewできる。
- [x] live modeはloopback gatewayから設定済みnative C++ binaryだけを起動する。
- [x] browserからarbitrary command/path/environmentを指定できない。
- [x] 3D/2D clickからGaussian depth initial editを作成・取消できる。
- [x] C++ frame 0がedit後のauthoritative initial stateとして表示される。
- [x] 表示されるwaveはsynthetic animationではなくC++が積分したfieldである。
- [x] run/progress/cancel/fail/reconnectが一貫したstate machineに従う。
- [x] 3D球体と2D全球地図が同じfield/frame/cellを表示する。
- [x] panel seam/全cell gridを両viewで表示・非表示にできる。
- [x] mass policy、positivity、tangency、conservation gateが通る。
- [x] request、resolved config、edit、build metadata、diagnostics、framesをbundle化できる。
- [x] existing Phase 0–2 CLI/CSV/checkpoint/testが回帰しない。
- [ ] native live E2E、3 browser、security、`N=48/96` performance結果が報告されている。

最終項目のうち native live E2E、security、performance は
[Phase 3 検証報告](validation/phase-3.md)で確認済みである。3 browser matrix と
gatewayからbuilt UIを直接配信する部分は既知 gap として次Phaseへ残す。browser entrypointは
URL fragmentでloopback gatewayとsession tokenを受け取るlive modeに対応している。

## 13. 技術資料

- [Node.js child process](https://nodejs.org/api/child_process.html): shellを介さない非同期process起動、stdout stream、termination。
- [Node.js HTTP](https://nodejs.org/api/http.html): local HTTP serverとstreaming response。
- [Vite: Building for Production](https://vite.dev/guide/build): UIの独立static build。
- [Three.js: BufferGeometry](https://threejs.org/docs/pages/BufferGeometry.html): typed geometry attributeとbatched mesh。
- [Three.js: Raycaster](https://threejs.org/docs/pages/Raycaster.html): 3D pointer picking。
- [D3: Geographic projections](https://d3js.org/d3-geo/projection): spherical projection、clipping、inverse projection。
- [Playwright: Test](https://playwright.dev/docs/intro): cross-browser/offline/live integration test。

# MyPlanetSim 開発・検証計画

## 1. 目的と前提

MyPlanetSim の目的は、観測値が十分にない架空惑星にも適用できる、再現可能で検証可能な全球大気シミュレーターを C++ で構築することである。最終的な基本構成は次のとおりとする。

- 水平格子: 6 面の equiangular gnomonic cubed-sphere
- 鉛直座標: 地表付近で地形追随し、上層で等圧面へ移る hybrid sigma-pressure 座標
- 初期力学系: 乾燥・静水圧 primitive equations
- 離散化の出発点: 保存型有限体積法、倍精度、明示時間積分
- 惑星依存量: 半径、自転角速度、重力、気体定数、比熱、基準気圧、日長などを設定値として分離

「動いて見える」ことではなく、単純な系で解析解、収束率、保存則を確認してから複雑さを一つずつ増やす。各段階は前段階の診断とテストを再利用し、完了条件を満たすまで次へ進まない。

初期スコープには湿潤物理、雲、詳細放射、化学、海洋、非静水圧運動を含めない。ただし後から追加できるよう、力学・物理・惑星設定・I/O の境界は分離する。

## 2. 数値設計の初期方針

### 2.1 二段階の数値方式

まず比較基準となる堅牢な方式を作り、その後に長時間積分向けの方式を評価する。

1. 基準実装: cell-centered finite volume、MUSCL 型再構築、Rusanov/HLL 型数値流束、SSP-RK3。拡散は大きいが、保存性と実装検証を優先する。
2. 目標実装: 面法線質量流束を共有し、質量・トレーサ整合性、エネルギー、渦位の性質を改善した mimetic/compatible な方式。Phase 2 の比較結果を ADR（Architecture Decision Record）にして採否を決める。

最初から高次・複雑な方式だけを実装すると、格子計量、面接続、再構築、時間積分の誤差を切り分けにくい。低次の基準実装は、後の高次方式に対する独立した比較器として残す。

### 2.2 cubed-sphere

- 各面を論理座標 `(alpha, beta)` の構造格子として扱う。
- 座標写像、Jacobian、計量テンソル、面積、辺長、面法線を geometry モジュールで一度だけ定義する。
- 面をまたぐ物理ベクトルは 3 次元 Cartesian 接ベクトルを共通表現にして変換し、局所成分の単純コピーを禁止する。
- 共有辺の数値流束は両面で同一値・逆符号になるよう一か所で計算または同一規約で参照し、全球質量保存を構造的に保証する。
- 8 個の cube corner は通常辺と別のテスト対象とし、corner を含む stencil の向きと重複更新を検査する。

### 2.3 hybrid sigma-pressure 座標

層境界（half level）の圧力を

```text
p[k+1/2] = A[k+1/2] + B[k+1/2] * ps
```

とする。本プロジェクトでは `A` の単位を Pa、`B` を無次元、`ps` を地表気圧と定義し、別流儀の `a * p0` との混同を防ぐ。上端では `B = 0`、下端では `A = 0, B = 1` を基本とし、全層で圧力が単調になることを設定読込時に検査する。

初期の鉛直配置は Lorenz grid とし、水平風、温位、トレーサを full level、圧力と鉛直質量流束を half level に置く。予報量は層質量、水平運動量、質量重み付き温位、質量重み付きトレーサとし、`ps` は鉛直積算された質量連続式から更新する。鉛直差分は、Simmons–Burridge 型のエネルギー・角運動量整合性を参照して設計する。

### 2.4 診断量と誤差規約

全実験が同じ診断 API から次を出力する。

- `L1`, `L2`, `Linf` 誤差と観測収束次数
- 全球または領域の質量、トレーサ質量、全エネルギー、可能ならポテンシャルエンストロフィー・軸角運動量
- 最小・最大値、負の層厚/圧力/密度、NaN/Inf の有無
- CFL 数、時間刻み、解像度、惑星定数、ビルド情報
- 面別誤差と seam/corner 近傍誤差

解析解がある場合、解像度 `h`, `h/2`, `h/4` の少なくとも 3 点で

```text
p_observed = log(error(h) / error(h/2)) / log(2)
```

を求める。原則として `p_observed >= p_design - 0.3` を合格の目安とするが、不連続解や limiter 作動時には適用せず、問題ごとの基準を記録する。回帰テストは画像比較ではなく数値診断を主とし、コンパイラ差を考慮した許容誤差を明記する。

## 3. 開発ロードマップ

### Phase 0 — 科学・ソフトウェア基盤

コミット単位の実装順、変更対象、受け入れ条件は [Phase 0 実装計画](phase-0-plan.md) に定める。
[検証報告](validation/phase-0.md)に示すローカルゲートは完了済みである。

実装するもの:

- C++20、CMake、テストランナー、警告設定、Debug/Release、sanitizer 構成
- SI 単位系、強い型または明示的な単位規約、惑星パラメータ、実験設定
- field/array、境界条件、時間積分器、診断、checkpoint の最小インターフェース
- 実験ごとの設定・出力 metadata と乱数 seed の固定
- 数式、変数配置、保存量、符号規約を記した設計文書

検証するもの:

- clean build と unit test が GCC/Clang の双方で通ること
- Debug で ASan/UBSan が通ること
- 定数場、ゼロ場、配列境界、設定値の単位・範囲エラーを検出できること
- 同じ入力・同じスレッド数で診断値が再現すること
- restart 前後で状態と診断値が許容誤差内で一致すること

完了ゲート: 空の時間積分ループではなく、製造解を用いた小さな ODE が設計次数で収束し、CI で上記項目を自動実行できる。

### Phase 1 — cubed-sphere 幾何と球面輸送

コミット単位の実装順、数値方式、変更対象、受け入れ条件は [Phase 1 実装計画](phase-1-plan.md) に定める。
[検証報告](validation/phase-1.md)に示すローカルゲートは完了済みである。

実装するもの:

- 6 面の格子生成、面 adjacency/orientation、halo exchange
- Cartesian と各面の局所座標間の scalar/vector 変換
- 面積・辺・法線・gradient/divergence/curl/Laplacian の離散演算子
- 球面上の prescribed-wind conservative tracer transport

検証するもの:

- 全球セル面積和が `4*pi*R^2` に収束し、全セル面積・辺長が正であること
- 座標の往復変換、共有辺位置、接ベクトル、法線向きが許容誤差内で一致すること
- 定数の gradient と curl(gradient) がゼロ、閉球面上の divergence 積分がゼロになること
- solid-body rotation を軸方向を変えて実行し、面境界を横切っても誤差が局在しないこと
- Williamson test 1 と Lauritzen et al. の非発散/発散変形流で収束、filament、mixing、極値、質量を評価すること
- seam/corner セルの誤差分布が面内部に比べて異常に増幅しないこと

完了ゲート: 全球質量保存を保ち、少なくとも二次の滑らかな輸送試験で設計収束率を示す。回転軸を変えた結果が同等の誤差範囲に入る。

### Phase 2 — cubed-sphere 全球 shallow-water

コミット単位の実装順、基準 finite-volume 法と compatible 候補の比較方法、
受け入れ条件は [Phase 2 実装計画](phase-2-plan.md) に定める。実装は完了し、
結果は [Phase 2 検証報告](validation/phase-2.md) に、採用する水平離散化は
[ADR 0003](adr/0003-shallow-water-horizontal-discretization.md) に記録した。

実装するもの:

- 球面計量と惑星自転を含む shallow-water 力学
- 面間で共有される質量・運動量流束
- 明示的な diffusion/hyperdiffusion とその budget
- 基準実装と conservation 改良候補を切替可能にする実験インターフェース

検証するもの:

- Williamson test 2（定常地衡流）と test 6（Rossby–Haurwitz wave）。山岳を含む test 5 は Phase 6 へ回す
- Galewsky–Scott–Polvani の不安定 jet を高解像度または公開参照解と比較
- 質量、全エネルギー、ポテンシャルエンストロフィー、軸角運動量の時系列
- grid imprinting、4-grid-wave、seam 反射、cube corner noise の有無
- 解像度・時間刻みを別々に細分化し、空間誤差と時間誤差を分離
- 基準方式と改良方式を同一条件で比較し、精度・保存性・コストを ADR に記録

完了ゲート: 標準テストの数値解・収束曲線・保存量が再生成可能で、採用する水平離散化を根拠付きで固定できる。地形項そのものは Phase 6 まで正式機能にしない。

### Phase 3 — C++ 制御付き独立 Web 可視化 UI

アーキテクチャ、data contract、画面・操作仕様、コミット単位の実装順、受け入れ条件は
[Phase 3 実装計画](phase-3-plan.md) に定める。

実装するもの:

- C++ source/headerへ直接依存しない `web/` の TypeScript UI、local control gateway、独立 CI
- version付きrun/event/frame protocolと、native C++ processの安全な開始・停止
- 3D/2D clickによるGaussian depth初期条件editとC++ shallow-water再積分
- Phase 1/2 の CSV snapshot とlive binary frameを取り込む表示 data model
- cubed-sphere field を表示する interactive 3D globe と 2D 全球地図
- 3D/2D の双方で panel seam/全 cell edge を切り替える cubed-sphere grid overlay
- field/color/selection を共有する inspector と view control
- C++が出力するinitial/intermediate/final frame、diagnostics、provenance、run bundle

検証するもの:

- panel、seam、corner、antimeridianを含むgeometry・投影・pickingの一致
- grid overlayのunique edge数と3D/2D表示modeの一致
- 3D/2D間のfield、range、selection、initial edit、current C++ frameの同期
- initial editのmass policy、positivity、tangencyと、C++で伝播するgravity wave
- protocol version、process lifecycle、cancel、path/command injection防止
- invalid CSV/frame、非有限値、欠落/重複cell、WebGL failureの明示的な処理
- keyboard、reduced motion、responsive layoutの基本操作
- `N=48/96` datasetのload、buffer/node数、frame timeとbrowser互換性

完了ゲート: `web/` はC++なしでもoffline build/test/previewでき、live modeではloopback gateway
から設定済みnative C++だけを制御できる。3D/2D clickで作ったdepth perturbationがC++のframe 0へ
一致し、後続frameで実際のshallow-water waveとして伝播することを自動検証できる。

### Phase 4 — 鉛直 1D と hybrid sigma-pressure

数値方式、state、C++ column core のコミット列は
[Phase 4 実装計画](phase-4-plan.md) に定める。鉛直座標・mass flux の規約は
[ADR 0005](adr/0005-hybrid-vertical-coordinate-and-column-state.md) に記録する。
C++ column core の検証を先に完了し、その後で C++ が出力する `column_profile.csv` を読む
offline profile visualizer を更新する。live 実行と 3D 連動は Phase 5 で更新する。

実装するもの:

- `A/B` 係数、full/half level、層質量 `Delta p/g`
- 静水圧積分、温位/温度/Exner 関数変換
- 固定 `A/B` 上の鉛直移流と鉛直 mass flux（conservative remap は導入しない）
- 上端・下端境界条件と鉛直診断
- C++ の最終 profile CSV を再計算なしで表示する offline column view

検証するもの:

- pure pressure (`B=0`) と pure sigma に近い極限がそれぞれ既知の式へ一致すること
- 任意の妥当な `ps` で圧力境界が単調、層質量が正、層質量和が `(ps-p_top)/g` と一致すること
- 等温・乾燥断熱の静水圧柱で geopotential が解析解に一致し収束すること
- 静止柱の state が変化せず、鉛直 mass flux が質量・定数トレーサを保存すること
- `ps` 変化時にも定数場保持、positivity、上下境界 flux が整合すること

完了ゲート: まず鉛直演算子だけの unit/convergence テストが通り、水平力学から独立して質量・
熱力学 budget を閉じられる。その後、strict CSV parser と column profile view の Web test/build が
通り、既存 shallow-water visualizer が回帰しない。

### Phase 5 — 地形なし乾燥 3D 静水圧力学コア

state、数値結合、C++ gate、最後に visualizer を更新するコミット順は
[Phase 5 実装計画](phase-5-plan.md) に定める。3D state と水平・鉛直 mass coupling は
[ADR 0006](adr/0006-dry-hydrostatic-state-and-coupling.md) に固定する。

実装するもの:

- cubed-sphere と hybrid 座標を直積した 3D field
- 乾燥 hydrostatic primitive equations、熱力学式、3D tracer 輸送
- surface pressure、coordinate-relative pressure velocity、hydrostatic geopotential の診断
- 水平・鉛直 flux の整合、外部重力波に対する時間積分
- 3D state 確定後の version 付き frame/data contract と gateway 対応
- 既存 globe/map の model-level 表示と、選択した水平 cell の鉛直 profile panel

検証するもの:

- 静止等温大気と静止成層大気を維持すること
- 3D solid-body tracer transport、DCMIP の 3D deformational flow / Hadley-like circulation
- UMJS14 の平坦地表・定常基本場を保ち、摂動時の baroclinic wave を参照解と比較
- 3D の乾燥質量、トレーサ質量、全エネルギー、軸角運動量 budget
- 水平・鉛直解像度を独立に上げた収束、鉛直モード・層間 decoupling の検査
- 長時間積分で上端反射、負圧、温度逸走、面境界ノイズがないこと
- globe/map の field・level・time・選択 cell と鉛直 profile が同じ C++ frame に一致すること
- 3D frame の size、転送、decode、描画 budget と旧 shallow-water UI の回帰

完了ゲート: 地形なしの標準 3D テストが複数解像度で再現し、shallow-water 段階の水平誤差と
新しい鉛直誤差を区別できる。C++ gate の完了後、同じ 3D state を既存 visualizer の level map と
選択 column profile で表示し、native run から UI まで自動検証できる。

### Phase 6 — 地形と下部境界

固定地形の data contract、数値式、C++ gate、コミット順は
[Phase 6 実装計画](phase-6-plan.md) に定める。地形は予報 state へ追加せず、判断は
[ADR 0008](adr/0008-fixed-orography-and-lower-boundary.md) に記録する。Phase 5 の既知の定量的
validation gap を閉じてから実装を開始する。

実装するもの:

- geopotential/orography field と地表境界条件
- 傾いた hybrid 面上の水平圧力傾度項
- 地形データ補間、平滑化、解像度別フィルタ
- 山岳トルクとエネルギー交換の診断

検証するもの:

- 平坦地形が Phase 5 と数値的に同一であること
- DCMIP 2-0-x 型の、山岳上で静止する静水圧大気に偽風が発生しないこと
- 山の高さ・幅・鉛直解像度に対する pressure-gradient error の収束
- Williamson test 5 の山岳付き shallow-water を回帰試験として再利用
- 低振幅の線形山岳波について波長・位相・運動量 flux を理論または高解像度解と比較
- 地形が panel seam/corner を横切る位置へ回転しても結果が同等であること
- Jablonowski–Williamson の smooth topography variant を比較すること

完了ゲート: 静止大気の偽風と pressure-gradient error が解像度とともに減少し、山岳位置を変えても保存量と主要診断が同等になる。

### Phase 7 — 理想化乾燥物理と気候統計

forcing、physics--dynamics coupling、気候統計、production matrix の実装順と境界は
[Phase 7 実装計画](phase-7-plan.md) に定める。Phase 5/6 の production validation と
ADR 0006 の SSP-RK3 実装一致を閉じてから開始し、汎用 physics framework や Web 拡張は含めない。
限定した physics、coupling、online statistics と CI gate は実装済みである。時間のかかる Phase 5/6
entry gate と Phase 7 の 1200 日 production matrix は保留されており、結果と未完了範囲は
[Phase 7 検証報告](validation/phase-7.md) に記録する。

実装するもの:

- Held–Suarez 型 Newtonian relaxation と Rayleigh friction
- physics–dynamics coupling と加熱・摩擦 budget
- ensemble/時間平均、帯状平均、spin-up 除外期間の解析

検証するもの:

- Held–Suarez の帯状平均温度・東西風、Hadley cell、eddy statistics を公開結果と比較
- 時間刻み、解像度、拡散に対する気候統計の限定的な一要因感度
- 強制・散逸を含むエネルギー budget が閉じること
- restart が瞬時値だけでなく長時間統計を変えないこと
- 少なくとも複数 seed/微小摂動で内部変動の幅を見積もること

完了ゲート: 単一 snapshot ではなく平均値とばらつきを伴うベンチマーク報告を自動生成できる。

### Phase 8 — 架空惑星化

惑星・軌道、解析的 benchmark forcing、Earth を含む固定 surface geography、fractional land--ocean
surface の詳細設計と実装順は [Phase 8 実装計画](phase-8-plan.md) に定める。海洋は深さや輸送を持たず、
陸海で異なる面積熱容量だけを持つ。海岸 cell は二値化せず `0<f_land<1` を保持する。主要な判断は
[ADR 0010](adr/0010-planetary-forcing-and-fractional-surface.md) に記録する。

実装するもの:

- 惑星半径 `R`、重力 `g`、自転 `Omega`、気体定数 `Rd`、`cp`、基準気圧、熱的強制、軌道/恒星日の設定化
- Earth の elevation と fractional land field、および `0<f_land<1` の海岸 cell
- 深さ・輸送を持たず、陸海で面積熱容量だけが異なる surface energy reservoir
- Earth-like、slow rotator、rapid rotator、tidally locked の実験 preset
- 次元量と無次元量（Rossby 数、Burger 数、変形半径、放射/回転時間尺度など）の metadata

検証するもの:

- Earth-like preset が Phase 7 の基準結果を変えないこと
- 単位だけを整合的に変換した相似問題が同じ無次元解を与えること
- `Omega -> 0`、弱い強制、半径変更などの極限で挙動が物理的期待と一致すること
- パラメータ sweep で CFL、層厚、温度、風速が安全範囲を外れたとき明示的に停止すること
- tidally locked terrestrial planet は Heng et al. の乾燥ベンチマーク、将来の物理追加後は THAI protocol と比較すること

完了ゲート: コード変更なしに惑星設定を切替でき、各 preset が設定・診断・再現手順を伴う。

### Phase 9 — 正しさの回復、性能、平衡ベンチマーク

**完了（2026-09-06）。** 結果と既知の DCMIP pressure-gradient gap は
[Phase 9 validation](validation/phase-9.md) に記録した。

Phase 8 完了時点の全体レビューが見つけた欠陥と、Phase 5--7 検証報告および ADR 0009 に記録済みの
production gap を、新しい物理を追加せずに閉じる。詳細な設計とコミット順は
[Phase 9 実装計画](phase-9-plan.md) に定める。

実装するもの:

- 球面ベクトル Laplacian の符号と曲率項の修正（[ADR 0011](adr/0011-vector-operators-and-time-step-normalization.md)）
- 乾燥・shallow-water・輸送の CFL 定義を cell 単位へ統一し、地表 reservoir に時間刻み制約を課す
- 乾燥コアへの水平拡散の配線。ADR 0009 の 6-run matrix の blocker を外す
- 静的格子ジオメトリ cache、RHS workspace、observer 間隔化と性能ゲート（[ADR 0012](adr/0012-static-grid-cache-and-performance-gates.md)）
- 平衡 benchmark 初期場と DCMIP 座標の修復（[ADR 0013](adr/0013-balanced-benchmark-initial-states.md)）

検証するもの:

- `steady` を名乗る preset の初期残差が解像度とともに減少すること
- 性能コミットの bit-exact 性と、数値を動かす唯一のコミットの登録許容差
- 暗黙 Rusanov 散逸と陽的拡散を同じ単位で比較した測定

完了ゲート: production matrix が MPI/OpenMP なしで実行可能であり、1 run の pilot コストが実測されている。

### Phase 10 — 重力波 semi-implicit 時間積分

乾燥静水圧コアで CFL を支配する水平 Lamb 波と高速な内部重力波を reference-linear
semi-implicit 法で解き、既存 SSP-RK3 を比較基準として残したまま通常 1200--1800 s の時間刻みを
可能にする。詳細な方式、段階的実装、精度・地形・長時間・性能 gate は
[Phase 10 実装計画](phase-10-plan.md) に定める。

実装するもの:

- full nonlinear RHS と保存的 reference-linear fast-wave operator の add--subtract 分割
- iterative Crank--Nicolson、鉛直重力波 mode 分解、2-D Helmholtz solve
- matrix-free GMRES と preconditioner、solver residual・iteration・retry 診断
- explicit/implicit/vertical/surface CFL の分離と conditional config
- flat wave、DCMIP terrain rest、mountain wave、baroclinic、Held--Suarez の段階的 gate

完了ゲート: `N=12`, `K=20` の Held--Suarez で median accepted step が 1200 s 以上、同じモデル期間の
explicit 200 s run より 3 倍以上速く、1 本の 1200 日 run が保存量・solver 診断付きで完走する。

状況: **完了（2026-09-07）**。反復停止則を残差許容値から固定反復数へ変更した結果、Held--Suarez は
通常1800 sで進み、explicit 200 s の5.78倍で走る。1200日pilotは20分57秒で完走した。day 412の
一時的な強風では陽的移流CFLに従う短縮とretryが生じたが直後に1800 sへ復帰し、solver反復の
経時悪化、非有限値、質量driftはなかった。実測値と完了判定は
[Phase 10 検証報告](validation/phase-10.md) を参照。

### Phase 11 — 乾燥大気の灰色放射

短波・長波を分けた灰色放射を各鉛直 column で解き、大気の吸収・熱放射と地表の熱収支を結合する。
数値方式、既存 surface/半陰解法との接続、実装順序、検証ゲートは
[Phase 11 実装計画](phase-11-plan.md) と
[ADR 0017](adr/0017-gray-radiation-column-and-energy-attribution.md) に定める。状態は実装中。

- 短波の直接吸収・地表反射、長波の上下二流輸送と地表への下向き放射
- 同一の界面 flux による各層の加熱率と地表蓄熱の更新
- 単独 column の解析解・保存・収束を検証した後、乾燥 3D core へ結合
- 放射・顕熱の時間刻み制約、30 分刻みの再検証、restart と放射 budget の診断

完了ゲート: column の解析・保存・収束試験、全球結合と長時間 pilot、旧経路の回帰を通過し、
放射の局所保存と既存力学のエネルギー残差を区別して報告できる。湿潤・雲・散乱・乾燥対流調節と
詳細放射は後続とする。

### Phase 12 — 乾燥対流調節・境界層

既存モデルの調査、数値方式、実装順序と受け入れ条件は
[Phase 12 実装計画](phase-12-plan.md) に定める（計画段階）。

実装するもの:

- 固定質量で柱の enthalpy を保存する乾燥対流調節
- 安定度・粗度に依存する地表顕熱/応力と、境界層の鉛直乱流拡散
- 大気・地表を結合した陰的 column solve、散逸熱、過程ごとの収支
- 既存灰色放射との責務分離、step 全体の retry、診断・checkpoint 再開

検証するもの:

- 調節後安定性・保存・冪等性、拡散解析解・tracer 最大値原理、地表交換の解析的緩和
- 地表付近を細分化した鉛直格子、結合の時間精度、物理過程の有無を分けた全球比較
- 30日比較と1200日 pilot、局所過程収支と既存全系エネルギー残差の区別

完了ゲート: 独立 column と全球結合の保存・安定性・restart を検証し、鉛直解像度・時間刻みの
適用範囲を登録できる。Phase 11 の全系残差は別途追跡し、対流導入だけで気候平衡とは判定しない。

### Phase 13 — 蒸発・水蒸気・降水

設計判断、保存量、実装順序と受け入れ条件は
[Phase 13 実装計画](phase-13-plan.md) に定める（計画段階）。

- 希薄水近似の予報水蒸気と複数 tracer、保存的・非負な輸送
- 境界層と地表に結合した蒸発・露・潜熱、陸面 bucket と流出の水収支
- operator splitting による瞬時飽和調節と即時降水
- Simple Betts–Miller / Frierson 型の深い・浅い対流
- 過程別の水・エネルギー台帳、restart、独立 column から全球 pilot への段階検証

完了ゲート: 多成分輸送・独立湿潤 column・両 integrator の結合で保存、非負性、再開を検証し、
30日比較と1200日 aquaplanet pilot の適用範囲を報告する。水蒸気の放射効果、雲・氷・降水再蒸発、
既存全系エネルギー残差の解消はこの完了条件に含めない。

### Phase 14 以後 — 物理拡張と並列性能

候補:

- 水蒸気の放射効果・簡易雲・氷過程、その後に詳細放射
- 非静水圧方程式（必要な惑星・解像度が明確になった場合）
- OpenMP/MPI、並列 halo exchange、並列 I/O、accelerator 対応
- advective CFL が実測上の次の支配項になった場合の semi-Lagrangian または別の長時間刻み輸送

各物理過程は単一 column test と process budget を先に作る。並列化は serial 結果との許容誤差内一致、
領域分割不変性、strong/weak scaling を検証し、科学的回帰試験を通過させる。

## 4. 横断的な品質基準

### 4.1 テスト階層

- Unit: 座標変換、計量、EOS、flux、境界条件、`A/B` 係数
- Property: 定数保持、正値性、回転軸不変性、共有 flux の反対称性
- Convergence: 解析解/製造解を 3 解像度以上で評価
- Regression: 標準ケースのスカラー診断と少数の場データ checksum
- Long-run: 保存量 drift、統計平衡、restart
- Performance: cell-updates/s、メモリ/cell、通信比率。正しさのゲートとは分離

### 4.2 実験の再現性

各出力に Git commit、dirty flag、compiler、浮動小数点設定、設定ファイル全文、格子、`A/B`、時間刻み、フィルタ係数を保存する。参照結果は生成スクリプトと出典を持ち、巨大バイナリを無条件に Git へ入れない。

### 4.3 段階を進めない条件

以下が一つでも残る場合は、次の物理機能を追加しない。

- 質量不整合、NaN/Inf、負の圧力/層厚を「フィルタで隠している」
- seam/corner だけに説明不能な誤差が集中する
- 解像度を上げても誤差が減らない
- 保存量 drift の発生項を budget で特定できない
- 初期条件、設定、参照値、実行コマンドが記録されていない

## 5. 直近の実装順

Phase 9 までの基盤修正と性能改善は完了した。次は、新しい物理を積む前に explicit Lamb CFL による
200 s 制約を外し、長時間 climate integration の step 数そのものを減らす。順序と根拠は
[Phase 10 実装計画](phase-10-plan.md) に定める。

1. ADR 0016 と explicit baseline測定で fast/slow split、保存、solver、精度・性能閾値を固定する。
2. total RHSを変えずに componentとCourant診断を分離する。
3. GMRES/Helmholtz solverを線形 shallow-water fixtureで独立に検証する。
4. hydrostatic reference column、鉛直mode、保存的 tangent-linear fast operatorを実装する。
5. external modeから始め、必要なinternal modesへ iterative Crank--Nicolson solveを拡張する。
6. flat、terrain、baroclinic、Held--Suarezの順に1200--1800 sを検証する。
7. 30日比較で3倍以上のspeedupを確認してから、1本の1200日pilotを実行する。
8. 放射・湿潤・非静水圧・本格的並列化はPhase 11以後とする。

## 6. 調査資料と計画への反映

- C. Ronchi, R. Iacono, P. S. Paolucci (1996), [The “Cubed Sphere”: A New Method for the Solution of Partial Differential Equations in Spherical Geometry](https://doi.org/10.1006/jcph.1996.0047): 6 面分割、great-circle coordinate、面結合という cubed-sphere の基礎。
- D. L. Williamson et al. (1992), [A Standard Test Set for Numerical Approximations to the Shallow Water Equations in Spherical Geometry](https://doi.org/10.1016/S0021-9991(05)80016-6): 全球 shallow-water の標準テスト群。
- P. H. Lauritzen et al. (2012), [A standard test case suite for two-dimensional linear transport on the sphere](https://doi.org/10.5194/gmd-5-887-2012): 収束、filament preservation、mixing、非発散/発散輸送の診断。
- J. Galewsky, R. K. Scott, L. M. Polvani (2004), [An initial-value problem for testing numerical models of the global shallow-water equations](https://doi.org/10.3402/tellusa.v56i5.14436): 高速重力波と遅い不安定渦発達を同時に試す jet case。
- J. Thuburn, C. J. Cotter, T. Dubos (2014), [A mimetic, semi-implicit, forward-in-time, finite volume shallow water model](https://doi.org/10.5194/gmd-7-909-2014): cubed-sphere 上の質量・エネルギー・PV 整合性と mimetic 演算子の評価観点。
- A. J. Simmons, D. M. Burridge (1981), [An Energy and Angular-Momentum Conserving Vertical Finite-Difference Scheme and Hybrid Vertical Coordinates](https://doi.org/10.1175/1520-0493(1981)109%3C0758:AEAAMC%3E2.0.CO;2): hybrid 座標と鉛直差分の保存性。
- ECMWF, [Vertical coordinate](https://confluence.ecmwf.int/plugins/viewsource/viewpagesrc.action?pageId=85406187): 実運用される half-level の `p = A + B ps` 定義と上下層の性質。
- P. A. Ullrich et al. (2014), [A proposed baroclinic wave test case for deep- and shallow-atmosphere dynamical cores](https://doi.org/10.1002/qj.2241): constant surface pressure と zero surface geopotential を持つ Phase 5 の平坦地表 baroclinic test。
- C. Jablonowski, D. L. Williamson (2006), [A baroclinic instability test case for atmospheric model dynamical cores](https://doi.org/10.1256/qj.06.12): smooth surface geopotential を含む Phase 6 の baroclinic test。
- DCMIP (2012), [Dynamical Core Model Intercomparison Project test case document](https://public.websites.umich.edu/~cjablono/DCMIP-2012_TestCaseDocument_v1.7.pdf): 3D 輸送、静水圧静止大気、地形性 pressure-gradient error などの段階的テスト。
- I. M. Held, M. J. Suarez (1994), [A Proposal for the Intercomparison of the Dynamical Cores of Atmospheric General Circulation Models](https://www.gfdl.noaa.gov/bibliography/related_files/ih9401.pdf): 乾燥力学コアの長時間統計ベンチマーク。
- K. Heng, K. Menou, P. J. Phillipps (2011), [Atmospheric circulation of tidally locked exoplanets: a suite of benchmark tests for dynamical solvers](https://doi.org/10.1111/j.1365-2966.2011.18315.x): 架空惑星・同期回転惑星への拡張時の比較ケース。
- E. M. Sergeev et al. (2023), [Simulations of idealised 3D atmospheric flows on terrestrial planets using LFRic-Atmosphere](https://doi.org/10.5194/gmd-16-5601-2023): cubed-sphere を用いた地球型惑星 3D 力学と惑星ベンチマークの現代的な実例。

この計画の閾値は初期値である。各 Phase の最初に、選んだ離散化次数と標準ケースの公開参照値に基づいて数値許容値を確定し、途中で緩和する場合は理由を ADR と検証記録に残す。

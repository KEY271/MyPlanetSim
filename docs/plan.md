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

実装するもの:

- 球面計量と惑星自転を含む shallow-water 力学
- 面間で共有される質量・運動量流束
- 明示的な diffusion/hyperdiffusion とその budget
- 基準実装と conservation 改良候補を切替可能にする実験インターフェース

検証するもの:

- Williamson test 2（定常地衡流）と test 6（Rossby–Haurwitz wave）。山岳を含む test 5 は Phase 5 へ回す
- Galewsky–Scott–Polvani の不安定 jet を高解像度または公開参照解と比較
- 質量、全エネルギー、ポテンシャルエンストロフィー、軸角運動量の時系列
- grid imprinting、4-grid-wave、seam 反射、cube corner noise の有無
- 解像度・時間刻みを別々に細分化し、空間誤差と時間誤差を分離
- 基準方式と改良方式を同一条件で比較し、精度・保存性・コストを ADR に記録

完了ゲート: 標準テストの数値解・収束曲線・保存量が再生成可能で、採用する水平離散化を根拠付きで固定できる。地形項そのものは Phase 5 まで正式機能にしない。

### Phase 3 — 鉛直 1D と hybrid sigma-pressure

実装するもの:

- `A/B` 係数、full/half level、層質量 `Delta p/g`
- 静水圧積分、温位/温度/Exner 関数変換
- 鉛直移流、鉛直 mass flux、必要なら conservative remap
- 上端・下端境界条件と鉛直診断

検証するもの:

- pure pressure (`B=0`) と pure sigma に近い極限がそれぞれ既知の式へ一致すること
- 任意の妥当な `ps` で圧力境界が単調、層質量が正、層質量和が `(ps-p_top)/g` と一致すること
- 等温・乾燥断熱の静水圧柱で geopotential が解析解に一致し収束すること
- 静止柱が加速せず、鉛直 remap が質量・定数トレーサを保存すること
- `ps` 変化時にも定数場保持、positivity、上下境界 flux が整合すること

完了ゲート: 鉛直演算子だけの unit/convergence テストが通り、水平力学から独立して質量・熱力学 budget を閉じられる。

### Phase 4 — 地形なし乾燥 3D 静水圧力学コア

実装するもの:

- cubed-sphere と hybrid 座標を直積した 3D field
- 乾燥 hydrostatic primitive equations、熱力学式、3D tracer 輸送
- surface pressure、鉛直速度、hydrostatic geopotential の診断
- 水平・鉛直 flux の整合、外部重力波に対する時間積分

検証するもの:

- 静止等温大気と静止成層大気を維持すること
- 3D solid-body tracer transport、DCMIP の 3D deformational flow / Hadley-like circulation
- Jablonowski–Williamson の地形なし定常基本場を保ち、摂動時の baroclinic wave を参照解と比較
- 3D の乾燥質量、トレーサ質量、全エネルギー、軸角運動量 budget
- 水平・鉛直解像度を独立に上げた収束、鉛直モード・層間 decoupling の検査
- 長時間積分で上端反射、負圧、温度逸走、面境界ノイズがないこと

完了ゲート: 地形なしの標準 3D テストが複数解像度で再現し、shallow-water 段階の水平誤差と新しい鉛直誤差を区別できる。

### Phase 5 — 地形と下部境界

実装するもの:

- geopotential/orography field と地表境界条件
- 傾いた hybrid 面上の水平圧力傾度項
- 地形データ補間、平滑化、解像度別フィルタ
- 山岳トルクとエネルギー交換の診断

検証するもの:

- 平坦地形が Phase 4 と数値的に同一であること
- DCMIP 2-0-x 型の、山岳上で静止する静水圧大気に偽風が発生しないこと
- 山の高さ・幅・鉛直解像度に対する pressure-gradient error の収束
- Williamson test 5 の山岳付き shallow-water を回帰試験として再利用
- 低振幅の線形山岳波について波長・位相・運動量 flux を理論または高解像度解と比較
- 地形が panel seam/corner を横切る位置へ回転しても結果が同等であること
- Jablonowski–Williamson の smooth topography variant を比較すること

完了ゲート: 静止大気の偽風と pressure-gradient error が解像度とともに減少し、山岳位置を変えても保存量と主要診断が同等になる。

### Phase 6 — 理想化乾燥物理と気候統計

実装するもの:

- Held–Suarez 型 Newtonian relaxation と Rayleigh friction
- physics–dynamics coupling と加熱・摩擦 budget
- ensemble/時間平均、帯状平均、spin-up 除外期間の解析

検証するもの:

- Held–Suarez の帯状平均温度・東西風、Hadley cell、eddy statistics を公開結果と比較
- 時間刻み、解像度、拡散、coupling interval に対する気候統計の感度
- 強制・散逸を含むエネルギー budget が閉じること
- restart が瞬時値だけでなく長時間統計を変えないこと
- 少なくとも複数 seed/微小摂動で内部変動の幅を見積もること

完了ゲート: 単一 snapshot ではなく平均値とばらつきを伴うベンチマーク報告を自動生成できる。

### Phase 7 — 架空惑星化

実装するもの:

- 惑星半径 `R`、重力 `g`、自転 `Omega`、気体定数 `Rd`、`cp`、基準気圧、熱的強制、軌道/恒星日の設定化
- Earth-like、slow rotator、rapid rotator、tidally locked の実験 preset
- 次元量と無次元量（Rossby 数、Burger 数、変形半径、放射/回転時間尺度など）の metadata

検証するもの:

- Earth-like preset が Phase 6 の基準結果を変えないこと
- 単位だけを整合的に変換した相似問題が同じ無次元解を与えること
- `Omega -> 0`、弱い強制、半径変更などの極限で挙動が物理的期待と一致すること
- パラメータ sweep で CFL、層厚、温度、風速が安全範囲を外れたとき明示的に停止すること
- tidally locked terrestrial planet は Heng et al. の乾燥ベンチマーク、将来の物理追加後は THAI protocol と比較すること

完了ゲート: コード変更なしに惑星設定を切替でき、各 preset が設定・診断・再現手順を伴う。

### Phase 8 — 物理拡張と性能（基本コア完成後）

候補:

- 灰色放射、境界層、乾燥対流調節
- 水蒸気・凝結・簡易雲、その後に詳細放射
- 非静水圧方程式（必要な惑星・解像度が明確になった場合）
- OpenMP/MPI、並列 halo exchange、並列 I/O、accelerator 対応

各物理過程は単一 column test と process budget を先に作る。並列化は serial 結果との許容誤差内一致、領域分割不変性、strong/weak scaling を検証し、科学的回帰試験を通過させる。

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

Phase 0 は完了済みであり、次の実装スプリントは Phase 1 だけを対象にする。

1. cubed-sphere の面番号、向き、equiangular gnomonic 写像を ADR で固定する。
2. 座標変換、panel topology、セル面積・辺長・接法線を実装して幾何恒等式を検証する。
3. scalar halo exchange と共有辺 flux を実装し、全球 flux cancellation を検証する。
4. prescribed wind による一次風上輸送と CFL 制御を実装する。
5. MUSCL 再構築と limiter を追加し、滑らかな tracer の二次収束と極値保存を検証する。
6. Williamson test 1 と Lauritzen et al. の標準ケースを CI/回帰テスト化し、Phase 1 report を作る。

## 6. 調査資料と計画への反映

- C. Ronchi, R. Iacono, P. S. Paolucci (1996), [The “Cubed Sphere”: A New Method for the Solution of Partial Differential Equations in Spherical Geometry](https://doi.org/10.1006/jcph.1996.0047): 6 面分割、great-circle coordinate、面結合という cubed-sphere の基礎。
- D. L. Williamson et al. (1992), [A Standard Test Set for Numerical Approximations to the Shallow Water Equations in Spherical Geometry](https://doi.org/10.1016/S0021-9991(05)80016-6): 全球 shallow-water の標準テスト群。
- P. H. Lauritzen et al. (2012), [A standard test case suite for two-dimensional linear transport on the sphere](https://doi.org/10.5194/gmd-5-887-2012): 収束、filament preservation、mixing、非発散/発散輸送の診断。
- J. Galewsky, R. K. Scott, L. M. Polvani (2004), [An initial-value problem for testing numerical models of the global shallow-water equations](https://doi.org/10.3402/tellusa.v56i5.14436): 高速重力波と遅い不安定渦発達を同時に試す jet case。
- J. Thuburn, C. J. Cotter, T. Dubos (2014), [A mimetic, semi-implicit, forward-in-time, finite volume shallow water model](https://doi.org/10.5194/gmd-7-909-2014): cubed-sphere 上の質量・エネルギー・PV 整合性と mimetic 演算子の評価観点。
- A. J. Simmons, D. M. Burridge (1981), [An Energy and Angular-Momentum Conserving Vertical Finite-Difference Scheme and Hybrid Vertical Coordinates](https://doi.org/10.1175/1520-0493(1981)109%3C0758:AEAAMC%3E2.0.CO;2): hybrid 座標と鉛直差分の保存性。
- ECMWF, [Vertical coordinate](https://confluence.ecmwf.int/plugins/viewsource/viewpagesrc.action?pageId=85406187): 実運用される half-level の `p = A + B ps` 定義と上下層の性質。
- C. Jablonowski, D. L. Williamson (2006), [A baroclinic instability test case for atmospheric model dynamical cores](https://doi.org/10.1256/qj.06.12): 地形なし 3D 定常場と baroclinic wave の二段階評価。
- DCMIP (2012), [Dynamical Core Model Intercomparison Project test case document](https://public.websites.umich.edu/~cjablono/DCMIP-2012_TestCaseDocument_v1.7.pdf): 3D 輸送、静水圧静止大気、地形性 pressure-gradient error などの段階的テスト。
- I. M. Held, M. J. Suarez (1994), [A Proposal for the Intercomparison of the Dynamical Cores of Atmospheric General Circulation Models](https://www.gfdl.noaa.gov/bibliography/related_files/ih9401.pdf): 乾燥力学コアの長時間統計ベンチマーク。
- K. Heng, K. Menou, P. J. Phillipps (2011), [Atmospheric circulation of tidally locked exoplanets: a suite of benchmark tests for dynamical solvers](https://doi.org/10.1111/j.1365-2966.2011.18315.x): 架空惑星・同期回転惑星への拡張時の比較ケース。
- E. M. Sergeev et al. (2023), [Simulations of idealised 3D atmospheric flows on terrestrial planets using LFRic-Atmosphere](https://doi.org/10.5194/gmd-16-5601-2023): cubed-sphere を用いた地球型惑星 3D 力学と惑星ベンチマークの現代的な実例。

この計画の閾値は初期値である。各 Phase の最初に、選んだ離散化次数と標準ケースの公開参照値に基づいて数値許容値を確定し、途中で緩和する場合は理由を ADR と検証記録に残す。

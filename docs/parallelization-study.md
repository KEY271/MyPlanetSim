# 並列化調査

調査日: 2026-09-06  
対象commit: `d5ced1f95fd6`

## 結論

- 複数runはコードを変更せず、独立プロセスとして並列実行する。
- 単一runの短縮が必要な場合だけ、dry coreの面状態再構築を鉛直レベル単位で並列化する。
- flux集約、鉛直coupling、診断、forcing、I/Oは現時点では並列化しない。
- production matrixではプロセス並列とスレッド並列を重ねない。

## 測定条件

- `configs/phase7_held_suarez.cfg`（864 cells x 20 levels）
- Release + interprocedural optimization
- Apple M1、4 performance + 4 efficiency cores、16 GB
- 現行直列版: 約22.32 ms/step（300 steps、中央値）

5秒間のsamplingでは、面状態再構築、scalar/vector勾配、limiterが合計約60.7%を占めた。主な残りは鉛直輸送7.2%、`pow` 5.7%、Rusanov flux 4.1%、source 3.1%、Held--Suarez forcing 1.9%だった。

## 複数runのプロセス並列

各プロセスを200 steps実行した。

| processes | wall time | aggregate speedup |
|---:|---:|---:|
| 1 | 4.47 s | 1.00x |
| 2 | 4.66 s | 1.92x |
| 4 | 5.17 s | 3.46x |
| 6 | 7.70 s | 3.48x |
| 8 | 9.94 s | 3.60x |

このホストでは4 processes付近が効率上の飽和点である。6本を一括実行すればwall timeは約3.5倍短縮できるが、4本を超えるとaggregate throughputはほぼ伸びない。

## 単一runの試作

作業ツリー外の一時コピーで、面状態再構築のlevel loopだけをOpenMP並列化した。各threadは永続scratch workspaceを持ち、出力は異なる`(edge, level)`へ書き込む。

| threads | median ms/step | whole-run speedup |
|---:|---:|---:|
| 1 | 22.32 | 1.00x |
| 2 | 16.42 | 1.36x |
| 4 | 13.13 | 1.70x |
| 8 | 13.48 | 1.66x |

4 threadsが最良だった。300 stepsのCPU timeは1 threadで6.66 s、4 threadsで7.07 sであり、CPU time約6%増に対してwall timeは約41%減少した。8 threadsは4 threadsより遅く、CPU timeも増える。

試作はCMake、reconstruction header、reconstruction実装の3ファイルに限定できた。warm-up後のRHS allocationは0を維持し、以下のfocused testsは4 threadsで通過した。

- `unit.dry_hydrostatic_driver`
- `unit.rhs_allocation`
- `phase7.forced_dry_core`
- `phase5.dry_hydrostatic_baseline`
- `phase9.balanced_benchmarks`

## 実装する場合の境界

`reconstruct_dry_hydrostatic_face_states`の外側のlevel loopだけを`parallel for`にする。threadごとのworkspaceと、整数の`limiter_activations` reductionを使い、内側のcell/edge loopは並列化しない。OpenMPは任意のCMake optionとし、1-thread gateも残す。

次は効果または複雑度の理由から対象外とする。

- edge flux集約: 同じcellへの更新が競合し、約4%のためatomic等に見合わない。
- vertical coupling: thread別scratchが必要で、追加効果は小さい。
- diagnosis/forcing: 比率が小さく、浮動小数点reductionの再現性が複雑になる。
- SSP-RK3 stages: stage間に依存関係がある。
- I/O/checkpoint: wall timeの5%を超える測定根拠がない。

6 processesを同時実行した比較では、1 thread/processの7.00 sに対し4 threads/processは6.56 sだった。約7%のwall短縮に対してCPU timeが約15%増えたため、複数runでは1 thread/processを使う。

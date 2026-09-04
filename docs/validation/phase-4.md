# Phase 4 検証報告 — hybrid 鉛直 column

Phase 4 は水平格子と結合しない C++ column core として実装した。`vertical_column` は
`A/B` hybrid interface、log-pressure full level、dry Exner/温位変換、静水圧
geopotential、mass-flux recurrence、保存型 scalar transport を所有する。Web、gateway、
既存 shallow-water frame は変更していない。

## 再現手順

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev -R 'phase4|unit.experiment_config' --output-on-failure
./build/dev/my_planet_sim --config configs/phase4_isothermal.cfg
./build/dev/my_planet_sim --config configs/phase4_dry_adiabatic.cfg
./build/dev/my_planet_sim --config configs/phase4_moving_surface_pressure.cfg
./build/dev/my_planet_sim --config configs/phase4_manufactured_transport.cfg
```

`phase4.vertical_column`、`phase4.csv_output`、`phase4.checkpoint_restart` は以下を確認する。

- hybrid mass の telescoping と full pressure の layer 内部性
- Exner/温位の round trip と layer-wise hydrostatic residual
- top/bottom impermeable mass-flux recurrence と donor-cell constant preservation
- moving-`ps` の restart 一致
- 静止 column の予報 state が bitwise に不変であること
- donor-cell の一次、unlimited linear の二次空間収束
- SSP-RK3 の三次以上の時間収束
- minmod の保存性、boundedness、nonnegative tracer
- 非ゼロ `run.start_time_s` と pressure-bound stage の自動細分化
- source-integrated dry/theta/tracer budget、上下端 flux、continuity、実 CFL
- version 2 checkpoint と profile/diagnostics CSV schema
- 実ファイル checkpoint からの CLI restart と連続実行結果の bitwise 一致

`nz=32/64` の smooth manufactured flux に対する L2 error と observed order は次のとおり。

| Scheme | `nz=32` | `nz=64` | Observed order |
|---|---:|---:|---:|
| donor-cell | `8.03573e-2` | `4.04676e-2` | `0.98966` |
| unlimited linear | `3.85135e-3` | `9.56647e-4` | `2.00931` |

quarter-period surface-pressure solutionを `dt=2,1,0.5 s` で比較した error は
`3.42718e-2, 2.11244e-3, 1.31573e-4 Pa`、observed order は
`4.02004, 4.00498` だった。time-only forcing では SSP-RK3 の stage quadrature が
この問題に対して四次精度となり、三次 gate を満たす。

標準 moving-surface-pressure run では最大実 CFL は `1.2626262626262648e-4`、
continuity residual と上下端 flux は `0` だった。最終 source-integrated residual は
dry mass `-3.09e-11 kg m-2`、theta mass `-5.12e-9 K kg m-2`、tracer mass
`-2.50e-12 kg m-2` である。dev の全 CTest、GCC 15 と ASan/UBSan の Phase 4 test、
format-check も成功した。

各 CLI run は `column_profile.csv`、`column_diagnostics.csv`、必要に応じて
`vertical_column_hybrid_v1` checkpoint を出力する。Phase 5 で水平との結合、3D state、
browser-facing frame は別途検証する。

# Phase 4 検証報告 — hybrid 鉛直 column

Phase 4 は水平格子と結合しない C++ column core として実装した。`vertical_column` は
`A/B` hybrid interface、log-pressure full level、dry Exner/温位変換、静水圧
geopotential、mass-flux recurrence、保存型 scalar transport を所有する。Web、gateway、
既存 shallow-water frame は変更していない。

## 再現手順

```sh
cmake --preset dev
cmake --build --preset dev
ctest --preset dev -R 'phase4.vertical_column|smoke.phase4_vertical_column|unit.experiment_config' --output-on-failure
./build/dev/my_planet_sim --config configs/phase4_isothermal.cfg
./build/dev/my_planet_sim --config configs/phase4_moving_surface_pressure.cfg
```

`phase4.vertical_column` は以下を確認する。

- hybrid mass の telescoping と full pressure の layer 内部性
- Exner/温位の round trip と layer-wise hydrostatic residual
- top/bottom impermeable mass-flux recurrence と donor-cell constant preservation
- moving-`ps` の restart 一致
- 静止 column の予報 state が bitwise に不変であること

各 CLI run は `column_profile.csv`、`column_diagnostics.csv`、必要に応じて
`vertical_column_hybrid_v1` checkpoint を出力する。Phase 5 で水平との結合、3D state、
browser-facing frame は別途検証する。

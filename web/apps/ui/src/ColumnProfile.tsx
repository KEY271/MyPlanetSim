import { FrameV2Field, FrameV2FieldDescriptor, selectedColumn, VisualFrameV2 } from "@myplanetsim/protocol";
import { formatValue } from "./format";

interface ColumnProfileProps {
  readonly frame: VisualFrameV2;
  readonly field: FrameV2FieldDescriptor;
  readonly cell: number;
  readonly level: number;
  readonly onSelectLevel: (level: number) => void;
}

const width = 300;
const height = 320;
const plot = { left: 62, right: 286, top: 20, bottom: 272 };

// Pressure is plotted on a downward-increasing logarithmic axis, which is the standard
// reading order for an atmospheric column: the model top is at the top of the chart.
function pressureToY(pressurePa: number, topPa: number, bottomPa: number): number {
  const span = Math.log(bottomPa) - Math.log(topPa);
  const fraction = span > 0 ? (Math.log(pressurePa) - Math.log(topPa)) / span : 0.5;
  return plot.top + fraction * (plot.bottom - plot.top);
}

export function ColumnProfile({ frame, field, cell, level, onSelectLevel }: ColumnProfileProps) {
  const column = selectedColumn(frame, field.id as FrameV2Field, cell);
  const minimum = column.values.reduce((low, value) => Math.min(low, value), Infinity);
  const maximum = column.values.reduce((high, value) => Math.max(high, value), -Infinity);
  // A constant column has no horizontal extent, so it is centred rather than divided by
  // zero; the axis labels still report the single value.
  const span = maximum - minimum;
  const valueToX = (value: number) => span > 0
    ? plot.left + ((value - minimum) / span) * (plot.right - plot.left)
    : (plot.left + plot.right) / 2;
  const topPa = column.pressurePa[0];
  const bottomPa = column.pressurePa[column.pressurePa.length - 1];
  const points = Array.from(column.values, (value, index) =>
    `${valueToX(value)},${pressureToY(column.pressurePa[index], topPa, bottomPa)}`);
  const tickStride = Math.max(1, Math.ceil(frame.levels / 6));
  const levelTicks = Array.from(column.pressurePa, (pressurePa, index) => ({ index, pressurePa }))
    .filter((tick) => tick.index % tickStride === 0 || tick.index === frame.levels - 1);
  return <article className="view-panel column-profile">
    <h2>Column profile</h2>
    <div className="column-profile-body">
      <svg viewBox={`0 0 ${width} ${height}`} role="img" preserveAspectRatio="xMidYMid meet"
        aria-label={`${field.label} against pressure for cell ${cell}: ${frame.levels} model levels between ${formatValue(minimum)} and ${formatValue(maximum)} ${field.unit}`}>
        <line x1={plot.left} y1={plot.top} x2={plot.left} y2={plot.bottom} className="axis" />
        <line x1={plot.left} y1={plot.bottom} x2={plot.right} y2={plot.bottom} className="axis" />
        {levelTicks.map((tick) => <g key={tick.index}>
          <line x1={plot.left - 4} x2={plot.left} className="axis"
            y1={pressureToY(tick.pressurePa, topPa, bottomPa)}
            y2={pressureToY(tick.pressurePa, topPa, bottomPa)} />
          <text x={plot.left - 7} y={pressureToY(tick.pressurePa, topPa, bottomPa) + 3}
            textAnchor="end" className="tick">{Math.round(tick.pressurePa / 100)}</text>
        </g>)}
        <text x={plot.left} y={plot.bottom + 16} textAnchor="start" className="tick">{formatValue(minimum)}</text>
        <text x={plot.right} y={plot.bottom + 16} textAnchor="end" className="tick">{formatValue(maximum)}</text>
        <text x={(plot.left + plot.right) / 2} y={height - 4} textAnchor="middle" className="tick">
          {field.label} ({field.unit})
        </text>
        <text x={12} y={(plot.top + plot.bottom) / 2} textAnchor="middle" className="tick"
          transform={`rotate(-90 12 ${(plot.top + plot.bottom) / 2})`}>Pressure (hPa)</text>
        <polyline className="profile-line" fill="none" points={points.join(" ")} />
        {points.map((point, index) => {
          const [x, y] = point.split(",");
          return <circle key={index} cx={x} cy={y} r={index === level ? 5 : 2.6}
            className={index === level ? "profile-point selected" : "profile-point"}
            onClick={() => onSelectLevel(index)}>
            <title>{`Level ${index}: ${formatValue(column.values[index])} ${field.unit} at ${formatValue(column.pressurePa[index])} Pa`}</title>
          </circle>;
        })}
      </svg>
      <p className="profile-status">
        Cell {cell} · level {level} · {formatValue(column.values[level])} {field.unit} at{" "}
        {formatValue(column.pressurePa[level])} Pa
      </p>
    </div>
  </article>;
}

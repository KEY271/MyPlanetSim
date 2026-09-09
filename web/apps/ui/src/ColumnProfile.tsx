import { formatValue } from "./format";
import { profileGeometry, profilePlot, profileSize } from "./profile";

interface ColumnProfileProps {
  readonly values: Float64Array;
  readonly pressurePa: Float64Array;
  readonly label: string;
  readonly unit: string;
  readonly cell: number;
  readonly level: number;
  readonly onSelectLevel: (level: number) => void;
}
export function ColumnProfile({ values, pressurePa, label, unit, cell, level,
  onSelectLevel }: ColumnProfileProps) {
  const { points, ticks, minimum, maximum } = profileGeometry({ values, pressurePa });
  const centreX = (profilePlot.left + profilePlot.right) / 2;
  const centreY = (profilePlot.top + profilePlot.bottom) / 2;
  return <article className="view-panel column-profile">
    <h2>Column profile</h2>
    <div className="column-profile-body">
      <svg viewBox={`0 0 ${profileSize.width} ${profileSize.height}`} role="img" preserveAspectRatio="xMidYMid meet"
        aria-label={`${label} against pressure for cell ${cell}: ${values.length} model levels between ${formatValue(minimum)} and ${formatValue(maximum)} ${unit}`}>
        <line className="axis" x1={profilePlot.left} y1={profilePlot.top} x2={profilePlot.left} y2={profilePlot.bottom} />
        <line className="axis" x1={profilePlot.left} y1={profilePlot.bottom} x2={profilePlot.right} y2={profilePlot.bottom} />
        {ticks.map((tick) => <g key={tick.index}><line className="axis" x1={profilePlot.left - 4} x2={profilePlot.left} y1={tick.y} y2={tick.y} />
          <text className="tick" x={profilePlot.left - 7} y={tick.y + 3} textAnchor="end">{Math.round(tick.pressurePa / 100)}</text></g>)}
        <text className="tick" x={profilePlot.left} y={profilePlot.bottom + 16} textAnchor="start">{formatValue(minimum)}</text>
        <text className="tick" x={profilePlot.right} y={profilePlot.bottom + 16} textAnchor="end">{formatValue(maximum)}</text>
        <text className="tick" x={centreX} y={profileSize.height - 4} textAnchor="middle">{label} ({unit})</text>
        <text className="tick" x={12} y={centreY} textAnchor="middle" transform={`rotate(-90 12 ${centreY})`}>Pressure (hPa)</text>
        <polyline className="profile-line" fill="none" points={points.map((point) => `${point.x},${point.y}`).join(" ")} />
        {points.map((point) => <circle key={point.index} cx={point.x} cy={point.y} r={point.index === level ? 5 : 2.6}
          className={point.index === level ? "profile-point selected" : "profile-point"} onClick={() => onSelectLevel(point.index)}>
          <title>{`Level ${point.index}: ${formatValue(values[point.index])} ${unit} at ${formatValue(pressurePa[point.index])} Pa`}</title>
        </circle>)}
      </svg>
      <p className="profile-status">Cell {cell} · level {level} · {formatValue(values[level])} {unit} at {formatValue(pressurePa[level])} Pa</p>
    </div>
  </article>;
}

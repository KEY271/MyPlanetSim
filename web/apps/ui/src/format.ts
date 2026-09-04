// Frame values are raw C++ doubles that span surface pressure (1e5 Pa) and tracer
// mixing ratios (1e-3), so the inspector and the profile share one significant-figure
// formatter instead of rounding each field to a hand-picked precision.
export function formatValue(value: number, significantDigits = 5): string {
  if (!Number.isFinite(value)) return String(value);
  if (value === 0) return "0";
  const magnitude = Math.abs(value);
  if (magnitude >= 1e6 || magnitude < 1e-3) return value.toExponential(significantDigits - 1);
  return String(Number(value.toPrecision(significantDigits)));
}

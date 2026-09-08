"""Exercise comparison-fixture import against a real simulator checkpoint."""
import csv
import math
from pathlib import Path
import subprocess
import sys
import tempfile

simulator, benchmark, preset = sys.argv[1:]


def run(command):
    return subprocess.run(command, text=True, capture_output=True, check=True)


with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    config = root / 'source.cfg'
    text = Path(preset).read_text().replace('grid.cells_per_panel = 6',
                                          'grid.cells_per_panel = 2')
    text = text.replace('output.directory = output/phase13_moist_pilot',
                        f'output.directory = {root / "output"}')
    config.write_text(text)
    checkpoint = root / 'source.chk'
    run([simulator, '--config', str(config), '--stop-after-step', '1',
         '--checkpoint', str(checkpoint), '--progress-interval-s', '0'])
    snapshot = root / 'snapshot.csv'
    command = [benchmark, str(config), str(1800 / 86400), '2', '20', '1800',
               '2', '20000', '5', str(snapshot), str(checkpoint), str(config)]
    output = run(command).stdout
    metrics = dict(line.split('=', 1) for line in output.splitlines())
    assert float(metrics['initial_time_s']) == 1800
    assert float(metrics['accepted_seconds']) == 1800
    assert float(metrics['full_rhs_per_step']) == 4
    rows = list(csv.DictReader(Path(str(snapshot) + '.samples.csv').open()))
    assert len(rows) == 2
    area = 4 * math.pi * 6371220**2
    evaporation = (float(rows[-1]['cumulative_evaporation_kg']) -
                   float(rows[0]['cumulative_evaporation_kg'])) / (area * 1800)
    assert math.isclose(float(metrics['evaporation_kg_m2_s']), evaporation,
                        rel_tol=1e-12)
    first = snapshot.read_bytes()
    run(command)
    assert snapshot.read_bytes() == first
    candidate = root / 'changed.cfg'
    candidate.write_text(text.replace('surface.albedo = 0.3', 'surface.albedo = 0.4'))
    command[1] = str(candidate)
    result = subprocess.run(command, text=True, capture_output=True)
    assert result.returncode != 0 and 'physical configuration differs' in result.stderr
    command[1] = str(config)
    command[8] = '4'
    result = subprocess.run(command, text=True, capture_output=True)
    assert result.returncode != 0 and 'top/refinement differ' in result.stderr

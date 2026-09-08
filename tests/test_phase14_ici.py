"""Check mass weighting and reject incompatible/nonfinite comparison snapshots."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location(
    'ici', Path(__file__).resolve().parents[1] / 'tools/compare_phase14_ici.py')
ici = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ici)


class SnapshotComparison(unittest.TestCase):
    def test_reference_mass_weighted_rms(self):
        with tempfile.TemporaryDirectory() as directory:
            candidate, reference = (Path(directory) / name for name in ('a.csv', 'b.csv'))
            header = 'cell,level,area_m2,air_mass_kg_m2,temperature_k\n'
            reference.write_text(header + '0,0,2,1,280\n0,1,2,3,290\n')
            candidate.write_text(header + '0,0,2,7,282\n0,1,2,8,290\n')
            self.assertEqual(ici.temperature_rms(candidate, reference), 1)
            self.assertEqual(ici.temperature_rms(reference, reference), 0)
            for rows in ('0,0,2,1,nan\n0,1,2,3,290\n',
                         '0,1,2,1,280\n0,0,2,3,290\n',
                         '0,0,2,1,280\n'):
                candidate.write_text(header + rows)
                with self.assertRaises(ValueError):
                    ici.temperature_rms(candidate, reference)


if __name__ == '__main__':
    unittest.main()

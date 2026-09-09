"""Check the Phase 14 final comparison metrics and their rejections."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location(
    'final', Path(__file__).resolve().parents[1] / 'tools/compare_phase14_final.py')
final = importlib.util.module_from_spec(spec)
spec.loader.exec_module(final)

HEADER = ('cell,level,area_m2,air_mass_kg_m2,temperature_k,zonal_wind_m_s,'
          'latitude_deg\n')
REFERENCE = ('0,0,2,1,280,4,-45\n'
             '0,1,2,3,290,6,-45\n'
             '1,0,2,1,300,10,45\n'
             '1,1,2,3,310,12,45\n')


def rows(text):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / 'snapshot.csv'
        path.write_text(HEADER + text)
        return final.snapshot_rows(path)


class Metrics(unittest.TestCase):
    def test_mass_weighted_rms_uses_reference_mass(self):
        reference = rows(REFERENCE)
        candidate = rows('0,0,2,7,282,4,-45\n0,1,2,8,290,6,-45\n'
                         '1,0,2,1,300,10,45\n1,1,2,3,310,12,45\n')
        # Only the first row differs, by 2 K, and carries mass 2*1 of the total 2*8.
        self.assertAlmostEqual(final.mass_weighted_temperature_rms(candidate, reference),
                               (4 * 2 / 16) ** 0.5)
        self.assertEqual(final.mass_weighted_temperature_rms(reference, reference), 0)

    def test_global_mean_temperature_is_mass_weighted(self):
        reference = rows(REFERENCE)
        self.assertAlmostEqual(final.global_mean_temperature(reference),
                               (280 + 3 * 290 + 300 + 3 * 310) / 8)

    def test_zonal_rms_compares_bands_of_equal_level(self):
        reference = rows(REFERENCE)
        candidate = rows('0,0,2,1,281,4,-45\n0,1,2,3,290,6,-45\n'
                         '1,0,2,1,300,13,45\n1,1,2,3,310,12,45\n')
        # Four (band, level) pairs; one differs by 1 K.
        self.assertAlmostEqual(final.zonal_rms(candidate, reference, 10, 'temperature_k'),
                               (1 / 4) ** 0.5)
        self.assertAlmostEqual(final.zonal_rms(candidate, reference, 10,
                                               'zonal_wind_m_s'), (9 / 4) ** 0.5)
        self.assertEqual(final.zonal_rms(reference, reference, 10, 'zonal_wind_m_s'), 0)

    def test_incompatible_snapshots_are_rejected(self):
        reference = rows(REFERENCE)
        for text in ('0,0,2,1,nan,4,-45\n0,1,2,3,290,6,-45\n'
                     '1,0,2,1,300,10,45\n1,1,2,3,310,12,45\n',
                     '0,1,2,1,280,4,-45\n0,0,2,3,290,6,-45\n'
                     '1,0,2,1,300,10,45\n1,1,2,3,310,12,45\n',
                     '0,0,2,1,280,4,-45\n'):
            candidate = rows(text)
            with self.assertRaises(ValueError):
                final.mass_weighted_temperature_rms(candidate, reference)

    def test_block_statistics_use_offsets_and_ledger_differences(self):
        header = ('day,mean_temperature_k,mean_surface_temperature_k,'
                  'precipitable_water_kg_m2,cumulative_evaporation_kg,'
                  'cumulative_convective_precipitation_kg,'
                  'cumulative_grid_scale_precipitation_kg,cumulative_runoff_kg\n')
        # A run that starts at day 1200 uses the same block as one starting at zero.
        body = ('1200,280,290,60,0,0,0,0\n'
                '1201,281,291,62,86400,43200,43200,0\n'
                '1202,283,293,66,259200,129600,129600,0\n')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'run.csv.samples.csv'
            path.write_text(header + body)
            block = final.block_statistics(path, 1, 2)
            self.assertEqual(block['block_days'], 1)
            self.assertEqual(block['mean_mean_temperature_k'], 283)
            self.assertEqual(block['mean_precipitable_water_kg_m2'], 66)
            # (259200 - 86400) kg over one day.
            self.assertAlmostEqual(block['evaporation_kg_s'], 2.0)
            self.assertAlmostEqual(block['precipitation_kg_s'], 2.0)
            self.assertAlmostEqual(block['convective_fraction'], 0.5)
            with self.assertRaises(ValueError):
                final.block_statistics(path, 2, 2)

    def test_missing_benchmark_output_is_rejected(self):
        with self.assertRaises(ValueError):
            final.metrics('accepted_steps=48\n')


if __name__ == '__main__':
    unittest.main()

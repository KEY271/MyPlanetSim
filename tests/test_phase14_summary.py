"""Regression checks for right-closed integrated diagnostic intervals."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location(
    'summary', Path(__file__).resolve().parents[1] / 'tools/summarize_moist_pilot.py')
summary = importlib.util.module_from_spec(spec)
spec.loader.exec_module(summary)


def row(day):
    return dict(time_s=str(day * 86400), interval_seconds='86400',
                net_evaporation_kg='1', convective_precipitation_kg='1',
                grid_scale_precipitation_kg='0', precipitable_water_kg_m2='10')


class Intervals(unittest.TestCase):
    def test_spinup_and_blocks(self):
        selected = summary.select_intervals([row(d) for d in range(1, 1201)],
                                            200 * 86400, 1200 * 86400)
        self.assertEqual(len(selected), 1000)
        blocks = summary.block_statistics(selected, 1, 0, 100)
        self.assertEqual(len(blocks), 10)
        self.assertEqual(blocks[0]['start_day'], 200)
        self.assertEqual(blocks[-1]['end_day'], 1200)
        self.assertTrue(all(b['seconds'] == 100 * 86400 for b in blocks))

    def test_straddling_interval_is_not_invented(self):
        with self.assertRaises(ValueError):
            summary.select_intervals([row(201)], 200.5 * 86400)


if __name__ == '__main__':
    unittest.main()

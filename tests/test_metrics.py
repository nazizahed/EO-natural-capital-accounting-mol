import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nca_prototype import calculate_extent_account, exploratory_condition_index, safe_normalized_difference


class MetricTests(unittest.TestCase):
    def test_normalized_difference(self):
        a = np.array([3.0, 0.0, np.nan])
        b = np.array([1.0, 0.0, 1.0])
        result = safe_normalized_difference(a, b)
        self.assertAlmostEqual(float(result[0]), 0.5)
        self.assertTrue(np.isnan(result[1]))
        self.assertTrue(np.isnan(result[2]))

    def test_extent_account_area_and_share(self):
        land_cover = np.array([[10, 10], [40, 50]], dtype="uint8")
        table = calculate_extent_account(land_cover, pixel_area_ha=0.01)
        tree = table.loc[table["worldcover_code"] == 10].iloc[0]
        self.assertAlmostEqual(float(tree["extent_ha"]), 0.02)
        self.assertAlmostEqual(float(tree["share_of_mapped_area_pct"]), 50.0)

    def test_condition_index_bounds(self):
        ndvi = np.array([-1.0, 0.4, 1.0, np.nan])
        ndmi = np.array([-1.0, 0.2, 1.0, 0.0])
        score = exploratory_condition_index(ndvi, ndmi)
        self.assertTrue(np.nanmin(score) >= 0)
        self.assertTrue(np.nanmax(score) <= 100)
        self.assertTrue(np.isnan(score[-1]))


if __name__ == "__main__":
    unittest.main()

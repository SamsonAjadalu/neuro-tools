#!/usr/bin/env python3
import unittest
from pathlib import Path
import sys

import numpy as np


sys.path.insert(0, str(Path(__file__).parent))
from make_nifti_mosaic import automatic_axial_views


class AutomaticAxialViewTests(unittest.TestCase):
    def test_selected_slices_are_eligible_with_and_without_mask(self):
        overlay = np.zeros((2, 2, 10), dtype=float)
        overlay[:, :, [1, 4, 7, 9]] = 1
        affine = np.eye(4)

        for mask in (None, np.where(overlay > 0, 1, 0)):
            views = automatic_axial_views(
                overlay,
                mask,
                affine,
                n_slices=3,
                min_active_percent=50,
                min_active_voxels=0,
            )
            selected = [view.voxel_index for view in views]
            self.assertEqual(selected, [1, 7, 9])
            self.assertTrue(set(selected).issubset({1, 4, 7, 9}))

    def test_fewer_eligible_slices_are_not_duplicated(self):
        overlay = np.zeros((2, 2, 5), dtype=float)
        overlay[:, :, 2] = 1

        views = automatic_axial_views(
            overlay,
            None,
            np.eye(4),
            n_slices=4,
            min_active_percent=50,
            min_active_voxels=0,
        )

        self.assertEqual([view.voxel_index for view in views], [2])


if __name__ == "__main__":
    unittest.main()

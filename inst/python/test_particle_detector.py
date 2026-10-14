"""Tests for particle_detector.py on synthetic LDIR-like images with known answers.

Run from the repository root:
    python -m unittest discover -s inst/python -p "test_*.py"
(also run by tests/testthat/test-python-detector.R through reticulate).

The synthetic image mimics an LDIR particle map: a near-black frame with
tile shading, and bright round particles at known positions.

Known limitation (not tested as correct): detect_scan_circle() thresholds at
the 80th brightness percentile, so it can only isolate a scan disk covering
less than ~20% of the frame. On examples/dataset_A/LDIR_image.png (a square
crop with no visible circle) it returns a wrong circle; the R side rejects it
(sanity check) and falls back. The pipeline uses the analyzed overlay there.
"""
import os
import sys
import tempfile
import unittest

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import particle_detector as pd  # noqa: E402

H, W = 600, 640
# Scan area used for the circle-mask test: (cx, cy, r) px
DISK = (330.0, 300.0, 260.0)
# Well separated particles, all inside DISK: (cx, cy, radius px)
BLOBS = [(200, 150, 6), (330, 120, 5), (450, 170, 7), (150, 300, 4),
         (260, 260, 8), (400, 290, 5), (520, 320, 6), (210, 420, 5),
         (330, 400, 7), (450, 440, 4), (300, 520, 6), (380, 210, 5)]
OUTSIDE = (30, 30, 6)                 # bright blob outside the scan disk


def make_image(path, blobs=BLOBS, extra=()):
    """Near-black frame with slightly different shading per quarter (the LDIR
    mosaic tiles, as in examples/dataset_A/LDIR_image.png) and bright
    particles."""
    yy, xx = np.mgrid[0:H, 0:W].astype(float)
    img = np.full((H, W), 4.0)
    img[:H // 2, W // 2:] += 4.0
    img[H // 2:, :W // 2] += 2.0
    img += 3.0 * xx / W                                 # gentle gradient
    for bx, by, br in list(blobs) + list(extra):
        img[(xx - bx) ** 2 + (yy - by) ** 2 <= br ** 2] = 230.0
    rgb = np.stack([img * 0.6, img * 0.8, img], axis=-1).clip(0, 255).astype(np.uint8)
    Image.fromarray(rgb).save(path)


def make_disk_image(path, cx, cy, r, h=400, w=400):
    """A bright scan disk on black covering under 20% of the frame -- the case
    detect_scan_circle() is built for (it thresholds at the 80th percentile)."""
    yy, xx = np.mgrid[0:h, 0:w].astype(float)
    img = np.where((xx - cx) ** 2 + (yy - cy) ** 2 <= r ** 2, 120.0, 3.0)
    Image.fromarray(np.stack([img] * 3, -1).astype(np.uint8)).save(path)


def centroids(result):
    p = result['particles']
    return np.column_stack([p['centroid_x'], p['centroid_y']])


def nearest_error(found, truth):
    d = np.sqrt(((found[:, None, :] - truth[None, :, :]) ** 2).sum(-1))
    return d.min(axis=0)               # for each true blob, nearest detection


class DetectorTest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, 'ldir.png')
        make_image(self.path)

    def tearDown(self):
        self.tmp.cleanup()

    def test_scan_circle(self):
        p = os.path.join(self.tmp.name, 'disk.png')
        make_disk_image(p, cx=210.0, cy=190.0, r=90.0)
        c = pd.detect_scan_circle(p)
        self.assertEqual(c['method'], 'python_cc')
        self.assertLess(abs(c['cx'] - 210), 2)
        self.assertLess(abs(c['cy'] - 190), 2)
        self.assertEqual((c['width'], c['height']), (400, 400))

    @unittest.expectedFailure
    def test_scan_circle_radius_unbiased(self):
        # KNOWN ISSUE: the radius comes out ~12% too large here. The disk is
        # smoothed (sigma = 1% of the image) and thresholded just above the
        # background, so the blurred halo counts as disk. On an optical LDIR
        # image the radius sets um/px, so this would scale LDIR sizes and
        # distances. Needs a real optical export with a visible scan circle
        # to confirm before changing the detector (docs/TODO.md).
        p = os.path.join(self.tmp.name, 'disk.png')
        make_disk_image(p, cx=210.0, cy=190.0, r=90.0)
        c = pd.detect_scan_circle(p)
        self.assertLess(abs(c['r'] / 90 - 1), 0.05)

    def test_scan_circle_fallback_without_disk(self):
        p = os.path.join(self.tmp.name, 'flat.png')
        Image.fromarray(np.full((100, 120, 3), 40, np.uint8)).save(p)
        c = pd.detect_scan_circle(p)
        self.assertEqual(c['method'], 'fallback_center')

    def test_detects_every_particle_at_its_position(self):
        r = pd.run_full_pipeline(self.path, merge_fibers=False)
        self.assertEqual(r['n_particles'], len(BLOBS))
        err = nearest_error(centroids(r), np.array([b[:2] for b in BLOBS], float))
        self.assertLess(err.max(), 1.0)          # sub-pixel centroids
        # Area follows the blob radius (pi r^2, within pixelation)
        found = centroids(r)
        areas = np.array(r['particles']['area_px'], float)
        for bx, by, br in BLOBS:
            i = int(np.argmin(((found - [bx, by]) ** 2).sum(1)))
            self.assertLess(abs(areas[i] / (np.pi * br ** 2) - 1), 0.35)

    def test_target_count_tunes_to_the_expected_number(self):
        r = pd.run_full_pipeline(self.path, target_count=len(BLOBS), merge_fibers=False)
        self.assertEqual(r['n_particles'], len(BLOBS))

    def test_circle_mask_excludes_particles_outside_the_scan(self):
        p = os.path.join(self.tmp.name, 'outside.png')
        make_image(p, extra=[OUTSIDE])
        open_ = pd.run_full_pipeline(p, merge_fibers=False)
        self.assertEqual(open_['n_particles'], len(BLOBS) + 1)
        cx, cy, r = DISK
        masked = pd.run_full_pipeline(p, circle_cx=cx, circle_cy=cy, circle_r=r,
                                      merge_fibers=False)
        self.assertEqual(masked['n_particles'], len(BLOBS))
        self.assertGreater(nearest_error(centroids(masked),
                                         np.array([OUTSIDE[:2]], float))[0], 50)

    def test_deterministic(self):
        a = pd.run_full_pipeline(self.path, target_count=len(BLOBS))
        b = pd.run_full_pipeline(self.path, target_count=len(BLOBS))
        np.testing.assert_array_equal(centroids(a), centroids(b))
        self.assertEqual(a['threshold_used'], b['threshold_used'])


if __name__ == '__main__':
    unittest.main()

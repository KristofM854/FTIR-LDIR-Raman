"""
particle_detector.py
Called from R via reticulate. All functions return numpy arrays or dicts
that reticulate can convert to R data structures.
"""

import numpy as np
from scipy import ndimage
from scipy.ndimage import label, gaussian_filter, center_of_mass, find_objects
from PIL import Image


def load_and_prepare(image_path):
    """Load image and convert to grayscale float64 array.

    Args:
        image_path: Path to the LDIR mosaic PNG.

    Returns:
        dict with 'rgb' (H,W,3 uint8), 'gray' (H,W float64),
        'height', 'width' keys.
    """
    img = Image.open(image_path)
    arr = np.array(img, dtype=np.float64)

    # Handle RGBA by dropping alpha
    if arr.ndim == 3 and arr.shape[2] == 4:
        arr = arr[:, :, :3]

    # Handle grayscale
    if arr.ndim == 2:
        gray = arr.copy()
        rgb = np.stack([arr, arr, arr], axis=-1).astype(np.uint8)
    else:
        gray = 0.299 * arr[:, :, 0] + 0.587 * arr[:, :, 1] + 0.114 * arr[:, :, 2]
        rgb = arr.astype(np.uint8)

    return {
        'rgb': rgb,
        'gray': gray,
        'height': int(arr.shape[0]),
        'width': int(arr.shape[1])
    }


def correct_background(gray, grid_rows=4, grid_cols=4, bg_sigma=30.0,
                        clip_sigma=3.0, max_iter=10):
    """Per-tile iterative sigma-clipping background subtraction.

    Args:
        gray: 2D float64 array (grayscale image).
        grid_rows, grid_cols: Tile grid dimensions.
        bg_sigma: Gaussian smoothing sigma for background estimation.
        clip_sigma: Number of standard deviations for clipping threshold.
        max_iter: Maximum sigma-clipping iterations.

    Returns:
        dict with:
          'corrected': 2D float64 array (background-subtracted, clipped >= 0)
          'background': 2D float64 array (estimated background)
          'tile_height': int
          'tile_width': int
    """
    h, w = gray.shape
    tile_h = h // grid_rows
    tile_w = w // grid_cols

    corrected = np.zeros_like(gray)
    background = np.zeros_like(gray)

    for ty in range(grid_rows):
        for tx in range(grid_cols):
            y0, y1 = ty * tile_h, (ty + 1) * tile_h
            x0, x1 = tx * tile_w, (tx + 1) * tile_w
            tile = gray[y0:y1, x0:x1].copy()

            # Scale sigma proportionally to tile size
            sigma_scaled = bg_sigma * (tile_h / 500.0)

            # Iterative sigma-clipping
            mask = np.ones_like(tile, dtype=bool)
            bg = gaussian_filter(tile, sigma=sigma_scaled)

            for iteration in range(max_iter):
                residual = tile - bg
                mu = np.mean(residual[mask])
                sigma_noise = np.std(residual[mask])
                new_mask = residual < (mu + clip_sigma * sigma_noise)

                if np.sum(new_mask) == np.sum(mask):
                    break

                mask = new_mask
                filled = tile.copy()
                filled[~mask] = bg[~mask]
                bg = gaussian_filter(filled, sigma=sigma_scaled)

            background[y0:y1, x0:x1] = bg
            corrected[y0:y1, x0:x1] = tile - bg

    corrected = np.clip(corrected, 0, None)

    return {
        'corrected': corrected,
        'background': background,
        'tile_height': int(tile_h),
        'tile_width': int(tile_w)
    }


def detect_particles(corrected, threshold=25.0, min_area=20):
    """Global thresholding + connected component analysis.

    IMPORTANT: Detection is GLOBAL (whole image), not per-tile,
    because ~3% of particles span tile boundaries.

    Args:
        corrected: 2D float64 background-subtracted image.
        threshold: Intensity threshold for detection.
        min_area: Minimum particle area in pixels.

    Returns:
        dict with:
          'labeled': 2D int array (labeled particles, 0=background)
          'n_particles': int
          'binary': 2D bool array (thresholded image before size filter)
    """
    binary = corrected > threshold
    labeled_raw, n_raw = label(binary)

    if n_raw == 0:
        return {
            'labeled': np.zeros_like(corrected, dtype=np.int32),
            'n_particles': 0,
            'binary': binary
        }

    # Size filter
    sizes = ndimage.sum(binary, labeled_raw, range(1, n_raw + 1))
    sizes = np.array(sizes)
    valid_labels = np.where(sizes >= min_area)[0] + 1

    clean = np.isin(labeled_raw, valid_labels)
    labeled_final, n_final = label(clean)

    return {
        'labeled': labeled_final.astype(np.int32),
        'n_particles': int(n_final),
        'binary': binary
    }


def extract_properties(corrected, labeled, n_particles,
                        grid_rows=4, grid_cols=4):
    """Extract comprehensive properties for all detected particles.

    Returns a dict-of-lists (one key per property, each value is a list
    of length n_particles). This converts cleanly to an R data.frame
    via reticulate.

    Args:
        corrected: 2D float64 background-subtracted image.
        labeled: 2D int labeled particle image.
        n_particles: Number of particles.
        grid_rows, grid_cols: Tile grid dimensions.

    Returns:
        dict with keys: particle_id, centroid_y, centroid_x, tile_row,
        tile_col, local_y, local_x, area_px, max_intensity,
        mean_intensity, integrated_intensity, bbox_width, bbox_height,
        equivalent_diameter, aspect_ratio
    """
    if n_particles == 0:
        return {k: [] for k in [
            'particle_id', 'centroid_y', 'centroid_x', 'tile_row',
            'tile_col', 'local_y', 'local_x', 'area_px',
            'max_intensity', 'mean_intensity', 'integrated_intensity',
            'bbox_width', 'bbox_height', 'equivalent_diameter',
            'aspect_ratio'
        ]}

    h, w = corrected.shape
    tile_h = h // grid_rows
    tile_w = w // grid_cols

    indices = range(1, n_particles + 1)

    centroids = np.array(center_of_mass(labeled > 0, labeled, indices))
    areas = np.array(ndimage.sum(labeled > 0, labeled, indices))
    max_ints = np.array(ndimage.maximum(corrected, labeled, indices))
    mean_ints = np.array(ndimage.mean(corrected, labeled, indices))
    sum_ints = np.array(ndimage.sum(corrected, labeled, indices))
    bboxes = find_objects(labeled)

    result = {
        'particle_id': [],
        'centroid_y': [],
        'centroid_x': [],
        'tile_row': [],
        'tile_col': [],
        'local_y': [],
        'local_x': [],
        'area_px': [],
        'max_intensity': [],
        'mean_intensity': [],
        'integrated_intensity': [],
        'bbox_width': [],
        'bbox_height': [],
        'equivalent_diameter': [],
        'aspect_ratio': []
    }

    for i in range(n_particles):
        cy, cx = centroids[i]
        ty = int(min(cy / tile_h, grid_rows - 1))
        tx = int(min(cx / tile_w, grid_cols - 1))

        bbox = bboxes[i]
        bbox_h = bbox[0].stop - bbox[0].start if bbox else 0
        bbox_w = bbox[1].stop - bbox[1].start if bbox else 0

        eq_diam = 2.0 * np.sqrt(areas[i] / np.pi)
        ar = max(bbox_w, bbox_h) / max(min(bbox_w, bbox_h), 1)

        result['particle_id'].append(i + 1)
        result['centroid_y'].append(round(float(cy), 2))
        result['centroid_x'].append(round(float(cx), 2))
        result['tile_row'].append(int(ty))
        result['tile_col'].append(int(tx))
        result['local_y'].append(round(float(cy - ty * tile_h), 2))
        result['local_x'].append(round(float(cx - tx * tile_w), 2))
        result['area_px'].append(int(areas[i]))
        result['max_intensity'].append(round(float(max_ints[i]), 2))
        result['mean_intensity'].append(round(float(mean_ints[i]), 2))
        result['integrated_intensity'].append(round(float(sum_ints[i]), 2))
        result['bbox_width'].append(int(bbox_w))
        result['bbox_height'].append(int(bbox_h))
        result['equivalent_diameter'].append(round(float(eq_diam), 2))
        result['aspect_ratio'].append(round(float(ar), 2))

    return result


def auto_tune_threshold(corrected, target_count, min_area=20,
                         low=10.0, high=100.0, max_iter=25, tol=3):
    """Binary search for threshold that gives closest to target_count.

    Args:
        corrected: 2D float64 background-subtracted image.
        target_count: Desired number of particles.
        min_area: Minimum particle area in pixels.
        low, high: Search bounds for threshold.
        max_iter: Maximum binary search iterations.
        tol: Accept if |detected - target| <= tol.

    Returns:
        dict with 'threshold', 'n_particles', 'det_result'.
    """
    best_thr = (low + high) / 2
    best_diff = float('inf')
    best_result = None

    for _ in range(max_iter):
        mid = (low + high) / 2
        det = detect_particles(corrected, threshold=mid, min_area=min_area)
        n = det['n_particles']
        diff = n - target_count

        if abs(diff) < abs(best_diff):
            best_thr = mid
            best_diff = diff
            best_result = det

        if abs(diff) <= tol:
            break
        elif diff > 0:  # too many -> raise threshold
            low = mid
        else:  # too few -> lower threshold
            high = mid

    return {
        'threshold': round(best_thr, 1),
        'n_particles': best_result['n_particles'] if best_result else 0,
        'det_result': best_result
    }


def run_full_pipeline(image_path, grid_rows=4, grid_cols=4,
                       bg_sigma=30.0, clip_sigma=3.0, max_iter=10,
                       threshold=25.0, min_area=20, target_count=0):
    """Convenience function: runs the entire pipeline in one call.

    This is the main entry point from R.

    Args:
        image_path: Path to LDIR mosaic PNG.
        target_count: If > 0, auto-tune threshold to match this count.
        All other args: algorithm parameters (see individual functions).

    Returns:
        dict with:
          'particles': dict-of-lists (-> R data.frame)
          'n_particles': int
          'threshold_used': float (actual threshold used)
          'image_height': int
          'image_width': int
          'tile_height': int
          'tile_width': int
    """
    data = load_and_prepare(image_path)

    bg_result = correct_background(
        data['gray'], grid_rows, grid_cols,
        bg_sigma, clip_sigma, max_iter
    )

    if target_count > 0:
        # Auto-tune threshold to match expected particle count
        tune = auto_tune_threshold(
            bg_result['corrected'], target_count, min_area
        )
        det_result = tune['det_result']
        threshold_used = tune['threshold']
    else:
        det_result = detect_particles(
            bg_result['corrected'], threshold, min_area
        )
        threshold_used = threshold

    particles = extract_properties(
        bg_result['corrected'], det_result['labeled'],
        det_result['n_particles'], grid_rows, grid_cols
    )

    return {
        'particles': particles,
        'n_particles': det_result['n_particles'],
        'threshold_used': threshold_used,
        'image_height': data['height'],
        'image_width': data['width'],
        'tile_height': bg_result['tile_height'],
        'tile_width': bg_result['tile_width']
    }

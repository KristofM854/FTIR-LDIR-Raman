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


def apply_circle_mask(arr, cx, cy, radius, fill_value=0.0):
    """Zero out pixels outside the scan circle.

    Args:
        arr: 2D float64 array.
        cx, cy: Circle centre in pixels (float).
        radius: Circle radius in pixels (float).
        fill_value: Value written outside the circle (default 0.0).

    Returns:
        2D float64 array with pixels outside circle set to fill_value.
    """
    h, w = arr.shape
    yy, xx = np.ogrid[:h, :w]
    outside = (xx - cx) ** 2 + (yy - cy) ** 2 > radius ** 2
    out = arr.copy()
    out[outside] = fill_value
    return out


def detect_particles(corrected, threshold=25.0, min_area=10):
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


def auto_tune_threshold(corrected, target_count, min_area=10,
                         low=5.0, high=250.0, max_iter=50, tol=0):
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


def detect_scan_circle(image_path):
    """Robust scan circle detection using connected-component analysis.

    Approach:
      1. Load grayscale and apply Gaussian smoothing (suppresses particles).
      2. Threshold at the 80th percentile to isolate the bright scan disk.
      3. Morphological opening (removes bright specks outside the disk).
      4. Keep only the largest connected component.
      5. Derive centre from component centroid; radius from the 95th
         percentile of boundary-pixel distances (robust to small notches).

    This is substantially more robust than the R algebraic edge-fit because
    it is not fooled by bright single particles near the image boundary.

    Args:
        image_path: Path to LDIR image file.

    Returns:
        dict with cx (float), cy (float), r (float),
        width (int), height (int), method (str).
        On failure returns image-centre defaults.
    """
    try:
        data = load_and_prepare(image_path)
        gray = data['gray']
        h, w = gray.shape

        # Smooth heavily so particles don't bias the disk mask
        smoothed = gaussian_filter(gray, sigma=max(h, w) * 0.01)

        # Threshold: pixels above 80th percentile are "scan area"
        thresh = float(np.percentile(smoothed, 80))
        disk_mask = smoothed > thresh

        # Morphological opening: remove small bright specks
        struct = ndimage.generate_binary_structure(2, 2)
        opened = ndimage.binary_opening(disk_mask, structure=struct,
                                         iterations=max(3, int(min(h, w) * 0.005)))

        # Keep largest connected component
        labeled_disk, n_comp = label(opened)
        if n_comp == 0:
            raise ValueError("no foreground components found")
        sizes = ndimage.sum(opened, labeled_disk, range(1, n_comp + 1))
        largest_label = int(np.argmax(sizes)) + 1
        disk_only = labeled_disk == largest_label

        # Centre from centroid of the largest component
        cy_c, cx_c = center_of_mass(disk_only)

        # Radius: 95th-percentile distance of boundary pixels from centre
        # (robust to the occasional notch or clipped edge)
        eroded = ndimage.binary_erosion(disk_only, structure=struct)
        boundary = disk_only & ~eroded
        by, bx = np.where(boundary)
        if len(bx) < 20:
            raise ValueError("too few boundary pixels")
        dists = np.sqrt((bx - cx_c) ** 2 + (by - cy_c) ** 2)
        radius = float(np.percentile(dists, 95))

        return {
            'cx': float(cx_c), 'cy': float(cy_c), 'r': float(radius),
            'width': int(w), 'height': int(h),
            'method': 'python_cc'
        }

    except Exception as e:
        # Graceful fallback to image-centre defaults
        try:
            data = load_and_prepare(image_path)
            h, w = data['height'], data['width']
        except Exception:
            h, w = 0, 0
        return {
            'cx': w / 2.0, 'cy': h / 2.0, 'r': min(w, h) / 2.0 * 0.95,
            'width': int(w), 'height': int(h),
            'method': 'fallback_center'
        }


def run_full_pipeline(image_path, grid_rows=4, grid_cols=4,
                       bg_sigma=30.0, clip_sigma=3.0, max_iter=10,
                       threshold=25.0, min_area=10, target_count=0,
                       circle_cx=-1.0, circle_cy=-1.0, circle_r=-1.0):
    """Convenience function: runs the entire pipeline in one call.

    This is the main entry point from R.

    Args:
        image_path: Path to LDIR mosaic PNG.
        target_count: If > 0, auto-tune threshold to match this count.
        circle_cx, circle_cy, circle_r: Scan-circle centre and radius in
            pixels.  When circle_r > 0, pixels outside the circle are
            zeroed after background correction so they cannot be detected
            as particles.
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

    corrected = bg_result['corrected']
    if circle_r > 0:
        corrected = apply_circle_mask(corrected, circle_cx, circle_cy, circle_r)
        bg_result = dict(bg_result)          # copy so we can update
        bg_result['corrected'] = corrected

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

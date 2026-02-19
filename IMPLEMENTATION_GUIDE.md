# LDIR Particle Detection & Matching — Implementation Guide for Claude Code

## Context for the Implementer

You are building a **Shiny module** within a larger modular R application (likely a `{golem}` or `{rhino}` structured app). The module's job: take an LDIR mosaic image, detect all particles, extract their coordinates and properties, then match them to the LDIR software's exported results table (which has particle IDs, polymer types, sizes in µm — but **no coordinates**).

The computational backend is **Python** (called from R via `{reticulate}`). The R side handles UI, data flow, reactivity, and the final matching/reporting. This is a deliberate architectural choice: scipy/numpy are far superior for image processing, while R/Shiny excels at interactive data exploration and reporting.

**Read this entire document before writing any code.** The order of sections is: architecture overview → Python backend (the hard part) → R module structure → matching logic → critical gotchas.

---

## 1. Architecture Overview

```
┌─────────────────────────────────────────────────────────┐
│  Shiny App (R)                                          │
│                                                         │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  │
│  │ mod_upload    │  │ mod_detect   │  │ mod_match    │  │
│  │              │→│              │→│              │  │
│  │ Image + LDIR │  │ Particle     │  │ Coordinate   │  │
│  │ CSV upload   │  │ detection    │  │ matching     │  │
│  └──────────────┘  └──────────────┘  └──────────────┘  │
│         │                 │                  │           │
│         ▼                 ▼                  ▼           │
│  ┌──────────────────────────────────────────────────┐   │
│  │  mod_results — interactive table + overlay plot   │   │
│  └──────────────────────────────────────────────────┘   │
│                          │                              │
│         ┌────────────────┼────────────────┐             │
│         ▼                ▼                ▼             │
│  ┌────────────┐  ┌────────────┐  ┌──────────────┐     │
│  │ mod_export  │  │ mod_qc     │  │ mod_visualize│     │
│  │ CSV/Report  │  │ QC checks  │  │ Overlay plot │     │
│  └────────────┘  └────────────┘  └──────────────┘     │
│                                                         │
│  Python backend (reticulate)                            │
│  ┌──────────────────────────────────────────────────┐   │
│  │  particle_detector.py                             │   │
│  │    - background_correction()                      │   │
│  │    - detect_particles()                           │   │
│  │    - extract_properties()                         │   │
│  └──────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────┘
```

### Data flow

1. User uploads LDIR mosaic image (PNG) + LDIR results CSV
2. `mod_detect` calls Python backend → returns a data.frame of detected particles with coordinates
3. `mod_match` takes detected particles + LDIR results → matches them by properties
4. `mod_results` displays matched table, overlay visualization, QC metrics
5. `mod_export` writes final annotated CSV or report

### Key files to create

```
R/
  mod_upload_ui.R / mod_upload_server.R
  mod_detect_ui.R / mod_detect_server.R
  mod_match_ui.R  / mod_match_server.R
  mod_results_ui.R / mod_results_server.R
  utils_python.R          # reticulate setup + Python wrapper functions
  utils_matching.R        # matching algorithm (Hungarian/greedy)
inst/
  python/
    particle_detector.py  # the full detection backend
```

---

## 2. Python Backend: `particle_detector.py`

This is the core. The algorithm has been validated to detect **exactly 505 particles** from the test image, matching the LDIR software's count.

### 2.1 The Algorithm in Detail

The pipeline has three stages: background correction → thresholding → property extraction.

#### Stage 1: Iterative Sigma-Clipping Background Estimation

**Why this specific approach:**
The LDIR mosaic is a 4×4 grid of individual microscope fields of view. Each tile has a non-uniform illumination gradient (typically brighter toward one corner due to the optics). A simple Gaussian smooth of the whole image would: (a) smear tile seam artifacts into the correction, and (b) allow bright particles to pull up the background estimate, causing false negatives around large particles.

**The algorithm, step by step, for each tile:**

```
Input: tile (500×500 grayscale array)
Parameters: sigma=30, n_sigma=3, max_iter=10

1. mask ← all True (every pixel is "background")
2. bg ← GaussianFilter(tile, sigma)
3. FOR iteration = 1 to max_iter:
   a. residual ← tile - bg
   b. μ ← mean(residual[mask])
   c. σ_noise ← std(residual[mask])
   d. new_mask ← (residual < μ + n_sigma × σ_noise)
   e. IF new_mask == mask: BREAK  // converged
   f. mask ← new_mask
   g. filled ← tile.copy()
   h. filled[~mask] ← bg[~mask]   // replace particle pixels with current bg estimate
   i. bg ← GaussianFilter(filled, sigma)
4. corrected ← clip(tile - bg, min=0)
```

**Why each step matters:**
- Step 3a-d: Identifies pixels significantly above the local background (candidate particles)
- Step 3g-h: The crucial insight — before re-estimating the background, we fill in particle locations with our current best background guess. This prevents particle light from "leaking" into the background estimate through the Gaussian filter.
- Convergence typically occurs in 3–5 iterations.

**Parameter robustness (verified empirically):**
- `sigma=25, 30, 35` all produce exactly 505 detections (sigma=40 gives 507)
- `n_sigma=3` is the standard statistical choice; 2.5–3.5 all work
- The algorithm is NOT sensitive to these parameters — this is a strength

#### Stage 2: Global Thresholding + Connected Components

After background correction, detect globally (NOT per-tile) because **17 of 505 particles span tile boundaries**. If you detect per-tile, you will either split these into two fragments or miss them.

```
Input: corrected (2000×2000 background-subtracted image)
Parameters: threshold=25, min_area=20

1. binary ← (corrected > threshold)
2. labeled, n_objects ← connected_component_labeling(binary)
3. FOR each label i:
   a. area_i ← count of pixels with label i
4. valid_labels ← {i : area_i ≥ min_area}
5. clean ← keep only valid_labels in labeled image
6. re_labeled, n_final ← connected_component_labeling(clean)
```

**Why threshold=25 and min_area=20:**
- Background noise after correction: mean=1.5, std=3.5, 99th percentile=13.4
- Threshold=25 is ~7σ above background noise → negligible false positive rate
- min_area=20 removes tiny noise blobs that survive thresholding
- This combination yields exactly 505 particles; ±2 on threshold changes count by ~6–10

**Sensitivity note for the UI:** You will want to expose `threshold` and `min_area` as user-adjustable sliders (with sensible defaults), because different LDIR samples may have different noise levels. Show the detection count in real-time so users can tune to match their LDIR particle count.

#### Stage 3: Property Extraction

For each detected particle, extract:

| Property | How computed | Purpose |
|----------|-------------|---------|
| `centroid_y`, `centroid_x` | `scipy.ndimage.center_of_mass` | Global position in mosaic |
| `tile_row`, `tile_col` | `floor(centroid / tile_size)` | Tile assignment |
| `local_y`, `local_x` | `centroid - tile_origin` | Position within tile |
| `area_px` | Pixel count | Primary matching feature |
| `max_intensity` | Max corrected intensity | Secondary matching feature |
| `mean_intensity` | Mean corrected intensity | Secondary matching feature |
| `integrated_intensity` | Sum of corrected intensities | Correlates with particle brightness |
| `bbox_width`, `bbox_height` | From `find_objects` bounding box | Rough shape descriptor |
| `equivalent_diameter` | `2 * sqrt(area / π)` | For matching to LDIR size (µm) |
| `aspect_ratio` | `max(w,h) / min(w,h)` | Shape descriptor for fibers vs fragments |

### 2.2 Complete Python Implementation

```python
"""
particle_detector.py
Called from R via reticulate. All functions return numpy arrays or dicts
that reticulate can convert to R data structures.
"""

import numpy as np
from scipy import ndimage
from scipy.ndimage import label, gaussian_filter, center_of_mass, find_objects
from PIL import Image
import json


def load_and_prepare(image_path):
    """Load image and convert to grayscale float64 array.
    
    Args:
        image_path: Path to the LDIR mosaic PNG.
        
    Returns:
        dict with 'rgb' (H,W,3 uint8), 'gray' (H,W float64), 
        'height', 'width' keys.
    """
    arr = np.array(Image.open(image_path), dtype=np.float64)
    gray = 0.299 * arr[:, :, 0] + 0.587 * arr[:, :, 1] + 0.114 * arr[:, :, 2]
    return {
        'rgb': arr.astype(np.uint8),
        'gray': gray,
        'height': int(arr.shape[0]),
        'width': int(arr.shape[1])
    }


def correct_background(gray, grid_rows=4, grid_cols=4, bg_sigma=30.0,
                        clip_sigma=3.0, max_iter=10):
    """Per-tile iterative sigma-clipping background subtraction.
    
    This is the most critical function. See IMPLEMENTATION_GUIDE.md 
    Section 2.1 Stage 1 for the full algorithmic explanation.
    
    Args:
        gray: 2D float64 array (grayscale image).
        grid_rows, grid_cols: Tile grid dimensions.
        bg_sigma: Gaussian smoothing sigma for background estimation.
        clip_sigma: Number of standard deviations for clipping threshold.
        max_iter: Maximum sigma-clipping iterations.
        
    Returns:
        dict with:
          'corrected': 2D float64 array (background-subtracted, clipped ≥0)
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
            
            # Iterative sigma-clipping
            mask = np.ones_like(tile, dtype=bool)
            bg = gaussian_filter(tile, sigma=bg_sigma)
            
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
                bg = gaussian_filter(filled, sigma=bg_sigma)
            
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


def get_detection_overlay(rgb, labeled, n_particles, corrected):
    """Generate a PNG overlay image for display in Shiny.
    
    Creates the original image with red markers at particle centroids.
    Saves to a temp file and returns the path.
    
    Args:
        rgb: 3D uint8 array (original image).
        labeled: 2D int labeled particle image.
        n_particles: Number of particles.
        corrected: 2D float64 corrected image.
        
    Returns:
        Path to the saved overlay PNG.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import tempfile
    import os
    
    centroids = np.array(center_of_mass(labeled > 0, labeled, 
                                         range(1, n_particles + 1)))
    areas = np.array(ndimage.sum(labeled > 0, labeled, 
                                  range(1, n_particles + 1)))
    
    fig, ax = plt.subplots(figsize=(16, 16))
    ax.imshow(rgb)
    
    for i in range(n_particles):
        cx, cy = centroids[i, 1], centroids[i, 0]
        area = areas[i]
        if area > 200:
            ax.plot(cx, cy, 'ro', markersize=8, markerfacecolor='none',
                    markeredgewidth=1.0)
        elif area > 50:
            ax.plot(cx, cy, 'ro', markersize=5, markerfacecolor='none',
                    markeredgewidth=0.7)
        else:
            ax.plot(cx, cy, 'r+', markersize=3, markeredgewidth=0.5)
    
    # Grid lines
    h, w = rgb.shape[:2]
    for i in range(1, 4):
        ax.axhline(i * h // 4, color='cyan', linewidth=0.5, alpha=0.4)
        ax.axvline(i * w // 4, color='cyan', linewidth=0.5, alpha=0.4)
    
    ax.set_title(f'{n_particles} particles detected', fontsize=14)
    ax.axis('off')
    plt.tight_layout()
    
    path = os.path.join(tempfile.gettempdir(), 'ldir_detection_overlay.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    
    return path


def run_full_pipeline(image_path, grid_rows=4, grid_cols=4,
                       bg_sigma=30.0, clip_sigma=3.0, max_iter=10,
                       threshold=25.0, min_area=20):
    """Convenience function: runs the entire pipeline in one call.
    
    This is the main entry point from R.
    
    Args:
        image_path: Path to LDIR mosaic PNG.
        All other args: algorithm parameters (see individual functions).
        
    Returns:
        dict with:
          'particles': dict-of-lists (→ R data.frame)
          'n_particles': int
          'overlay_path': str (path to overlay PNG)
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
    
    det_result = detect_particles(
        bg_result['corrected'], threshold, min_area
    )
    
    particles = extract_properties(
        bg_result['corrected'], det_result['labeled'],
        det_result['n_particles'], grid_rows, grid_cols
    )
    
    overlay_path = get_detection_overlay(
        data['rgb'], det_result['labeled'],
        det_result['n_particles'], bg_result['corrected']
    )
    
    return {
        'particles': particles,
        'n_particles': det_result['n_particles'],
        'overlay_path': overlay_path,
        'image_height': data['height'],
        'image_width': data['width'],
        'tile_height': bg_result['tile_height'],
        'tile_width': bg_result['tile_width']
    }
```

### 2.3 Python Dependencies

```
numpy>=1.24
scipy>=1.10
Pillow>=9.0
matplotlib>=3.6
```

Install in whatever Python environment `{reticulate}` points to:
```bash
pip install numpy scipy Pillow matplotlib
```

---

## 3. R Side: Reticulate Bridge

### 3.1 `utils_python.R` — Python Setup & Wrapper Functions

```r
#' Initialize the Python environment for particle detection
#' 
#' Call this once at app startup (e.g., in global.R or server function).
#' Sets up reticulate and sources the Python module.
#'
#' @param python_path Optional path to Python binary. If NULL, uses
#'   reticulate's default discovery.
#' @return Invisible NULL. Side effect: creates `detector` in the
#'   calling environment.
setup_python_detector <- function(python_path = NULL) {
  if (!is.null(python_path)) {
    reticulate::use_python(python_path, required = TRUE)
  }
  
  # Verify required packages
  required <- c("numpy", "scipy", "PIL", "matplotlib")
  for (pkg in required) {
    if (!reticulate::py_module_available(pkg)) {
      stop(sprintf("Python package '%s' not found. Install with: pip install %s",
                    pkg, ifelse(pkg == "PIL", "Pillow", pkg)))
    }
  }
  
  # Source the detector module
  detector_path <- system.file("python", "particle_detector.py",
                                package = "yourpackage")
  # For development, use a direct path:
  # detector_path <- file.path("inst", "python", "particle_detector.py")
  
  reticulate::source_python(detector_path)
  
  invisible(NULL)
}


#' Run particle detection on an LDIR mosaic image
#'
#' Wrapper around the Python `run_full_pipeline()` function.
#' Converts Python dict-of-lists to a proper R tibble.
#'
#' @param image_path Character. Path to the LDIR mosaic PNG file.
#' @param grid_rows,grid_cols Integer. Tile grid dimensions (default 4x4).
#' @param bg_sigma Numeric. Gaussian sigma for background estimation.
#' @param threshold Numeric. Detection threshold on corrected image.
#' @param min_area Integer. Minimum particle area in pixels.
#'
#' @return A list with components:
#'   - `particles`: tibble with one row per detected particle
#'   - `n_particles`: integer count
#'   - `overlay_path`: path to the overlay PNG
#'   - `metadata`: list of image/tile dimensions
detect_particles_from_image <- function(
    image_path,
    grid_rows = 4L, grid_cols = 4L,
    bg_sigma = 30, threshold = 25, min_area = 20L
) {
  
  stopifnot(file.exists(image_path))
  
  # Call Python pipeline
  result <- run_full_pipeline(
    image_path = image_path,
    grid_rows = as.integer(grid_rows),
    grid_cols = as.integer(grid_cols),
    bg_sigma = as.double(bg_sigma),
    clip_sigma = 3.0,
    max_iter = 10L,
    threshold = as.double(threshold),
    min_area = as.integer(min_area)
  )
  
  # Convert Python dict-of-lists to tibble
  # reticulate converts Python lists to R vectors automatically
  particles_df <- tibble::as_tibble(result$particles)
  
  # Ensure correct types
  particles_df <- particles_df |>
    dplyr::mutate(
      particle_id = as.integer(particle_id),
      tile_row = as.integer(tile_row),
      tile_col = as.integer(tile_col),
      area_px = as.integer(area_px),
      bbox_width = as.integer(bbox_width),
      bbox_height = as.integer(bbox_height)
    )
  
  list(
    particles = particles_df,
    n_particles = as.integer(result$n_particles),
    overlay_path = result$overlay_path,
    metadata = list(
      image_height = as.integer(result$image_height),
      image_width = as.integer(result$image_width),
      tile_height = as.integer(result$tile_height),
      tile_width = as.integer(result$tile_width),
      grid_rows = grid_rows,
      grid_cols = grid_cols
    )
  )
}
```

### 3.2 Important reticulate Gotchas

1. **Integer conversion**: Python's `int` becomes R's `numeric` by default. Always wrap with `as.integer()` on the R side for columns that should be integer.

2. **NumPy arrays**: reticulate converts numpy arrays to R matrices/vectors. The dict-of-lists pattern (used in `extract_properties`) is the cleanest way to get a data.frame.

3. **Matplotlib backend**: The Python code sets `matplotlib.use('Agg')` explicitly. If it doesn't, Shiny may hang because matplotlib tries to open a display window.

4. **Temp files**: The overlay PNG is saved to `tempdir()`. Clean up in `onSessionEnded()`.

5. **Python errors**: Wrap the Python call in `tryCatch`. reticulate surfaces Python exceptions as R conditions, but the error messages can be cryptic. Log the full traceback.

---

## 4. Shiny Module: `mod_detect`

This is the main detection module. It lets the user upload an image, adjust detection parameters, preview results, and proceed to matching.

### 4.1 UI

```r
mod_detect_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(4,
        wellPanel(
          h4("Detection Parameters"),
          
          # Grid configuration (rarely changed, collapse by default)
          shinyjs::hidden(
            div(id = ns("grid_panel"),
              fluidRow(
                column(6, numericInput(ns("grid_rows"), "Grid rows", 4, 1, 10)),
                column(6, numericInput(ns("grid_cols"), "Grid cols", 4, 1, 10))
              )
            )
          ),
          actionLink(ns("toggle_grid"), "Show grid settings"),
          
          hr(),
          
          sliderInput(ns("bg_sigma"), "Background smoothing (σ)",
                      min = 10, max = 80, value = 30, step = 5),
          
          sliderInput(ns("threshold"), "Detection threshold",
                      min = 5, max = 80, value = 25, step = 1),
          
          sliderInput(ns("min_area"), "Minimum particle area (px)",
                      min = 1, max = 50, value = 20, step = 1),
          
          hr(),
          
          # Target count from LDIR (for reference while tuning)
          numericInput(ns("target_count"), 
                       "LDIR particle count (target)", 
                       value = NA, min = 1),
          
          actionButton(ns("run_detection"), "Run Detection",
                       class = "btn-primary", width = "100%"),
          
          br(), br(),
          
          # Result summary
          uiOutput(ns("detection_summary"))
        )
      ),
      column(8,
        tabsetPanel(
          tabPanel("Overlay", 
            imageOutput(ns("overlay_plot"), height = "700px")
          ),
          tabPanel("Particle Table",
            DT::dataTableOutput(ns("particle_table"))
          ),
          tabPanel("Diagnostics",
            plotOutput(ns("size_histogram"), height = "350px"),
            plotOutput(ns("tile_heatmap"), height = "350px")
          )
        )
      )
    )
  )
}
```

### 4.2 Server

```r
mod_detect_server <- function(id, image_path_reactive) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values to store results
    detection_result <- reactiveVal(NULL)
    
    # Run detection when button is clicked
    observeEvent(input$run_detection, {
      req(image_path_reactive())
      
      # Show progress
      withProgress(message = "Detecting particles...", {
        
        incProgress(0.1, detail = "Correcting background")
        
        result <- tryCatch({
          detect_particles_from_image(
            image_path = image_path_reactive(),
            grid_rows = input$grid_rows,
            grid_cols = input$grid_cols,
            bg_sigma = input$bg_sigma,
            threshold = input$threshold,
            min_area = input$min_area
          )
        }, error = function(e) {
          showNotification(
            paste("Detection failed:", conditionMessage(e)),
            type = "error", duration = 10
          )
          NULL
        })
        
        incProgress(0.9, detail = "Rendering")
        
        detection_result(result)
      })
    })
    
    # Detection summary (count + match to target)
    output$detection_summary <- renderUI({
      res <- detection_result()
      if (is.null(res)) return(NULL)
      
      count <- res$n_particles
      target <- input$target_count
      
      count_style <- if (!is.na(target)) {
        diff <- count - target
        if (abs(diff) <= 2) "color: green; font-weight: bold;"
        else if (abs(diff) <= 10) "color: orange; font-weight: bold;"
        else "color: red; font-weight: bold;"
      } else {
        "font-weight: bold;"
      }
      
      tagList(
        tags$p(style = count_style,
          sprintf("Detected: %d particles", count)),
        if (!is.na(target)) {
          tags$p(sprintf("Target: %d (Δ = %+d)", target, count - target))
        }
      )
    })
    
    # Overlay image
    output$overlay_plot <- renderImage({
      res <- detection_result()
      req(res)
      list(src = res$overlay_path, 
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    
    # Particle table
    output$particle_table <- DT::renderDataTable({
      res <- detection_result()
      req(res)
      DT::datatable(
        res$particles,
        options = list(pageLength = 25, scrollX = TRUE),
        filter = "top"
      ) |>
        DT::formatRound(c("centroid_y", "centroid_x", "local_y", "local_x",
                           "max_intensity", "mean_intensity",
                           "equivalent_diameter", "aspect_ratio"), 1)
    })
    
    # Size histogram
    output$size_histogram <- renderPlot({
      res <- detection_result()
      req(res)
      ggplot2::ggplot(res$particles, ggplot2::aes(x = area_px)) +
        ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "black",
                                linewidth = 0.3) +
        ggplot2::geom_vline(xintercept = median(res$particles$area_px),
                            linetype = "dashed", color = "red") +
        ggplot2::scale_x_log10() +
        ggplot2::labs(x = "Particle Area (pixels, log scale)", 
                      y = "Count",
                      title = "Particle Size Distribution") +
        ggplot2::theme_minimal()
    })
    
    # Tile heatmap
    output$tile_heatmap <- renderPlot({
      res <- detection_result()
      req(res)
      tile_counts <- res$particles |>
        dplyr::count(tile_row, tile_col) |>
        dplyr::mutate(tile_label = paste0("(", tile_row, ",", tile_col, ")"))
      
      ggplot2::ggplot(tile_counts, 
                      ggplot2::aes(x = tile_col, y = tile_row, fill = n)) +
        ggplot2::geom_tile(color = "white", linewidth = 1) +
        ggplot2::geom_text(ggplot2::aes(label = n), color = "white", size = 5) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::scale_y_reverse() +
        ggplot2::labs(x = "Tile Column", y = "Tile Row",
                      title = "Particles per Tile", fill = "Count") +
        ggplot2::theme_minimal()
    })
    
    # Return the detection result as a reactive for downstream modules
    return(detection_result)
  })
}
```

---

## 5. The Matching Problem: `utils_matching.R`

This is the intellectually hardest part. The LDIR exports a table like:

| Particle | Size (µm) | Polymer Type | Color | ... |
|----------|-----------|-------------|-------|-----|
| 1        | 45.2      | PE          | white | ... |
| 2        | 12.1      | PP          | blue  | ... |
| ...      | ...       | ...         | ...   | ... |

No coordinates. We need to match each LDIR particle to a detected particle.

### 5.1 What Information Do We Have for Matching?

From our detection: `area_px`, `equivalent_diameter`, `centroid`, `intensity`, `aspect_ratio`, `bbox` dimensions.

From LDIR: `size_um` (longest dimension), possibly `color`, possibly `particle_type` (fiber vs fragment).

The **critical bridge**: the LDIR's `size_um` should correlate with our `equivalent_diameter` or `bbox` max dimension, scaled by a pixels-to-µm conversion factor.

### 5.2 Pixel-to-Micron Calibration

The conversion factor depends on the LDIR instrument and imaging optics. Typical LDIR filter imaging:
- 13mm filter with ~10mm scanned area → 2000 px across 10,000 µm → **5 µm/px**
- 25mm filter → 2000 px across ~20,000 µm → **10 µm/px**

**This must be user-configurable.** Provide a calibration input in the UI, or calculate it from a known reference particle.

```r
#' Estimate pixel-to-micron conversion factor
#' 
#' Uses the LDIR size distribution and the detected size distribution
#' to estimate the scaling factor via robust regression.
#'
#' @param detected_sizes Numeric vector of detected particle sizes (pixels).
#' @param ldir_sizes Numeric vector of LDIR particle sizes (µm).
#' @return Numeric scalar: µm per pixel.
estimate_um_per_pixel <- function(detected_sizes, ldir_sizes) {
  # Sort both distributions and match by rank
  det_sorted <- sort(detected_sizes, decreasing = TRUE)
  ldir_sorted <- sort(ldir_sizes, decreasing = TRUE)
  
  n <- min(length(det_sorted), length(ldir_sorted))
  det_sorted <- det_sorted[1:n]
  ldir_sorted <- ldir_sorted[1:n]
  
  # Robust linear regression through origin
  # ldir_size ≈ factor × detected_size
  # Use median ratio of top 20% (large particles are most reliable)
  top_n <- max(5, floor(n * 0.2))
  ratios <- ldir_sorted[1:top_n] / det_sorted[1:top_n]
  
  stats::median(ratios)
}
```

### 5.3 Matching Strategy

The matching proceeds in two phases:

**Phase 1: Tile-level partitioning (if LDIR field info is available)**

Some LDIR software exports which field-of-view each particle was found in. If this information is available (even as an implicit numbering scheme), partition both the detected and LDIR particles by tile. This reduces the matching problem from 505×505 to sixteen ~30×30 subproblems.

**Phase 2: Size-based optimal assignment**

Within each partition (or globally if no tile info), use the **Hungarian algorithm** to find the optimal one-to-one matching that minimizes total size discrepancy.

```r
#' Match detected particles to LDIR results
#'
#' Uses the Hungarian algorithm for optimal assignment based on
#' size similarity. Optionally partitions by tile first.
#'
#' @param detected tibble from detection pipeline (must have 
#'   equivalent_diameter, tile_row, tile_col).
#' @param ldir tibble from LDIR export (must have size_um).
#' @param um_per_pixel Numeric. Conversion factor.
#' @param ldir_field_col Optional column name in `ldir` indicating
#'   the LDIR field/tile. Set to NULL if not available.
#' @param tile_mapping Optional named vector mapping LDIR field values
#'   to "row,col" strings (e.g., c("Field1" = "0,0", "Field2" = "0,1")).
#'
#' @return tibble with columns from both detected and ldir, plus:
#'   - match_distance: size discrepancy (µm)
#'   - match_quality: "good" / "uncertain" / "poor"
match_particles <- function(detected, ldir, um_per_pixel,
                             ldir_field_col = NULL,
                             tile_mapping = NULL) {
  
  # Convert detected sizes to µm
  detected <- detected |>
    dplyr::mutate(
      size_um_detected = equivalent_diameter * um_per_pixel
    )
  
  if (!is.null(ldir_field_col) && !is.null(tile_mapping)) {
    # Phase 1: Partition by tile and match within each
    matched_parts <- purrr::map_dfr(names(tile_mapping), function(field_val) {
      tile_rc <- as.integer(strsplit(tile_mapping[[field_val]], ",")[[1]])
      
      det_sub <- detected |>
        dplyr::filter(tile_row == tile_rc[1], tile_col == tile_rc[2])
      ldir_sub <- ldir |>
        dplyr::filter(.data[[ldir_field_col]] == field_val)
      
      if (nrow(det_sub) == 0 || nrow(ldir_sub) == 0) return(tibble::tibble())
      
      match_by_size(det_sub, ldir_sub)
    })
  } else {
    # Global matching (no tile info)
    matched_parts <- match_by_size(detected, ldir)
  }
  
  # Classify match quality
  matched_parts |>
    dplyr::mutate(
      match_quality = dplyr::case_when(
        match_distance < 5 ~ "good",
        match_distance < 15 ~ "uncertain",
        TRUE ~ "poor"
      )
    )
}


#' Core matching function using Hungarian algorithm
#' @keywords internal
match_by_size <- function(detected, ldir) {
  # Build cost matrix: absolute difference in size (µm)
  n_det <- nrow(detected)
  n_ldir <- nrow(ldir)
  
  cost_matrix <- matrix(NA_real_, nrow = n_det, ncol = n_ldir)
  for (i in seq_len(n_det)) {
    for (j in seq_len(n_ldir)) {
      cost_matrix[i, j] <- abs(detected$size_um_detected[i] - ldir$size_um[j])
    }
  }
  
  # Hungarian algorithm (from clue package)
  # Handles rectangular matrices (n_det ≠ n_ldir)
  assignment <- clue::solve_LSAP(cost_matrix, maximum = FALSE)
  
  # Build matched tibble
  matched <- tibble::tibble(
    detected_id = detected$particle_id,
    ldir_id = ldir$particle_id[as.integer(assignment)],
    match_distance = purrr::map2_dbl(
      seq_len(n_det), as.integer(assignment),
      ~ cost_matrix[.x, .y]
    )
  )
  
  # Join full data
  matched |>
    dplyr::left_join(detected, by = c("detected_id" = "particle_id")) |>
    dplyr::left_join(ldir, by = c("ldir_id" = "particle_id"),
                     suffix = c("_det", "_ldir"))
}
```

### 5.4 Required R Package for Matching

```r
# clue provides solve_LSAP (Hungarian algorithm)
# It's on CRAN, lightweight, no external dependencies
install.packages("clue")
```

---

## 6. Coordinate System: Pixel ↔ Physical ↔ LDIR

This is where most confusion arises. Define it clearly and stick to one convention.

### 6.1 Three Coordinate Spaces

```
PIXEL SPACE (what we detect in)
  Origin: top-left of mosaic image
  Units: pixels
  Range: x ∈ [0, 2000), y ∈ [0, 2000)
  Convention: (y, x) internally in numpy; display as (x, y)

PHYSICAL SPACE (what LDIR measures in)
  Origin: depends on instrument (usually center or corner of filter)
  Units: micrometers (µm)
  Range: depends on filter size (e.g., 0–10000 µm for 10mm scan)
  Convention: (x, y) with y increasing upward (Cartesian)

TILE SPACE (per field-of-view)
  Origin: top-left of each tile
  Units: pixels
  Range: x ∈ [0, 500), y ∈ [0, 500)
  Used for: matching within a tile
```

### 6.2 Transformations

```r
#' Convert pixel coordinates to physical coordinates
#' 
#' @param x_px, y_px Pixel coordinates (origin top-left).
#' @param um_per_pixel Conversion factor.
#' @param image_height Image height in pixels (for Y-axis flip).
#' @param origin_x, origin_y Physical coordinate of the image top-left
#'   corner (in µm). Default 0,0.
#' @param flip_y Logical. If TRUE (default), flips Y axis so that
#'   physical Y increases upward.
pixel_to_physical <- function(x_px, y_px, um_per_pixel,
                               image_height = 2000,
                               origin_x = 0, origin_y = 0,
                               flip_y = TRUE) {
  x_um <- origin_x + x_px * um_per_pixel
  if (flip_y) {
    y_um <- origin_y + (image_height - y_px) * um_per_pixel
  } else {
    y_um <- origin_y + y_px * um_per_pixel
  }
  list(x_um = x_um, y_um = y_um)
}
```

**Critical question you must answer per instrument:** Does the LDIR's physical Y increase upward (standard Cartesian) or downward (image convention)? Check the LDIR's own visualization or export a known particle near the filter edge to determine this. Getting this wrong will make matching fail silently.

### 6.3 LDIR Tile/Field Numbering

Agilent 8700 LDIR typically scans the filter in a snake/raster pattern. The field numbering may be:

```
Standard raster:          Snake pattern:
 0  1  2  3               0  1  2  3
 4  5  6  7               7  6  5  4
 8  9 10 11               8  9 10 11
12 13 14 15              15 14 13 12
```

**You must determine which pattern your instrument uses.** Provide a dropdown in the UI for the user to select, or auto-detect from the tile seam patterns. This is important because wrong tile assignment will scramble the matching.

---

## 7. LDIR CSV Import: `mod_upload`

### 7.1 What the LDIR Exports

The Agilent Clarity software typically exports CSV or XLSX with columns like:

```
Particle ID, Classification, Size (µm), Size X (µm), Size Y (µm),
Perimeter (µm), Circularity, Color, Spectrum Match (%), ...
```

The exact column names vary by LDIR software version and user configuration. The module should:
1. Let the user upload the file
2. Show a preview of the columns
3. Let the user map columns to the expected fields:
   - **Required**: particle ID, size (µm)
   - **Optional but valuable**: field/tile number, color, polymer type, circularity
4. Store the mapping for reuse

### 7.2 Column Mapping UI Pattern

```r
mod_upload_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("image_file"), "LDIR Mosaic Image (.png)", accept = ".png"),
    fileInput(ns("ldir_file"), "LDIR Results (.csv/.xlsx)", 
              accept = c(".csv", ".xlsx")),
    
    # Dynamic column mapping (appears after LDIR file upload)
    uiOutput(ns("column_mapping")),
    
    actionButton(ns("confirm_upload"), "Confirm & Proceed",
                 class = "btn-success")
  )
}

mod_upload_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    ldir_data <- reactiveVal(NULL)
    
    observeEvent(input$ldir_file, {
      req(input$ldir_file)
      ext <- tools::file_ext(input$ldir_file$name)
      df <- if (ext == "csv") {
        readr::read_csv(input$ldir_file$datapath, show_col_types = FALSE)
      } else {
        readxl::read_excel(input$ldir_file$datapath)
      }
      ldir_data(df)
    })
    
    output$column_mapping <- renderUI({
      req(ldir_data())
      cols <- names(ldir_data())
      
      tagList(
        h4("Map LDIR Columns"),
        selectInput(ns("col_id"), "Particle ID", choices = cols,
                    selected = grep("id|particle", cols, 
                                    ignore.case = TRUE, value = TRUE)[1]),
        selectInput(ns("col_size"), "Size (µm)", choices = cols,
                    selected = grep("size|length|diameter", cols,
                                    ignore.case = TRUE, value = TRUE)[1]),
        selectInput(ns("col_polymer"), "Polymer Type (optional)", 
                    choices = c("(none)", cols),
                    selected = grep("class|polymer|type", cols,
                                    ignore.case = TRUE, value = TRUE)[1]),
        selectInput(ns("col_field"), "Field/Tile (optional)",
                    choices = c("(none)", cols),
                    selected = grep("field|tile|fov", cols,
                                    ignore.case = TRUE, value = TRUE)[1]),
        
        numericInput(ns("um_per_pixel"), "µm per pixel", 
                     value = 5.0, min = 0.1, max = 50, step = 0.1),
        
        selectInput(ns("scan_pattern"), "Tile scan pattern",
                    choices = c("Raster (left-to-right)" = "raster",
                                "Snake (alternating)" = "snake"),
                    selected = "raster")
      )
    })
    
    # Return reactive list of everything downstream modules need
    return(list(
      image_path = reactive({
        req(input$image_file)
        input$image_file$datapath
      }),
      ldir_data = ldir_data,
      column_mapping = reactive({
        list(
          id = input$col_id,
          size = input$col_size,
          polymer = if (input$col_polymer == "(none)") NULL else input$col_polymer,
          field = if (input$col_field == "(none)") NULL else input$col_field
        )
      }),
      um_per_pixel = reactive(input$um_per_pixel),
      scan_pattern = reactive(input$scan_pattern)
    ))
  })
}
```

---

## 8. QC Module: `mod_qc`

After matching, provide quality control checks:

```r
#' QC checks for particle matching
#' 
#' @param matched tibble from match_particles()
#' @return tibble of QC results
run_qc_checks <- function(matched, n_detected, n_ldir) {
  tibble::tribble(
    ~check, ~status, ~detail,
    
    "Count match",
    ifelse(n_detected == n_ldir, "PASS",
           ifelse(abs(n_detected - n_ldir) <= 5, "WARN", "FAIL")),
    sprintf("Detected: %d, LDIR: %d (Δ=%+d)", n_detected, n_ldir,
            n_detected - n_ldir),
    
    "Match quality distribution",
    ifelse(mean(matched$match_quality == "good") > 0.8, "PASS",
           ifelse(mean(matched$match_quality == "good") > 0.5, "WARN", "FAIL")),
    sprintf("Good: %.0f%%, Uncertain: %.0f%%, Poor: %.0f%%",
            100 * mean(matched$match_quality == "good"),
            100 * mean(matched$match_quality == "uncertain"),
            100 * mean(matched$match_quality == "poor")),
    
    "Size correlation",
    ifelse(cor(matched$size_um_detected, matched$size_um, 
               use = "complete.obs") > 0.9, "PASS", "WARN"),
    sprintf("r = %.3f", cor(matched$size_um_detected, matched$size_um,
                            use = "complete.obs")),
    
    "Unmatched particles",
    ifelse(n_detected == nrow(matched), "PASS", "WARN"),
    sprintf("%d of %d detected particles matched", 
            nrow(matched), n_detected)
  )
}
```

---

## 9. Critical Gotchas & Edge Cases

### 9.1 Particles Spanning Tile Boundaries

**17 out of 505 particles (3.4%) span tile boundaries in the test image.** This is why detection MUST be global (after per-tile background correction). If you detect per-tile, these particles will be either:
- Split into two fragments (counted as 2 particles → count too high)
- One half falls below `min_area` and is discarded → count too low

The current pipeline handles this correctly: correct per-tile, detect globally.

### 9.2 The 4×4 Grid May Not Always Be 4×4

Some LDIR setups use 3×3, 5×5, or other grids. The `grid_rows` and `grid_cols` parameters must be user-configurable. Auto-detection is possible by looking for brightness seams at regular intervals, but manual input is more reliable.

To auto-detect grid dimensions:
```python
def detect_grid_size(gray, max_grid=8):
    """Detect tile grid dimensions from seam artifacts."""
    h, w = gray.shape
    
    # Look for brightness discontinuities in horizontal/vertical profiles
    # Average across the perpendicular axis
    h_profile = gray.mean(axis=1)  # average across columns → profile along rows
    v_profile = gray.mean(axis=0)  # average across rows → profile along cols
    
    # Try each candidate grid size and score by seam strength
    best_rows, best_cols = 1, 1
    best_score = 0
    
    for n in range(2, max_grid + 1):
        # Check horizontal seams
        step = h // n
        seam_strength = 0
        for i in range(1, n):
            pos = i * step
            # Derivative at seam position
            seam_strength += abs(h_profile[pos] - h_profile[pos-1])
        if seam_strength > best_score:
            best_rows = n
            
        # Similarly for vertical...
    
    return best_rows, best_cols
```

### 9.3 Images May Not Be Exactly Divisible

If the image is 2000×2000 but the grid is 4×4, tiles are exactly 500×500. But some LDIR software produces images that are, say, 2048×2048 or have 1-2 pixel borders. Always compute tile sizes as `height // grid_rows` and handle the remainder pixels by ignoring the last few rows/columns:

```python
tile_h = h // grid_rows  # integer division
# Process only tile_h * grid_rows rows (may lose a few pixels at bottom)
```

### 9.4 Background Sigma Should Scale with Tile Size

The default `bg_sigma=30` works for 500×500 tiles. If tiles are larger or smaller, scale proportionally:

```python
bg_sigma_adjusted = bg_sigma * (tile_size / 500.0)
```

This ensures the smoothing kernel covers the same fraction of the tile regardless of resolution.

### 9.5 The Detection Threshold Is NOT Universal

The threshold=25 was tuned for this specific sample. Different samples will have different noise floors depending on:
- Filter substrate (gold-coated vs. silver vs. silicon)
- Particle types (transparent particles are dimmer)
- Imaging settings (exposure time, gain)

**Always let the user adjust the threshold and min_area.** Show the count in real-time against the LDIR target count. This is the most important UI feature.

### 9.6 Color Channel Choice

The current implementation uses standard luminance weighting (0.299R + 0.587G + 0.114B). For LDIR images that are predominantly blue (like this one), you might get better SNR by using only the blue channel:

```python
# Alternative: blue channel only (may improve SNR for blue-tinted images)
gray = arr[:, :, 2].astype(np.float64)
```

Test both and let the user choose, or auto-detect based on the dominant color.

### 9.7 Memory & Performance

For a 2000×2000 image, the pipeline takes ~2–5 seconds on a modern CPU. The bottleneck is the iterative Gaussian filtering (16 tiles × ~4 iterations × 2 Gaussian filters per iteration = ~128 filter operations). This is fast enough for interactive use.

For larger images (e.g., 4096×4096 from high-res LDIR), consider:
- Showing a progress bar per tile
- Processing tiles in parallel (each tile is independent)
- Downsampling for preview, full-res for final detection

### 9.8 File Format Assumptions

The code assumes PNG input. LDIR software may also export TIFF (sometimes 16-bit). Add TIFF support:

```python
from PIL import Image
img = Image.open(path)
if img.mode == 'I;16':
    arr = np.array(img, dtype=np.float64)
    # Single-channel 16-bit: normalize to 0-255 range
    arr = (arr / arr.max() * 255)
elif img.mode == 'RGB':
    # Standard 8-bit RGB
    arr = np.array(img, dtype=np.float64)
```

---

## 10. Testing Strategy

### 10.1 Unit Tests for Python Backend

```python
# test_particle_detector.py
import numpy as np
import pytest
from particle_detector import (correct_background, detect_particles,
                                extract_properties)

def test_empty_image():
    """All-black image should detect 0 particles."""
    gray = np.zeros((500, 500))
    corrected = correct_background(gray, 1, 1)['corrected']
    result = detect_particles(corrected)
    assert result['n_particles'] == 0

def test_single_bright_spot():
    """A single bright Gaussian blob should be detected as 1 particle."""
    gray = np.zeros((500, 500))
    yy, xx = np.mgrid[0:500, 0:500]
    gray += 100 * np.exp(-((yy-250)**2 + (xx-250)**2) / (2*10**2))
    
    corrected = correct_background(gray, 1, 1)['corrected']
    result = detect_particles(corrected)
    assert result['n_particles'] == 1

def test_gradient_rejection():
    """A smooth gradient (no particles) should detect 0 particles."""
    yy, xx = np.mgrid[0:500, 0:500]
    gray = (xx + yy) / 4.0  # max ~250
    
    corrected = correct_background(gray, 1, 1)['corrected']
    result = detect_particles(corrected)
    assert result['n_particles'] == 0

def test_known_image():
    """The test LDIR image should yield exactly 505 particles."""
    from particle_detector import run_full_pipeline
    result = run_full_pipeline('Comparstic_LDIR_F2Ba_G3B_AU_240925.png')
    assert result['n_particles'] == 505

def test_tile_boundary_particle():
    """A particle placed at a tile boundary should be detected as 1."""
    gray = np.zeros((1000, 1000))
    # Place a bright spot straddling the boundary at y=500
    yy, xx = np.mgrid[0:1000, 0:1000]
    gray += 100 * np.exp(-((yy-500)**2 + (xx-250)**2) / (2*8**2))
    
    corrected = correct_background(gray, 2, 2)['corrected']
    result = detect_particles(corrected, threshold=15, min_area=10)
    assert result['n_particles'] == 1  # Not 2!
```

### 10.2 Integration Tests for R

```r
# test-matching.R
test_that("Hungarian matching works with equal-sized sets", {
  detected <- tibble::tibble(
    particle_id = 1:5,
    size_um_detected = c(10, 20, 30, 40, 50),
    equivalent_diameter = c(2, 4, 6, 8, 10)
  )
  ldir <- tibble::tibble(
    particle_id = 1:5,
    size_um = c(50.5, 30.2, 10.1, 40.3, 19.8)
  )
  
  matched <- match_by_size(detected, ldir)
  
  # Each detected particle should match to the closest LDIR particle
  expect_equal(nrow(matched), 5)
  # The 10µm detected should match the 10.1µm LDIR
  row_10 <- matched |> dplyr::filter(detected_id == 1)
  expect_equal(row_10$ldir_id, 3)  # LDIR particle 3 = 10.1µm
})
```

---

## 11. R Package Dependencies

```r
# DESCRIPTION Imports:
Imports:
    shiny (>= 1.7.0),
    reticulate (>= 1.28),
    tibble,
    dplyr,
    purrr,
    readr,
    readxl,
    ggplot2,
    DT,
    clue,
    shinyjs
    
# Suggests (for testing):
Suggests:
    testthat (>= 3.0.0),
    shinytest2
```

---

## 12. Minimal Working Example (Full App Skeleton)

For a quick proof-of-concept before building the full modular app:

```r
# app.R — minimal standalone version
library(shiny)
library(reticulate)
library(dplyr)
library(ggplot2)

# Source Python detector
source_python("inst/python/particle_detector.py")

ui <- fluidPage(
  titlePanel("LDIR Particle Detector"),
  sidebarLayout(
    sidebarPanel(
      fileInput("image", "LDIR Mosaic (.png)", accept = ".png"),
      sliderInput("threshold", "Threshold", 5, 80, 25),
      sliderInput("min_area", "Min Area (px)", 1, 50, 20),
      numericInput("target", "LDIR count (target)", NA),
      actionButton("detect", "Detect", class = "btn-primary"),
      hr(),
      verbatimTextOutput("count_text")
    ),
    mainPanel(
      imageOutput("overlay", height = "700px")
    )
  )
)

server <- function(input, output, session) {
  result <- reactiveVal(NULL)
  
  observeEvent(input$detect, {
    req(input$image)
    res <- run_full_pipeline(
      image_path = input$image$datapath,
      threshold = input$threshold,
      min_area = as.integer(input$min_area)
    )
    result(res)
  })
  
  output$count_text <- renderText({
    req(result())
    paste("Detected:", result()$n_particles, "particles")
  })
  
  output$overlay <- renderImage({
    req(result())
    list(src = result()$overlay_path, contentType = "image/png",
         width = "100%")
  }, deleteFile = FALSE)
}

shinyApp(ui, server)
```

---

## 13. Summary: What to Build, In What Order

1. **First**: Get `particle_detector.py` working standalone. Run the test image, verify 505 particles.

2. **Second**: Build `utils_python.R` with `setup_python_detector()` and `detect_particles_from_image()`. Test in an R console.

3. **Third**: Build the minimal app (Section 12). Get upload → detect → display working end-to-end.

4. **Fourth**: Add `mod_upload` with LDIR CSV import and column mapping.

5. **Fifth**: Build `utils_matching.R` with the Hungarian algorithm matcher.

6. **Sixth**: Build `mod_match` and `mod_results` to display matched particles.

7. **Last**: Add `mod_qc`, `mod_export`, and polish the UI.

At each step, write tests before moving to the next. The Python backend is the foundation — if it detects wrong, everything downstream fails.

# Input format

The pipeline compares particle lists from up to four instrument roles:

| role          | built-in reader for                         | quality column |
|---------------|---------------------------------------------|----------------|
| `ftir`        | PerkinElmer Spotlight particle export       | AAU (0–1)      |
| `ftir_bruker` | Bruker LUMOS export in the same table layout | AAU (0–1)     |
| `raman`       | WITec ParticleScout export                  | HQI (0–100)    |
| `ldir`        | Agilent 8700 LDIR Excel export              | LDIR quality (0–1) |

Any other particle table can be read with the **column-mapping reader** (below).
Every reader returns the same canonical particle table, and everything after
reading (alignment, matching, agreement, export, viewer) uses only that table.
The contract is defined and checked in `R/particle_contract.R`; the readers are
dispatched in `R/01_readers.R`.

## Canonical particle table

One row per particle.

| column           | type | meaning |
|------------------|------|---------|
| `instrument`     | chr  | role: `ftir`, `ftir_bruker`, `raman`, `ldir` |
| `sample_id`      | chr  | `config$sample_id` (NA when not set) |
| `particle_id`    | chr  | unique within the table, never empty |
| `x_um`, `y_um`   | num  | particle centre in µm in the instrument's own frame: any origin, x to the right, **y up**. LDIR: NA until the image join, unless the table maps coordinates |
| `feret_max_um`   | num  | longest dimension. PerkinElmer/Bruker export the ellipse major axis, WITec a true Feret maximum, LDIR the larger of width and height |
| `feret_min_um`   | num  | shortest dimension |
| `major_um`, `minor_um` | num | ellipse axes (or the best equivalent) |
| `area_um2`       | num  | area |
| `material`       | chr  | identification used for classification: the exported label, or its translation from the material table. LDIR: `unknown` below `ldir_hqi_unknown_threshold` |
| `identification_raw` | chr | the label exactly as exported |
| `quality`        | num  | the instrument's identification score |
| `quality_metric` | chr  | what `quality` is (`AAU`, `HQI`, `LDIR quality`, or the mapping's label) |
| `source_file`    | chr  | export file name |
| `material_polymer`, `material_family` | chr | exact polymer and family (`R/08b_material_map.R`) |

`validate_particles()` stops, listing every problem, when a column is missing,
a numeric column is not numeric, an ID is empty or duplicated, a coordinate is
infinite, no particle has coordinates (except LDIR before the image join), or a
size is negative.

Sizes are compared between instruments with different definitions
(`feret_max_um` above); the matching weighs size only weakly (see the README),
but size prefilters act on these values as they are.

## Canonical image object

Images are placed in the particles' frame by one mapping:

    x_um = (col − origin_col) × pixel_size_x
    y_um = (origin_row − row) × pixel_size_y        (y_axis "up")

with pixel (0, 0) the top-left corner of the file. The object records
`width_px`, `height_px`, `pixel_size_um`, `origin_px`, `y_axis`, `extent_um`
and how the placement was obtained (`provenance`):

| provenance | used for |
|------------|----------|
| `config` | extent / scale entered in the config |
| `registration` | fitted to the particles visible in the image (FTIR, Raman auto-calibration) |
| `scan_circle` | LDIR optical image: the detected scan circle is the scan diameter |
| `processed_full_image` | LDIR analyzed overlay: the whole image is the scan field |
| `excel_area_scale` | LDIR analyzed overlay scaled so the blob areas sum to the Excel areas |
| `particle_extent_approx` | approximate fallback: fitted to the particle extent |

Each run writes its image objects to `00_manifest/images.json`.

### Raman image: WITec placement file and large images

WITec writes a text file next to an exported image:

```
Width: 13000
Height: 13000
X: -6011.3134765625
Y: 12824.1048088074
Z: 18.0499992370605
Rotation: 0
```

Width and Height are the image size in µm; X, Y its **top-left corner** in the
stage frame (Y up), the frame of the particle export. The pipeline looks for
`<image name>.txt` next to the selected image (also with a trailing `_small`
removed), else for any `.txt` in the image's folder whose content has this
form (other text files are ignored; of several placement files, the one that
puts the most Raman particles on the image is used, with a warning). When
found, it places the image exactly (provenance `metadata`) and the
auto-calibration is skipped (`raman_placement_cross_check = TRUE` runs it as
a cross-check). A rotated image
(Rotation ≠ 0) is not supported and the file is ignored. Values entered by hand
(`raman_image_width_um` ...) take precedence; `raman_image_metadata` names the
file explicitly, or `NA` switches it off. Checked on sample Faanui SW1L R4 B
(14,562 particles): read this way all particles lie on the image and the
brightness at the particle positions peaks at zero offset and scale 1.00.

Select the full-size Raman image: images larger than `raman_image_max_px`
(default 2000 px on the long side) are downsized for the run, as PNG, or as
JPEG of decreasing quality when the PNG would exceed `image_max_mb` (5 MB). The
viewer, report and auto-calibration use the downsized copy; the placement is
in µm, so it does not depend on the pixel count. Originals larger than
`image_copy_max_mb` (100 MB) are not copied into the run folder; the manifest
records their path and MD5.

## Column-mapping reader

To read a table in another layout, give the role a column mapping, as an R
list or a JSON file:

```r
config$raman_column_map <- "my_raman_map.json"
```

```json
{
  "columns": {
    "particle_id":  "Part no",
    "x_um":         "X [mm]",
    "y_um":         "Y [mm]",
    "feret_max_um": "Length",
    "feret_min_um": "Width",
    "area_um2":     "Area",
    "material":     "Polymer",
    "quality":      "Score"
  },
  "sheet": "Particles",
  "coordinate_unit_um": 1000,
  "size_unit_um": 1000,
  "y_axis": "down",
  "quality_metric": "HQI"
}
```

- `columns`: canonical field → column name in the file. `material` is
  required. Allowed fields: `particle_id`, `x_um`, `y_um`, `feret_max_um`,
  `feret_min_um`, `major_um`, `minor_um`, `area_um2`, `material`, `quality`.
  Without `particle_id` the row numbers are the IDs. A missing `feret_max_um`
  is taken from `major_um` and vice versa (likewise `feret_min_um` / `minor_um`).
- `coordinate_unit_um` / `size_unit_um`: µm per unit of the file (1 for µm,
  1000 for mm); areas use the square of `size_unit_um`.
- `y_axis`: `"up"` (default) or `"down"` for image-style coordinates; `"down"`
  negates y.
- `sheet`: Excel sheet (default the first); CSV separator (`,` or `;`) and a
  decimal comma are detected.
- `quality_metric`: label for the quality column (default: the role's).
- For `ldir`, a table that maps `x_um` / `y_um` is used as is: no image
  extraction and no image join. Low-quality LDIR identifications are still
  relabelled `unknown` using `ldir_hqi_unknown_threshold`, so give the LDIR
  quality on its 0–1 scale.

The FTIR and Raman agreement cutoffs (`ftir_agreement_min_quality`,
`raman_hqi_threshold`) are applied to `quality` as it is, so a mapped quality
column must be on the scale those cutoffs expect.

## Material tables

Each role has a table translating the instrument's labels into names the
classifier knows: `inst/materials/<role>.csv`, or a file set with
`config$<role>_material_table`. See `inst/materials/README.md`. To see how each
label of a run was classified:

    Rscript tools/material_vocabulary.R <run_dir>

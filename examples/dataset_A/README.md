# Example dataset A

One filter measured with all four supported instruments. It is the example
the quick start uses (`Rscript examples/run_example_A.R`) and the regression
dataset of the test suite (`tests/testthat/test-example-dataset-A.R`).

| File | Instrument | Content |
|---|---|---|
| `FTIR_perkin.csv` | PerkinElmer Spotlight FTIR | 33 particles |
| `FTIR_image.png` | PerkinElmer Spotlight FTIR | chemical image, 2497 × 2497 px (≈ 2.0 µm/px) |
| `FTIR_bruker.csv` | Bruker (Lumos) FTIR | 42 particles, no image |
| `Raman_data.csv` | WITec Raman | 43 particles (all PET) |
| `Raman_image_small.png` | WITec Raman | micrograph, 896 × 883 px |
| `LDIR_data.xlsx` | Agilent 8700 LDIR | 40 particles (no coordinates) |
| `LDIR_image.png` | Agilent 8700 LDIR | particle map |
| `LDIR_image_analyzed.png` | Agilent 8700 LDIR | the LDIR software's processed overlay |
| `LDIR_coord_swaps.csv` | — | manual corrections of the LDIR image join (10 rows) |

The FTIR and Raman images were exported as BMP and converted losslessly to
PNG for the repository (pixel-identical; the pipeline gives byte-identical
results on both).

Reference results (current code) are recorded in
`tests/testthat/fixtures/dataset_A_fingerprint.csv`: matched pairs and their
distances, transforms, and agreement tier counts.

## Licence

The data files in this folder are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) (see `LICENSE`);
cite the repository when you use them. The code is under the MIT licence.

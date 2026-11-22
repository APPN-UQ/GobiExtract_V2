# GobiExtract_V2

Extract **per-plot reflectance and vegetation indices** from Gobi (Headwall-style) hyperspectral mosaics. Point it at one or more cubes and a plot shapefile, and it writes per-plot band means, indices, a vegetation mask and preview images for each cube.

Cleaned-up successor to the original GobiExtract notebooks: one YAML config, a resumable command-line pipeline, and documented index formulas.

## What it does

For each hyperspectral cube in the config:
1. Optionally copies the cube to a local **buffer** for speed (and safely deletes the buffer afterwards — never the original on RDM).
2. Reads the ENVI/mosaic raster once and clips it to each plot polygon from the shapefile.
3. Computes **per-plot band means** and a set of **vegetation indices** from those means (see `indices.py`).
4. Builds an **OSAVI** vegetation mask per plot (per-plot Otsu threshold, with a configurable fallback) to separate canopy from background.
5. Renders a **false-colour preview** (configurable NIR/red/green wavelengths) and a combined reflectance-spectra plot.

## Inputs (`config.yaml`)

- `hyperspec_paths` — list of 1..N cubes (ENVI `.hdr`/`.dat` or a mosaic). Raw cubes live on the RDM.
- `shp_path` — one plot-boundary shapefile for all cubes.
- `id_field` (+ optional `id_field2`) — plot ID column(s) in the shapefile.
- `out_root` — output base; results go to `<out_root>/<project_name>/<date>/`.
- `fc_targets_nm` — NIR/red/green wavelengths for the false-colour preview; `nearest_wl_tolerance_nm` — band-match tolerance.
- `plot_scale`, `png_dpi` — preview sizing; `osavi_otsu_fallback` — threshold fallback.
- `copy_to_buffer` / `delete_buffer_when_done` / `buffer_dirname` — local-buffer behaviour.
- `debug` / `debug_max_plots` — limit plots for a quick test run.

## Run

```bash
conda env create -f environment.yml   # geopandas, rasterio, scikit-image, matplotlib, tqdm, pyyaml, ...
conda activate <env>
python hyperspec_pipeline.py -c config.yaml
```

Resumable: each finished cube writes a success marker keyed to a signature of the cube + shapefile + relevant config, so re-running skips cubes that are already done and unchanged.

`Extraction_V6.ipynb` is the interactive/exploratory version of the same workflow.

## Outputs (per cube, under `<out_root>/<project>/<date>/`)

- `band_means.csv` — per-plot mean reflectance per band.
- `veg_indices_from_means.csv` — per-plot vegetation indices.
- false-colour PNG, OSAVI GeoTIFF + PNG (vegetation mask), and a combined reflectance-spectra PNG.

## Indices (`indices.py`)

`veg_index_dict` holds each index as `{formula, bands (target nm), doi}` — NDVI, GNDVI, NDVI705, EVI, OSAVI, MCARI and more (~14), each with a literature DOI. Add an index by adding an entry; bands are matched to the nearest wavelength within `nearest_wl_tolerance_nm`.

## Notes

- Originals on the RDM are never modified or deleted — only the local buffer is cleaned up.
- See the brain `registries/repos.md` (project `dairy-hyperspectral`) for context; raw cubes live under the project's RDM `Gobi/<date>/...gpro/products/`.

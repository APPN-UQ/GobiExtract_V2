#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Hyperspectral → per-plot reflectance & indices
- Supports multiple cubes via a YAML config
- Uses a local buffer for speed and safely deletes it (never the original)
- Single-pass raster read per cube (avoid reopening repeatedly)
- Saves false-color PNG, OSAVI GeoTIFF + PNG, band_means.csv, veg_indices_from_means.csv
- Saves combined reflectance plot PNG
- Pretty nested progress bars (3 levels) with clear handover between cubes
"""

from __future__ import annotations

import os
import re
import sys
import shutil
import inspect
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence, Dict, Any, List

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterio.enums import Resampling
from rasterio.transform import array_bounds
from affine import Affine

import matplotlib.pyplot as plt
from tqdm.auto import tqdm
from skimage.filters import threshold_otsu

# Local index formulas (expects {"INDEX": {"bands": {...}, "formula": callable}})
from indices import veg_index_dict

# -----------------------------
# Config dataclass & loader
# -----------------------------
@dataclass
class Config:
    project_name: str
    hyperspec_paths: List[str]  # list of cubes
    shp_path: str
    id_field: str
    id_field2: Optional[str] = None

    # I/O and runtime
    out_root: str = "."                 # base folder to create <project>/<date> subfolders
    copy_to_buffer: bool = True
    delete_buffer_when_done: bool = True
    buffer_dirname: str = "_buffer"
    debug: bool = False
    debug_max_plots: int = 10

    buffer_copy_policy: str = "mtime_size"   # or "exists"
    skip_if_complete: bool = True
    force_rerun: bool = False
    success_marker: str = "_SUCCESS.json"


    # Preview/processing options
    plot_scale: int = 4                 # downsample factor for previews
    nearest_wl_tolerance_nm: float = 8.0

    # False-color targets (nm)
    fc_targets_nm: Dict[str, float] = None  # {"nir": 840, "red": 660, "green": 560}

    # OSAVI settings (bands from veg_index_dict; formula pulled there)
    # Fallback threshold if Otsu fails for per-plot masking
    osavi_otsu_fallback: float = 0.2

    # PNG DPI for quicklooks / plots
    png_dpi: int = 200

    def __post_init__(self):
        if self.fc_targets_nm is None:
            self.fc_targets_nm = {"nir": 840, "red": 660, "green": 560}


def load_config(cfg_path: Path) -> Config:
    import yaml  # keep local to avoid import cost if not needed
    with open(cfg_path, "r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)
    return Config(**raw)


# -----------------------------
# Helpers
# -----------------------------
PIPELINE_VERSION = "2025-10-15-r1"

def _file_sig(p: Path) -> dict:
    try:
        st = p.stat()
        return {"path": str(p), "size": st.st_size, "mtime": int(st.st_mtime)}
    except Exception:
        return {"path": str(p), "size": None, "mtime": None}

def _cfg_sig_subset(cfg: Config) -> dict:
    # Only include settings that actually affect outputs
    return {
        "project_name": cfg.project_name,
        "id_field": cfg.id_field,
        "id_field2": cfg.id_field2,
        "plot_scale": cfg.plot_scale,
        "nearest_wl_tolerance_nm": cfg.nearest_wl_tolerance_nm,
        "fc_targets_nm": cfg.fc_targets_nm,
        "osavi_otsu_fallback": cfg.osavi_otsu_fallback,
        "png_dpi": cfg.png_dpi,
        # reading/runtime does not change output content:
        # "copy_to_buffer": cfg.copy_to_buffer, etc. (intentionally omitted)
    }

def make_signature(cube_path: Path, shp_path: Path, cfg: Config) -> dict:
    import hashlib, json
    sig_obj = {
        "pipeline_version": PIPELINE_VERSION,
        "cube": _file_sig(cube_path),
        "shapefile": _file_sig(shp_path),
        "cfg": _cfg_sig_subset(cfg),
    }
    blob = json.dumps(sig_obj, sort_keys=True).encode("utf-8")
    sig_obj["sha1"] = hashlib.sha1(blob).hexdigest()
    return sig_obj

def write_success_marker(marker_path: Path, signature: dict):
    import json, time
    marker_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"finished_at": int(time.time()), "signature": signature}
    marker_path.write_text(json.dumps(payload, indent=2))

def marker_matches(marker_path: Path, signature: dict) -> bool:
    import json
    try:
        payload = json.loads(marker_path.read_text())
        return payload.get("signature", {}).get("sha1") == signature.get("sha1")
    except Exception:
        return False


def buffer_up_to_date(src: Path, dst: Path, policy: str = "mtime_size") -> bool:
    """
    'mtime_size' -> skip copy if same size AND dst mtime >= src mtime (current behavior)
    'exists'     -> skip copy if destination exists (fastest; good when you don't expect source changes)
    """
    if not dst.exists():
        return False
    if policy == "exists":
        return True
    # default 'mtime_size'
    s1, s2 = src.stat(), dst.stat()
    return (s1.st_size == s2.st_size) and (s2.st_mtime >= s1.st_mtime)



def copy_file_with_progress(src: Path, dst: Path, chunk_size=16 * 1024 * 1024, pb_desc: str = ""):
    dst.parent.mkdir(parents=True, exist_ok=True)
    total = src.stat().st_size
    with open(src, "rb") as fsrc, open(dst, "wb") as fdst, tqdm(
        total=total, unit="B", unit_scale=True, unit_divisor=1024, desc=pb_desc, leave=False, position=2
    ) as pbar:
        while True:
            buf = fsrc.read(chunk_size)
            if not buf:
                break
            fdst.write(buf)
            pbar.update(len(buf))
    shutil.copystat(src, dst)

def delete_local_buffer_safe(buffer_dir: Path, original_src: Path):
    """
    Safely delete local buffer. Never touch original source.
    Extra guardrails + robust Windows delete with retry.
    """
    if not buffer_dir.exists():
        return

    # Guardrails
    if buffer_dir.name != "_buffer":
        tqdm.write(f"[SAFE-DEL] Refusing to delete: buffer_dir name is '{buffer_dir.name}', expected '_buffer'.")
        return
    if "APPF_Data-Q6940" in str(buffer_dir):
        tqdm.write("[SAFE-DEL] Buffer appears to be on the network — aborting deletion.")
        return

    try:
        orig_drive = original_src.drive.lower()
    except Exception:
        orig_drive = ""
    try:
        buf_drive = buffer_dir.drive.lower()
    except Exception:
        buf_drive = ""
    if orig_drive and buf_drive and (orig_drive == buf_drive):
        tqdm.write("[SAFE-DEL] Buffer and original appear on the same drive — skipping for safety.")
        return

    import stat, time

    def _onerror(func, path, exc_info):
        # Clear read-only and retry once
        try:
            os.chmod(path, stat.S_IWRITE)
            func(path)
        except Exception:
            pass  # swallow; we'll retry at top level

    # Try delete with small backoff
    for attempt in range(3):
        try:
            shutil.rmtree(buffer_dir, onerror=_onerror)
            tqdm.write(f"[SAFE-DEL] Deleted local buffer: {buffer_dir}")
            return
        except Exception as e:
            time.sleep(0.5)
    tqdm.write(f"[SAFE-DEL] Could not delete buffer (locked?): {buffer_dir}")


def data_path_for_envi(hdr_path: Path) -> Path:
    """
    Given an ENVI .hdr path, return the matching data file path that GDAL expects.
    Tries 'data file = ...' in the header; falls back to same-stem .bin/.img/.dat.
    """
    try:
        txt = hdr_path.read_text(errors="ignore")
    except Exception:
        txt = ""
    m = re.search(r"data\s*file\s*=\s*([^\n\r]+)", txt, flags=re.IGNORECASE)
    if m:
        name = m.group(1).strip().strip("{}").strip()
        cand = (hdr_path.parent / name)
        if cand.exists():
            return cand

    # Fallback: same stem with common ENVI data extensions
    for ext in (".bin", ".img", ".dat"):
        cand = hdr_path.with_suffix(ext)
        if cand.exists():
            return cand

    raise FileNotFoundError(
        f"Could not locate ENVI data file for header: {hdr_path}\n"
        "Check the 'data file =' entry inside the .hdr, or ensure a same-stem .bin/.img/.dat exists."
    )


def nearest_band_1based(wavelengths: np.ndarray, target_nm: float) -> int:
    return int(np.argmin(np.abs(wavelengths - target_nm))) + 1


def percentile_stretch(img: np.ndarray, p_low=2, p_high=98) -> np.ndarray:
    finite = np.isfinite(img)
    if not finite.any():
        return np.zeros_like(img, dtype=np.float32)
    lo, hi = np.percentile(img[finite], [p_low, p_high])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo, hi = np.nanmin(img), np.nanmax(img)
        if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
            lo, hi = 0.0, 1.0
    return np.clip((img - lo) / (hi - lo + 1e-12), 0, 1).astype("float32")

def resolve_envi_pair(p: Path) -> tuple[Path, Path]:
    """
    Return (hdr_path, data_path) for an ENVI dataset.
    Accepts .hdr OR data (.bin/.dat/.img).
    Looks up 'data file =' inside the .hdr if needed.
    """
    suffix = p.suffix.lower()

    # If user gave the header
    if suffix == ".hdr":
        hdr = p
        try:
            txt = hdr.read_text(errors="ignore")
        except Exception:
            txt = ""
        m = re.search(r"data\s*file\s*=\s*([^\n\r]+)", txt, flags=re.IGNORECASE)
        if m:
            name = m.group(1).strip().strip("{}").strip()
            data = hdr.parent / name
            if data.exists():
                return hdr, data
        # Fallback: same stem with common data extensions
        for ext in (".bin", ".dat", ".img"):
            cand = hdr.with_suffix(ext)
            if cand.exists():
                return hdr, cand
        raise FileNotFoundError(f"Could not locate ENVI data file for header: {hdr}")

    # If user gave a data file (.bin/.dat/.img)
    if suffix in {"", ".bin", ".dat", ".img"}:
        data = p
        hdr = p.with_suffix(".hdr")
        if hdr.exists():
            return hdr, data
        # Search for any header in the folder whose 'data file =' points to this data filename
        for hdr_cand in p.parent.glob("*.hdr"):
            try:
                txt = hdr_cand.read_text(errors="ignore")
                m = re.search(r"data\s*file\s*=\s*([^\n\r]+)", txt, flags=re.IGNORECASE)
                if m and m.group(1).strip().strip("{}").strip() == p.name:
                    return hdr_cand, data
            except Exception:
                pass
        raise FileNotFoundError(f"Could not find matching .hdr for data file: {p}")

    # Anything else: not an ENVI file we handle here
    raise ValueError(f"Unsupported ENVI path extension: {p.suffix}")


def detect_reflectance_scale(arr_like: np.ndarray) -> float:
    p99 = np.nanpercentile(arr_like, 99)
    return 10000.0 if p99 > 5.0 else 1.0


def series_for_param(index_name: str, param_name: str, ref_df: pd.DataFrame,
                     band_cols: List[str], available_wls: np.ndarray,
                     tol_nm: float) -> Optional[pd.Series]:
    """
    Map a formula parameter to a reflectance Series in 'ref_df'.
      1) If param looks like r###, use nearest
      2) Else look into veg_index_dict bands map
      3) Else None (non-spectral param with default)
    """
    # r###
    m = re.fullmatch(r"r(\d+)", param_name, flags=re.IGNORECASE)
    if m:
        target = float(m.group(1))
        idx = int(np.argmin(np.abs(available_wls - target)))
        wl = available_wls[idx]
        col = band_cols[idx]
        return ref_df[col] if abs(wl - target) <= tol_nm else None

    # semantic
    bands = veg_index_dict[index_name].get("bands", {})
    if param_name in bands:
        target = float(bands[param_name])
        idx = int(np.argmin(np.abs(available_wls - target)))
        wl = available_wls[idx]
        col = band_cols[idx]
        return ref_df[col] if abs(wl - target) <= tol_nm else None

    return None


# -----------------------------
# Core processing per-cube
# -----------------------------
def process_cube(cfg: Config, shp: gpd.GeoDataFrame, cube_path: Path, pbar_lvl2_desc: str = ""):
    """
    Process one hyperspectral cube...
    """
    # ---------- Resolve ENVI pair & output folders (by data filename) ----------
    hdr_path, data_path = resolve_envi_pair(cube_path)
    cube_name = data_path.stem  # <-- use the DATA file stem so names are consistent
    out_base  = Path(cfg.out_root) / cfg.project_name / cube_name
    out_data  = out_base / "data"
    out_plots = out_base / "plots"
    buffer_dir = out_base / cfg.buffer_dirname
    out_data.mkdir(parents=True, exist_ok=True)
    out_plots.mkdir(parents=True, exist_ok=True)

    # ---------- Level 2 progress bar ----------
    sub = tqdm(total=7, desc=pbar_lvl2_desc or "Processing steps", leave=False, position=1)

    # ---------- Optional local buffer copy ----------
    work_data = data_path  # this is what we'll open with rasterio
    try:
        if cfg.copy_to_buffer:
            buffer_dir.mkdir(parents=True, exist_ok=True)
            local_hdr  = buffer_dir / hdr_path.name
            local_data = buffer_dir / data_path.name

            # copy header
            if not buffer_up_to_date(hdr_path, local_hdr, cfg.buffer_copy_policy):
                copy_file_with_progress(hdr_path, local_hdr, pb_desc=f"Copy {hdr_path.name}")
            else:
                tqdm.write(f"[BUFFER] Up-to-date: {local_hdr}")

            # copy data (.bin/.dat/.img)
            if not buffer_up_to_date(data_path, local_data, cfg.buffer_copy_policy):
                copy_file_with_progress(data_path, local_data, pb_desc=f"Copy {data_path.name}")
            else:
                tqdm.write(f"[BUFFER] Up-to-date: {local_data}")

            # copy sidecars that share the HEADER stem
            src_stem = hdr_path.stem
            for sidecar in hdr_path.parent.glob(f"{src_stem}.*"):
                if sidecar.name in {hdr_path.name, data_path.name}:
                    continue
                target = buffer_dir / sidecar.name
                if not buffer_up_to_date(sidecar, target, cfg.buffer_copy_policy):
                    try:
                        copy_file_with_progress(sidecar, target, pb_desc=f"Copy {sidecar.name}")
                    except PermissionError:
                        shutil.copy2(sidecar, target)
                        tqdm.write(f"[BUFFER] Copied (no progress) → {target}")
                else:
                    tqdm.write(f"[BUFFER] Up-to-date: {target}")

            # use local copy for processing
            work_data = local_data

        sub.update(1)

        # ---------- Open once ----------
        sub.set_description_str("Open cube")
        with rasterio.open(work_data) as src:  # <-- open work_data (not work_path)
            hs_crs = src.crs
            shp_use = shp.to_crs(hs_crs) if shp.crs != hs_crs else shp

            # ---------- Wavelengths ----------
            sub.set_description_str("Read wavelengths")
            env = src.tags(ns="ENVI")
            if "wavelength" in env:
                try:
                    wavelengths = np.array([float(w) for w in env["wavelength"].strip("{}").split(",")])
                except Exception:
                    wavelengths = np.linspace(400, 1000, src.count)
                if wavelengths.size != src.count:
                    wavelengths = np.linspace(400, 1000, src.count)
            else:
                wavelengths = np.linspace(400, 1000, src.count)
            sub.update(1)

            # ---------- False-color quicklook ----------
            sub.set_description_str("False-color preview")
            idx_fc = [
                nearest_band_1based(wavelengths, cfg.fc_targets_nm["nir"]),   # R <- NIR
                nearest_band_1based(wavelengths, cfg.fc_targets_nm["red"]),   # G <- Red
                nearest_band_1based(wavelengths, cfg.fc_targets_nm["green"]), # B <- Green
            ]
            out_h = max(1, src.height // max(1, cfg.plot_scale))
            out_w = max(1, src.width  // max(1, cfg.plot_scale))
            fc_small = src.read(indexes=idx_fc, out_shape=(3, out_h, out_w), resampling=Resampling.average).astype("float32")
            R = percentile_stretch(fc_small[0]); G = percentile_stretch(fc_small[1]); B = percentile_stretch(fc_small[2])
            fc_rgb = np.dstack([R, G, B])

            scale_y = src.height / out_h
            scale_x = src.width  / out_w
            transform_scaled = src.transform * Affine.scale(scale_x, scale_y)
            xmin, ymin, xmax, ymax = array_bounds(out_h, out_w, transform_scaled)

            fig, ax = plt.subplots(figsize=(10, 10))
            ax.imshow(fc_rgb, extent=(xmin, xmax, ymin, ymax), origin="upper")
            ax.set_title(f"False-color (NIR-Red-Green)\nRaster CRS: {hs_crs} | Shapefile CRS: {shp_use.crs}")
            try:
                shp_use.boundary.plot(ax=ax, facecolor="none", edgecolor="yellow", linewidth=0.8)
            except Exception as e:
                tqdm.write(f"[WARN] Overlay shapefile on false-color failed: {e}")
            ax.set_aspect("equal"); ax.set_xlabel("Easting / X"); ax.set_ylabel("Northing / Y")
            falsecolor_png = out_plots / f"falsecolor_{cube_name}.png"
            plt.savefig(falsecolor_png, dpi=cfg.png_dpi, bbox_inches="tight"); plt.close(fig)
            sub.update(1)

            # ---------- OSAVI full-res + mask ----------
            sub.set_description_str("Compute OSAVI")
            red_nm = veg_index_dict["OSAVI"]["bands"]["red"]
            nir_nm = veg_index_dict["OSAVI"]["bands"]["nir"]
            osavi_fn = veg_index_dict["OSAVI"]["formula"]

            # rasterio is 1-based band indexing
            idx_red = nearest_band_1based(wavelengths, red_nm)
            idx_nir = nearest_band_1based(wavelengths, nir_nm)

            red = src.read(idx_red, masked=True).astype("float32")
            nir = src.read(idx_nir, masked=True).astype("float32")

            valid = ~(red.mask | nir.mask) & np.isfinite(red.data) & np.isfinite(nir.data)
            osavi = np.full(red.shape, np.nan, dtype="float32")
            if valid.any():
                scale_r = detect_reflectance_scale(red.data[valid])
                scale_n = detect_reflectance_scale(nir.data[valid])
                r_data = red.data / scale_r
                n_data = nir.data / scale_n
                osavi_vals = osavi_fn(n_data[valid], r_data[valid]).astype("float32")
                osavi[valid] = osavi_vals
            else:
                raise RuntimeError("OSAVI has no valid pixels — check inputs.")

            thr_global = threshold_otsu(osavi[np.isfinite(osavi)])
            veg_mask = np.zeros(osavi.shape, dtype=np.uint8)
            finite = np.isfinite(osavi)
            veg_mask[finite] = (osavi[finite] > thr_global).astype(np.uint8)

            # Save GeoTIFFs
            prof_mask = src.profile.copy(); prof_mask.update(count=1, dtype="uint8", nodata=0)
            mask_tif = out_plots / f"osavi_mask_{cube_name}.tif"
            with rasterio.open(mask_tif, "w", **prof_mask) as dst:
                dst.write(veg_mask, 1)

            prof_osavi = src.profile.copy(); prof_osavi.update(count=1, dtype="float32", nodata=np.float32(np.nan))
            osavi_tif = out_plots / f"osavi_{cube_name}.tif"
            with rasterio.open(osavi_tif, "w", **prof_osavi) as dst:
                dst.write(osavi.astype("float32"), 1)

            # OSAVI preview PNG
            red_small = src.read(idx_red, out_shape=(out_h, out_w), resampling=Resampling.average).astype("float32")
            nir_small = src.read(idx_nir, out_shape=(out_h, out_w), resampling=Resampling.average).astype("float32")
            valid_small = np.isfinite(red_small) & np.isfinite(nir_small)
            red_small /= detect_reflectance_scale(red_small[valid_small]) if valid_small.any() else 1.0
            nir_small /= detect_reflectance_scale(nir_small[valid_small]) if valid_small.any() else 1.0
            osavi_small = np.full((out_h, out_w), np.nan, dtype="float32")
            osavi_small[valid_small] = osavi_fn(nir_small[valid_small], red_small[valid_small]).astype("float32")
            veg_mask_small = np.zeros_like(osavi_small, dtype=bool)
            finite_small = np.isfinite(osavi_small)
            veg_mask_small[finite_small] = osavi_small[finite_small] > thr_global

            fig, ax = plt.subplots(figsize=(10, 10))
            ax.imshow(veg_mask_small, extent=(xmin, xmax, ymin, ymax), origin="upper")
            ax.set_title(f"OSAVI Vegetation Mask (> Otsu={thr_global:.4f})\nRaster CRS: {hs_crs} | Shapefile CRS: {shp_use.crs}")
            try:
                shp_use.boundary.plot(ax=ax, facecolor="none", edgecolor="yellow", linewidth=0.8)
            except Exception as e:
                tqdm.write(f"[WARN] Overlay shapefile on mask failed: {e}")
            ax.set_aspect("equal"); ax.set_xlabel("Easting / X"); ax.set_ylabel("Northing / Y")
            mask_png = out_plots / f"osavi_mask_{cube_name}.png"
            plt.savefig(mask_png, dpi=cfg.png_dpi, bbox_inches="tight"); plt.close(fig)
            sub.update(1)

            # ---------- Per-plot band means (FAST: windowed + chunked) ----------
            sub.set_description_str("Band means per plot")
            band_labels = [f"r{int(round(wl))}" for wl in wavelengths]
            records, skipped = [], []

            # one-based for rasterio.read
            red_idx1 = nearest_band_1based(wavelengths, red_nm)
            nir_idx1 = nearest_band_1based(wavelengths, nir_nm)
            osavi_fn_plot = veg_index_dict["OSAVI"]["formula"]

            from rasterio.windows import from_bounds
            from rasterio.windows import transform as win_transform
            from rasterio.features import geometry_mask as rio_geom_mask

            B_CHUNK = 32  # bands per read

            n_plots_total = len(shp_use)
            n_plots = min(n_plots_total, cfg.debug_max_plots) if cfg.debug else n_plots_total
            pbar_plots = tqdm(total=n_plots, desc=f"Plots 0/{n_plots}", leave=False, position=2)

            try:
                for i, (_, row) in enumerate(shp_use.iterrows()):
                    if cfg.debug and i >= n_plots:
                        break

                    pid = f"{row[cfg.id_field]}" if cfg.id_field2 is None else f"{row[cfg.id_field]}_{row[cfg.id_field2]}"
                    geom = row.geometry

                    # Window around this geometry
                    minx, miny, maxx, maxy = geom.bounds
                    win = from_bounds(minx, miny, maxx, maxy, transform=src.transform).round_offsets().round_lengths()

                    if win.width <= 0 or win.height <= 0:
                        skipped.append(pid)
                        pbar_plots.update(1); pbar_plots.set_description(f"Plots {min(i+1, n_plots)}/{n_plots}")
                        continue

                    w_transform = win_transform(win, src.transform)
                    out_shape = (int(win.height), int(win.width))

                    # Geometry mask (True=outside polygon). We want inside → invert later.
                    gmask_out = rio_geom_mask([geom], out_shape=out_shape, transform=w_transform,
                                              invert=False, all_touched=True)
                    gmask_in = ~gmask_out

                    # Read red & NIR window to build OSAVI mask
                    red_w = src.read(red_idx1, window=win, masked=True).astype("float32")
                    nir_w = src.read(nir_idx1, window=win, masked=True).astype("float32")

                    valid2 = ~(red_w.mask | nir_w.mask)
                    if not np.any(valid2):
                        skipped.append(pid)
                        pbar_plots.update(1); pbar_plots.set_description(f"Plots {min(i+1, n_plots)}/{n_plots}")
                        continue

                    r = red_w.data.copy()
                    n = nir_w.data.copy()
                    # autoscale if needed
                    if (np.nanpercentile(n[valid2], 99) > 5) or (np.nanpercentile(r[valid2], 99) > 5):
                        r /= 10000.0
                        n /= 10000.0

                    osavi_local = np.full(r.shape, np.nan, dtype="float32")
                    osavi_local[valid2] = osavi_fn_plot(n[valid2], r[valid2]).astype("float32")

                    finite_local = np.isfinite(osavi_local)
                    if np.count_nonzero(finite_local) < 50:
                        osavi_mask_2d = valid2
                    else:
                        vals = osavi_local[finite_local]
                        try:
                            thr_local = threshold_otsu(vals)
                        except Exception:
                            thr_local = cfg.osavi_otsu_fallback
                        osavi_mask_2d = np.zeros_like(finite_local, dtype=bool)
                        osavi_mask_2d[finite_local] = vals > thr_local
                        if osavi_mask_2d.sum() < 20:
                            osavi_mask_2d = valid2

                    # Final per-pixel mask to use for means:
                    sel_mask = osavi_mask_2d & gmask_in

                    rec = {"Plot_ID": pid}

                    # Read all bands in chunks over the same window and compute means using sel_mask
                    total_bands = src.count
                    label_idx = 0
                    for start in range(1, total_bands + 1, B_CHUNK):
                        end = min(start + B_CHUNK - 1, total_bands)
                        idxs = list(range(start, end + 1))
                        chunk = src.read(indexes=idxs, window=win, masked=True).astype("float32")  # (B, H, W)

                        # Build a 3D mask: selected pixels and not band.mask
                        chunk_mask = np.broadcast_to(sel_mask[None, :, :], chunk.shape) & (~chunk.mask)
                        data = chunk.data

                        for b_off, band_idx in enumerate(idxs):
                            vals = data[b_off][chunk_mask[b_off]]
                            rec[band_labels[label_idx]] = float(np.nanmean(vals)) if vals.size else None
                            label_idx += 1

                    records.append(rec)

                    pbar_plots.update(1)
                    pbar_plots.set_description(f"Plots {min(i+1, n_plots)}/{n_plots}")

            except KeyboardInterrupt:
                # Save partial progress so you don't lose work
                partial_csv = out_data / "band_means_partial.csv"
                if records:
                    pd.DataFrame(records).to_csv(partial_csv, index=False)
                    tqdm.write(f"[INTERRUPTED] Wrote partial results → {partial_csv}")
                raise

            pbar_plots.close()
            band_means_df = pd.DataFrame(records)
            band_means_csv = out_data / "band_means.csv"
            band_means_df.to_csv(band_means_csv, index=False)
            if skipped:
                tqdm.write(f"[INFO] Skipped {len(skipped)} plots (no overlap/valid pixels): {skipped[:6]}{' ...' if len(skipped) > 6 else ''}")
            sub.update(1)

            # ---------- Reflectance plot ----------
            sub.set_description_str("Reflectance plot")
            wl_cols, wls = [], []
            for c in band_means_df.columns:
                m = re.fullmatch(r"r(\d+)", str(c))
                if m:
                    wl_cols.append(c); wls.append(int(m.group(1)))
            if not wl_cols:
                raise RuntimeError("No r### columns found in band_means.csv.")

            order = np.argsort(wls)
            wls_sorted = np.asarray(wls)[order]
            cols_sorted = [wl_cols[i] for i in order]
            Y = band_means_df[cols_sorted].to_numpy(dtype=float)

            scale_all = detect_reflectance_scale(Y)
            Y_scaled = Y / scale_all

            fig, ax = plt.subplots(figsize=(11, 6))
            cmap = plt.cm.viridis
            n_lines = Y_scaled.shape[0]
            for i in range(n_lines):
                y = Y_scaled[i, :]
                valid = np.isfinite(y)
                if valid.sum() < 5:
                    continue
                color = cmap(i / max(1, n_lines - 1))
                ax.plot(wls_sorted[valid], y[valid], linewidth=1.2, alpha=0.6, color=color)

            ax.set_xlabel("Wavelength (nm)")
            ax.set_ylabel("Reflectance (unitless)")
            ax.set_title(f"Reflectance curves — {cfg.project_name} ({cube_name})")
            ax.grid(True, alpha=0.25)
            refl_png = out_plots / f"reflectance_curves_{cube_name}.png"
            plt.savefig(refl_png, dpi=cfg.png_dpi, bbox_inches="tight"); plt.close(fig)

            # ---------- Vegetation indices from mean reflectances ----------
            sub.set_description_str("Veg indices from means")
            df_means = band_means_df.copy()
            plot_id_col = "Plot_ID"
            band_cols = [c for c in df_means.columns if re.fullmatch(r"r\d+", str(c))]
            avail_wls = np.array([int(re.findall(r"\d+", c)[0]) for c in band_cols], dtype=float)

            p99_all = float(np.nanpercentile(pd.to_numeric(df_means[band_cols].values.ravel(), errors="coerce"), 99))
            needs_scale = p99_all > 5.0
            ref = df_means[band_cols].astype("float32").copy()
            if needs_scale:
                ref /= 10000.0

            out_df = pd.DataFrame({plot_id_col: df_means[plot_id_col]})
            for idx_name, spec in veg_index_dict.items():
                fn = spec["formula"]
                sig = inspect.signature(fn)

                args = []
                missing = False
                for p_name, p_info in sig.parameters.items():
                    if p_info.default is not inspect._empty:
                        continue
                    s = series_for_param(idx_name, p_name, ref, band_cols, avail_wls, cfg.nearest_wl_tolerance_nm)
                    if s is None:
                        missing = True
                        break
                    args.append(s.astype("float32"))

                if missing:
                    out_df[idx_name] = np.nan
                    continue

                try:
                    vals = fn(*args)
                except Exception as e:
                    tqdm.write(f"[WARN] {idx_name} failed ({e}); fill NaN.")
                    vals = np.full(len(out_df), np.nan, dtype="float32")

                vals = np.asarray(vals, dtype="float32")
                vals[~np.isfinite(vals)] = np.nan
                out_df[idx_name] = vals

            out_indices_csv = out_data / "veg_indices_from_means.csv"
            out_df.to_csv(out_indices_csv, index=False)
            sub.update(1)

        # --- with rasterio.open ends here ---
    
        sig = make_signature(cube_path, Path(cfg.shp_path), cfg)
        write_success_marker(out_base / cfg.success_marker, sig)

    finally:
        sub.close()
        if cfg.copy_to_buffer and cfg.delete_buffer_when_done:
            delete_local_buffer_safe(buffer_dir, cube_path)

# -----------------------------
# Main
# -----------------------------
def main():
    import argparse
    parser = argparse.ArgumentParser(description="Hyperspectral per-plot extraction & indices")
    parser.add_argument("-c", "--config", required=True, help="Path to YAML configuration file")
    args = parser.parse_args()

    cfg = load_config(Path(args.config))

    # Validate inputs
    hypers = [Path(p) for p in cfg.hyperspec_paths]
    for hp in hypers:
        if not hp.exists():
            raise FileNotFoundError(f"Hyperspectral cube not found: {hp}")

    shp_path = Path(cfg.shp_path)
    if not shp_path.exists():
        raise FileNotFoundError(f"Shapefile not found: {shp_path}")

    # Load shapefile once
    shp = gpd.read_file(shp_path)
    if cfg.id_field not in shp.columns:
        raise KeyError(f"id_field '{cfg.id_field}' not found in shapefile attributes.")
    if (cfg.id_field2 is not None) and (cfg.id_field2 not in shp.columns):
        raise KeyError(f"id_field2 '{cfg.id_field2}' not found in shapefile attributes.")

    # Level 1 progress over cubes
    n_cubes = len(hypers)
    pbar_cubes = tqdm(total=n_cubes, desc=f"Cubes 0/{n_cubes}", position=0)

    for i, cube_path in enumerate(hypers, start=1):
        # Clear/close any residual level 2/3 bars by creating fresh ones per cube (leave=False)
        pbar_cubes.set_description(f"Cubes {i}/{n_cubes}")
        tqdm.write(f"[START] {cube_path}")

        cube_name = cube_path.stem
        out_base = Path(cfg.out_root) / cfg.project_name / cube_name
        marker = out_base / cfg.success_marker
        sig = make_signature(cube_path, Path(cfg.shp_path), cfg)

        if cfg.skip_if_complete and (not cfg.force_rerun) and marker.exists() and marker_matches(marker, sig):
            tqdm.write(f"[SKIP] Complete and up-to-date → {out_base}")
            pbar_cubes.update(1)
            continue

        process_cube(cfg, shp, cube_path, pbar_lvl2_desc="Processing steps")

        pbar_cubes.update(1)

    pbar_cubes.close()
    tqdm.write("[DONE] All cubes processed.")


if __name__ == "__main__":
    main()

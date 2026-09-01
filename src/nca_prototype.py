"""SEEA EA-informed EO screening prototype for Mol, Belgium.

The script discovers public ESA WorldCover and Sentinel-2 data through STAC,
aligns them on a 10 m grid, derives transparent spectral indicators, and writes
prototype ecosystem extent and condition tables. It is a portfolio workflow,
not an official ecosystem account.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
DEFAULT_BBOX = (4.95, 51.15, 5.25, 51.30)
DEFAULT_DATE_RANGE = "2024-06-01/2024-08-31"

WORLD_COVER_CLASSES = {
    10: "Tree cover",
    20: "Shrubland",
    30: "Grassland",
    40: "Cropland",
    50: "Built-up",
    60: "Bare or sparse vegetation",
    70: "Snow and ice",
    80: "Permanent water bodies",
    90: "Herbaceous wetland",
    95: "Mangroves",
    100: "Moss and lichen",
}

WORLD_COVER_COLORS = {
    10: "#006400",
    20: "#ffbb22",
    30: "#ffff4c",
    40: "#f096ff",
    50: "#fa0000",
    60: "#b4b4b4",
    70: "#f0f0f0",
    80: "#0064c8",
    90: "#0096a0",
    95: "#00cf75",
    100: "#fae6a0",
}


@dataclass(frozen=True)
class AnalysisGrid:
    crs: object
    transform: object
    width: int
    height: int
    pixel_area_ha: float


def safe_normalized_difference(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Calculate (a-b)/(a+b), returning NaN where the denominator is zero."""
    a = np.asarray(a, dtype="float32")
    b = np.asarray(b, dtype="float32")
    result = np.full(a.shape, np.nan, dtype="float32")
    denominator = a + b
    valid = np.isfinite(a) & np.isfinite(b) & (np.abs(denominator) > 1e-6)
    result[valid] = (a[valid] - b[valid]) / denominator[valid]
    return result


def calculate_extent_account(
    land_cover: np.ndarray,
    pixel_area_ha: float,
) -> pd.DataFrame:
    """Summarise mapped WorldCover extent in an accounting-style table."""
    valid = np.isin(land_cover, list(WORLD_COVER_CLASSES))
    total_pixels = int(valid.sum())
    rows: list[dict[str, object]] = []
    for code, name in WORLD_COVER_CLASSES.items():
        count = int((land_cover == code).sum())
        if count == 0:
            continue
        rows.append(
            {
                "worldcover_code": code,
                "ecosystem_proxy": name,
                "pixel_count": count,
                "extent_ha": round(count * pixel_area_ha, 2),
                "share_of_mapped_area_pct": round(100 * count / total_pixels, 3),
            }
        )
    return pd.DataFrame(rows).sort_values("extent_ha", ascending=False).reset_index(drop=True)


def exploratory_condition_index(ndvi: np.ndarray, ndmi: np.ndarray) -> np.ndarray:
    """Create a transparent, uncalibrated 0-100 EO screening index."""
    ndvi_scaled = np.clip((ndvi + 0.10) / 0.90, 0, 1)
    ndmi_scaled = np.clip((ndmi + 0.40) / 1.00, 0, 1)
    score = 100 * (0.65 * ndvi_scaled + 0.35 * ndmi_scaled)
    score[~(np.isfinite(ndvi) & np.isfinite(ndmi))] = np.nan
    return score.astype("float32")


def calculate_condition_table(
    land_cover: np.ndarray,
    valid_observation: np.ndarray,
    ndvi: np.ndarray,
    ndmi: np.ndarray,
    condition_index: np.ndarray,
) -> pd.DataFrame:
    """Summarise EO indicators and observation support by mapped type."""
    rows: list[dict[str, object]] = []
    for code, name in WORLD_COVER_CLASSES.items():
        class_mask = land_cover == code
        class_pixels = int(class_mask.sum())
        if class_pixels == 0:
            continue
        indicator_mask = class_mask & valid_observation & np.isfinite(ndvi) & np.isfinite(ndmi)
        valid_pixels = int(indicator_mask.sum())
        if valid_pixels:
            ndvi_values = ndvi[indicator_mask]
            ndmi_values = ndmi[indicator_mask]
            score_values = condition_index[indicator_mask]
            metrics = {
                "median_ndvi": round(float(np.nanmedian(ndvi_values)), 4),
                "ndvi_p25": round(float(np.nanpercentile(ndvi_values, 25)), 4),
                "ndvi_p75": round(float(np.nanpercentile(ndvi_values, 75)), 4),
                "median_ndmi": round(float(np.nanmedian(ndmi_values)), 4),
                "ndmi_p25": round(float(np.nanpercentile(ndmi_values, 25)), 4),
                "ndmi_p75": round(float(np.nanpercentile(ndmi_values, 75)), 4),
                "exploratory_condition_index_0_100": round(float(np.nanmedian(score_values)), 2),
            }
        else:
            metrics = {
                "median_ndvi": np.nan,
                "ndvi_p25": np.nan,
                "ndvi_p75": np.nan,
                "median_ndmi": np.nan,
                "ndmi_p25": np.nan,
                "ndmi_p75": np.nan,
                "exploratory_condition_index_0_100": np.nan,
            }
        rows.append(
            {
                "worldcover_code": code,
                "ecosystem_proxy": name,
                "mapped_pixel_count": class_pixels,
                "valid_indicator_pixel_count": valid_pixels,
                "valid_observation_coverage_pct": round(100 * valid_pixels / class_pixels, 2),
                **metrics,
            }
        )
    return pd.DataFrame(rows).sort_values("mapped_pixel_count", ascending=False).reset_index(drop=True)


def _catalog() -> Any:
    import planetary_computer as pc
    from pystac_client import Client

    return Client.open(STAC_URL, modifier=pc.sign_inplace)


def _select_sentinel_item(catalog: Any, bbox: tuple[float, ...], date_range: str):
    items = list(
        catalog.search(
            collections=["sentinel-2-l2a"],
            bbox=bbox,
            datetime=date_range,
            query={"eo:cloud_cover": {"lt": 15}},
            max_items=100,
        ).items()
    )
    if not items:
        raise RuntimeError("No Sentinel-2 items met the search criteria")
    return min(items, key=lambda item: float(item.properties.get("eo:cloud_cover", 100)))


def _select_worldcover_item(catalog: Any, bbox: tuple[float, ...]):
    items = list(catalog.search(collections=["esa-worldcover"], bbox=bbox, max_items=20).items())
    preferred = [item for item in items if "2021_v200" in item.id]
    if preferred:
        return preferred[0]
    if not items:
        raise RuntimeError("No ESA WorldCover item intersects the study area")
    return items[0]


def _read_reference_band(href: str, bbox: tuple[float, ...]) -> tuple[np.ndarray, AnalysisGrid]:
    import rasterio
    from rasterio.warp import transform_bounds
    from rasterio.windows import from_bounds

    with rasterio.open(href) as src:
        projected_bounds = transform_bounds("EPSG:4326", src.crs, *bbox, densify_pts=21)
        window = from_bounds(*projected_bounds, transform=src.transform).round_offsets().round_lengths()
        data = src.read(1, window=window).astype("float32")
        transform = src.window_transform(window)
        pixel_area_ha = abs(transform.a * transform.e) / 10000
        grid = AnalysisGrid(
            crs=src.crs,
            transform=transform,
            width=data.shape[1],
            height=data.shape[0],
            pixel_area_ha=pixel_area_ha,
        )
    return data, grid


def _read_to_grid(href: str, grid: AnalysisGrid, resampling: Any) -> np.ndarray:
    import rasterio
    from rasterio.vrt import WarpedVRT

    with rasterio.open(href) as src:
        with WarpedVRT(
            src,
            crs=grid.crs,
            transform=grid.transform,
            width=grid.width,
            height=grid.height,
            resampling=resampling,
        ) as vrt:
            return vrt.read(1)


def _create_figure(
    output_path: Path,
    blue: np.ndarray,
    green: np.ndarray,
    red: np.ndarray,
    land_cover: np.ndarray,
    ndvi: np.ndarray,
    condition_index: np.ndarray,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap

    rgb = np.dstack([red, green, blue])
    finite_rgb = rgb[np.isfinite(rgb)]
    upper = float(np.nanpercentile(finite_rgb, 98)) if finite_rgb.size else 0.3
    rgb = np.clip(rgb / max(upper, 0.05), 0, 1) ** 0.85

    codes = [code for code in WORLD_COVER_CLASSES if np.any(land_cover == code)]
    colours = [WORLD_COVER_COLORS[code] for code in codes]
    positions = {code: index for index, code in enumerate(codes)}
    categorical = np.full(land_cover.shape, np.nan, dtype="float32")
    for code, index in positions.items():
        categorical[land_cover == code] = index

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    axes[0, 0].imshow(rgb)
    axes[0, 0].set_title("Sentinel-2 true-colour context")

    cmap = ListedColormap(colours)
    bounds = np.arange(-0.5, len(codes) + 0.5, 1)
    axes[0, 1].imshow(categorical, cmap=cmap, norm=BoundaryNorm(bounds, cmap.N))
    axes[0, 1].set_title("ESA WorldCover ecosystem proxies")
    legend_handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", color=colour, label=WORLD_COVER_CLASSES[code])
        for code, colour in zip(codes, colours)
    ]
    axes[0, 1].legend(handles=legend_handles, loc="lower left", fontsize=7, framealpha=0.85)

    ndvi_plot = axes[1, 0].imshow(ndvi, cmap="RdYlGn", vmin=-0.2, vmax=0.9)
    axes[1, 0].set_title("Sentinel-2 NDVI")
    fig.colorbar(ndvi_plot, ax=axes[1, 0], fraction=0.046, pad=0.04)

    condition_plot = axes[1, 1].imshow(condition_index, cmap="viridis", vmin=0, vmax=100)
    axes[1, 1].set_title("Exploratory EO condition index (0-100)")
    fig.colorbar(condition_plot, ax=axes[1, 1], fraction=0.046, pad=0.04)

    for axis in axes.flat:
        axis.set_xticks([])
        axis.set_yticks([])
    fig.suptitle("EO-based ecosystem extent and condition screening around Mol, Belgium", fontsize=15)
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def run(output_dir: Path, bbox: tuple[float, ...], date_range: str, unsafe_ssl: bool) -> None:
    import rasterio
    from rasterio.enums import Resampling

    output_dir.mkdir(parents=True, exist_ok=True)
    catalog = _catalog()
    sentinel_item = _select_sentinel_item(catalog, bbox, date_range)
    worldcover_item = _select_worldcover_item(catalog, bbox)

    env_options = {
        "GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR",
        "GDAL_HTTP_MULTIRANGE": "YES",
        "VSI_CACHE": "TRUE",
    }
    if unsafe_ssl:
        env_options["GDAL_HTTP_UNSAFESSL"] = "YES"

    with rasterio.Env(**env_options):
        red_raw, grid = _read_reference_band(sentinel_item.assets["B04"].href, bbox)
        green_raw = _read_to_grid(sentinel_item.assets["B03"].href, grid, Resampling.bilinear)
        blue_raw = _read_to_grid(sentinel_item.assets["B02"].href, grid, Resampling.bilinear)
        nir_raw = _read_to_grid(sentinel_item.assets["B08"].href, grid, Resampling.bilinear)
        swir_raw = _read_to_grid(sentinel_item.assets["B11"].href, grid, Resampling.bilinear)
        scene_class = _read_to_grid(sentinel_item.assets["SCL"].href, grid, Resampling.nearest)
        land_cover = _read_to_grid(worldcover_item.assets["map"].href, grid, Resampling.nearest).astype("uint8")

    red = red_raw / 10000.0
    green = green_raw.astype("float32") / 10000.0
    blue = blue_raw.astype("float32") / 10000.0
    nir = nir_raw.astype("float32") / 10000.0
    swir = swir_raw.astype("float32") / 10000.0

    excluded_scl = np.array([0, 1, 3, 8, 9, 10, 11])
    valid_observation = (~np.isin(scene_class, excluded_scl)) & (red_raw > 0) & (nir_raw > 0)
    ndvi = safe_normalized_difference(nir, red)
    ndmi = safe_normalized_difference(nir, swir)
    ndvi[~valid_observation] = np.nan
    ndmi[~valid_observation] = np.nan
    condition_index = exploratory_condition_index(ndvi, ndmi)

    extent_table = calculate_extent_account(land_cover, grid.pixel_area_ha)
    condition_table = calculate_condition_table(
        land_cover,
        valid_observation,
        ndvi,
        ndmi,
        condition_index,
    )
    extent_table.to_csv(output_dir / "ecosystem_extent_account.csv", index=False)
    condition_table.to_csv(output_dir / "ecosystem_condition_indicators.csv", index=False)
    _create_figure(
        output_dir / "summary_map.png",
        blue,
        green,
        red,
        land_cover,
        ndvi,
        condition_index,
    )

    metadata = {
        "project": "EO-Based Ecosystem Extent and Condition Screening Around Mol, Belgium",
        "run_date": date.today().isoformat(),
        "study_area_bbox_wgs84": list(bbox),
        "analysis_grid": {
            "crs": str(grid.crs),
            "width": grid.width,
            "height": grid.height,
            "pixel_area_ha": grid.pixel_area_ha,
        },
        "sentinel_2": {
            "collection": "sentinel-2-l2a",
            "item_id": sentinel_item.id,
            "acquisition_datetime": sentinel_item.datetime.isoformat() if sentinel_item.datetime else None,
            "catalogue_cloud_cover_pct": sentinel_item.properties.get("eo:cloud_cover"),
            "valid_observation_coverage_pct": round(100 * float(valid_observation.mean()), 2),
            "bands": ["B02", "B03", "B04", "B08", "B11", "SCL"],
        },
        "worldcover": {
            "collection": "esa-worldcover",
            "item_id": worldcover_item.id,
            "asset": "map",
        },
        "processing": {
            "date_search_range": date_range,
            "cloud_filter_pct": 15,
            "target_resolution_m": 10,
            "excluded_scl_classes": excluded_scl.tolist(),
            "ndvi_formula": "(NIR-Red)/(NIR+Red)",
            "ndmi_formula": "(NIR-SWIR1)/(NIR+SWIR1)",
            "exploratory_condition_index": "100 * (0.65 * clip((NDVI+0.10)/0.90) + 0.35 * clip((NDMI+0.40)/1.00))",
        },
        "limitations": [
            "WorldCover classes are land-cover proxies, not an authoritative ecosystem or habitat typology.",
            "One Sentinel-2 date does not represent annual ecosystem condition.",
            "NDVI and NDMI do not directly measure biodiversity, ecological integrity, or ecosystem services.",
            "The exploratory condition index is uncalibrated and must not be used for policy ranking.",
            "A production account requires authoritative typologies, multi-temporal indicators, reference data, and uncertainty assessment.",
        ],
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print("Selected Sentinel-2 item:", sentinel_item.id)
    print("Selected WorldCover item:", worldcover_item.id)
    print("Valid Sentinel-2 coverage (%):", metadata["sentinel_2"]["valid_observation_coverage_pct"])
    print("Mapped area (ha):", round(float(extent_table["extent_ha"].sum()), 2))
    print("Outputs written to:", output_dir.resolve())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parents[1] / "outputs")
    parser.add_argument("--bbox", type=float, nargs=4, default=DEFAULT_BBOX, metavar=("WEST", "SOUTH", "EAST", "NORTH"))
    parser.add_argument("--date-range", default=DEFAULT_DATE_RANGE)
    parser.add_argument("--unsafe-ssl", action="store_true", help="Allow GDAL HTTPS reads without certificate verification")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.output_dir, tuple(args.bbox), args.date_range, args.unsafe_ssl)

# EO-Based Ecosystem Extent and Condition Screening Around Mol, Belgium

This compact portfolio project demonstrates how open Earth Observation data can support the upstream geospatial work required for Natural Capital Accounting. It creates a prototype ecosystem extent account and a set of ecosystem-condition indicators for a study area around Mol, Belgium.

The workflow is informed by the structure of the UN System of Environmental-Economic Accounting - Ecosystem Accounting (SEEA EA), but it is **not an official ecosystem account**. ESA WorldCover classes are used as ecosystem-type proxies, and Sentinel-2 spectral indicators are used as transparent screening variables rather than validated ecological condition metrics.

## Questions addressed

1. What broad ecosystem and land-cover types occur in the study area, and what is their mapped extent?
2. What are the observed distributions of vegetation greenness and canopy moisture within each mapped type?
3. How complete and reliable is the valid Sentinel-2 observation support?
4. How can the workflow remain reproducible, quality-controlled, and ready for extension with authoritative habitat, field, biodiversity, or ecosystem-service data?

## Data

- **ESA WorldCover 2021 v2.0**, 10 m: broad land-cover classes used as ecosystem-type proxies.
- **Sentinel-2 Level-2A**, 10/20 m: blue, green, red, near-infrared, shortwave-infrared, and scene-classification bands.
- Both datasets are discovered through the Microsoft Planetary Computer STAC API. Exact item identifiers and acquisition metadata are recorded in `outputs/run_metadata.json`.

## Workflow

1. Search the STAC catalog for data intersecting the Mol study area.
2. Select a low-cloud summer Sentinel-2 observation.
3. Align all bands and WorldCover to one 10 m analysis grid.
4. Apply Sentinel-2 scene-classification quality masks.
5. Calculate NDVI and NDMI.
6. Produce an ecosystem extent table in hectares and percent.
7. Summarise condition indicators and valid-observation coverage by mapped type.
8. Generate a transparent, exploratory EO condition index and a four-panel summary figure.
9. Write provenance, parameters, and limitations to machine-readable metadata.

## Run

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python src/nca_prototype.py
python -m unittest discover -s tests -v
```

In managed environments where GDAL cannot locate the local certificate chain, the script can be run with `--unsafe-ssl`. This affects only public Cloud Optimized GeoTIFF reads and should not be used when a valid certificate configuration is available.

## Outputs

- `ecosystem_extent_account.csv`: mapped extent by WorldCover ecosystem proxy.
- `ecosystem_condition_indicators.csv`: NDVI, NDMI, observation coverage, and exploratory condition index by mapped type.
- `summary_map.png`: true-colour context, ecosystem proxy, NDVI, and condition-index panels.
- `run_metadata.json`: data lineage, study area, selected assets, processing parameters, and limitations.

The repository contains the derived tables, metadata, and summary figure only.
No downloaded Sentinel-2 or WorldCover raster files are committed; the script
reads the required public Cloud Optimized GeoTIFF windows through STAC.

## Completed run results

The validated run used Sentinel-2 Level-2A data acquired on **25 June 2024** (catalogue cloud cover: **0.81%**) and ESA WorldCover 2021 v2.0. The aligned 10 m grid covered **37,013.76 ha**, with **99.04%** valid Sentinel-2 observation support after scene-classification masking.

The largest mapped ecosystem proxies were:

| Ecosystem proxy | Extent (ha) | Share of mapped area |
| --- | ---: | ---: |
| Tree cover | 15,301.56 | 41.34% |
| Grassland | 8,579.07 | 23.18% |
| Cropland | 7,680.67 | 20.75% |
| Built-up | 4,192.51 | 11.33% |
| Permanent water bodies | 1,198.93 | 3.24% |

Tree cover had a median NDVI of **0.4945** and median NDMI of **0.1995**; grassland had a median NDVI of **0.4469** and median NDMI of **0.1249**. These are observation-specific spectral summaries, not policy-ready ecological condition assessments.

## Interpretation limits

- WorldCover is a land-cover product, not an authoritative ecosystem or habitat typology.
- A single low-cloud Sentinel-2 acquisition does not represent annual ecosystem condition.
- NDVI and NDMI are partial spectral indicators; they do not directly measure biodiversity, ecological integrity, or ecosystem-service supply.
- The exploratory condition index is intentionally simple and uncalibrated. It is included to demonstrate a reproducible indicator workflow, not to rank ecological quality for policy decisions.
- A production SEEA EA workflow should use nationally accepted ecosystem typologies, multi-temporal indicators, uncertainty assessment, field/reference data, and stakeholder-agreed condition variables.

## Extension path

The prototype is designed to be extended with EUNIS habitat data, Copernicus Land Monitoring Service products, protected-area boundaries, Essential Biodiversity Variables, field observations, ecosystem-service models, uncertainty propagation, and annual change accounts.

## License

The source code is released under the [MIT License](LICENSE).

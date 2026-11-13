# Data sources and distributed tables

## Original inventory databases

Download the field databases in Access format from MITECO:

- [Third National Forest Inventory (IFN3)](https://www.miteco.gob.es/es/biodiversidad/servicios/banco-datos-naturaleza/informacion-disponible/ifn3_bbdd_descargas_htm.html).
- [Fourth National Forest Inventory (IFN4)](https://www.miteco.gob.es/es/biodiversidad/temas/inventarios-nacionales/inventario-forestal-nacional/cuarto_inventario.html).

These official pages also provide documentation. Run `python scripts/download_inventories.py` to retrieve and extract the study’s database coverage automatically. `provenance/inventory_sources.json` lists all 42 source URLs, target filenames, sizes, province coverage, and SHA-256 hashes from the original study inputs. The downloader enforces those hashes and rejects changed publisher files. See REPRODUCING.md for options and storage requirements.

## Plot-level forecasts

The tables in `analysis/results/plot_level/` contain one row per held-out plot and case. Plot IDs are identifiers within the inventory, not geographic coordinates. Biomass columns are **Mg of aboveground dry biomass per hectare**, not carbon or CO₂.

### `age_only_forecasts.csv`

| Column | Meaning |
|---|---|
| `case` | Species and, for maritime pine, regional case |
| `plot_id` | Plot identifier; read as text |
| `province_code` | Province code or name, depending on the source analysis; read as text |
| `age` | Recorded stand age in years |
| `observed` | Inventory-derived standing aboveground biomass |
| `default` | Forecast using the default stand assumptions |
| `measured_inputs` | Forecast using measured height and density inputs |
| `reference` | Median biomass of the fitting-province plots, used as a simple reference |
| `direct_age` | Direct-age indicator where exported; blank for cases without this exported field |

### `measured_tree_forecasts.csv`

| Column | Meaning |
|---|---|
| `case`, `plot_id`, `province_code` | Case and plot identifiers |
| `years` | Interval between inventory measurements |
| `agb3` | IFN3 biomass, also the no-change reference forecast |
| `pred_no_entry` | Forecast without adding trees entering the measured set |
| `pred_total` | Forecast including the entry component |
| `agb4` | IFN4 reference biomass |

`young_mediterranean_reference.json` contains the fitting and held-out counts, fitting median biomass, and reference RMSE for the young Mediterranean maritime-pine comparison.

## Other outputs

The remaining CSV and JSON files under `results/` are evaluation outputs and intermediate summaries. Some Markdown notes describe earlier screens or exploratory analyses; the plot-level files above are the final figure inputs in this snapshot.

## Map background

`analysis/figures/data/spain_provinces.geojson` is a Spain subset of Natural Earth 1:10m administrative boundaries, in ETRS89 / UTM zone 30N. It is a contextual map, not inventory plot locations. Natural Earth data are [public domain](https://www.naturalearthdata.com/about/terms-of-use/).

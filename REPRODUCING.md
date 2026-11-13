# Reproducing the analyses

## From the distributed forecasts

Follow the README quick start. It runs independently of the original AI4DIS repository. The distributed CSV files are sufficient to recompute descriptive error metrics and regenerate the four figures. The scripts read their inputs relative to their own location.

## From the inventory databases

Install the Python dependencies in `requirements-lock.txt` and install `mdbtools` (for example, `sudo apt install mdbtools` on Ubuntu). The code reads Access databases using `mdb-export`; Microsoft Access is not required.

Set an external data directory:

```bash
export IFN_DATA_ROOT="$HOME/ifn-data"
# Optional when mdb-export is not on PATH:
# export MDB_EXPORT=/path/to/mdb-export
```

Download the exact study coverage before running the analyses:

```bash
python scripts/download_inventories.py
```

This script needs only Python’s standard library. It uses the URLs, destination filenames, sizes, and study-input SHA-256 hashes in `provenance/inventory_sources.json`. It downloads ZIP archives to `$IFN_DATA_ROOT/downloads/`, extracts only the selected database from each archive, normalizes the expected filename (including the Cataluña/Cataluna spelling), and verifies the database before installing it. Existing correct files are skipped. Existing mismatched files are preserved and reported as errors. Downloads retry up to three times; partial downloads are never installed.

The default covers 42 databases (19 IFN4 files and 23 IFN3 provinces), about 552 MB uncompressed plus ZIP cache space. Allow at least 1 GB for these downloads, and additional space for generated caches. The 42-file download was tested against MITECO on 13 November 2026; every downloaded database matched the original study hash.

```bash
# Inspect all URLs and expected destinations without downloading.
python scripts/download_inventories.py --list

# Check an existing download, without network access or writes.
python scripts/download_inventories.py --verify-only

# Download selected manifest entries or use a different destination.
python scripts/download_inventories.py --ids ifn3-p02 --data-root /path/to/ifn-data
```

Use the full default set for reproducing the published coverage. A download receipt is written to `$IFN_DATA_ROOT/inventory_download_receipt.json`. If the publisher updates a database, the script fails rather than silently substituting a new version; consult the recorded source page or obtain the matching archived input. Do not update expected hashes merely to bypass a mismatch.

The resulting layout is:

```text
$IFN_DATA_ROOT/
├── ifn3/
│   ├── p02/Ifn3p02.accdb
│   ├── p06/Ifn3p06.accdb
│   └── pXX/Ifn3pXX.accdb
├── ifn4/
│   ├── Ifn4p15-coruna.accdb
│   ├── Ifn4p27-lugo.accdb
│   ├── Ifn4p32-ourense.accdb
│   ├── Ifn4p36-pontevedra.accdb
│   ├── Ifn4_Albacete.accdb
│   ├── Ifn4_Extremadura.accdb
│   └── ...
└── ifn_validation/                 # generated linked-tree caches
```

The filenames above illustrate the layout; the manifest lists the complete download set. Use a dedicated data directory: the four-other-pine analyses scan the available IFN4 databases, so adding provinces can change the evaluation set.

Run these commands from the repository root:

```bash
# Build strict tree links and their measurement intervals.
python models/ifn3_ifn4_pairs.py

# Rebuild the four-other-pine age-screen table used by the evaluation.
python analysis/audit_other_pines.py

# Recompute held-out plot forecasts, using the original evaluation functions.
python analysis/export_plot_level.py

# Recompute descriptive metrics and draw the figures.
python scripts/summarize_results.py
MPLBACKEND=Agg python analysis/figures/make_figures.py
```

These commands overwrite generated results. The plot export also invokes evaluation functions that write JSON summaries. The original analysis scripts retain their fit/test definitions and model parameters. Saved supporting summaries are included under `results/`; this sequence is focused on the final plot forecasts, not every exploratory experiment in those historical outputs.

Some shared modules retain IFN2 readers because they also provide utilities used by the IFN3–IFN4 code. IFN2 is not required for the quick start or the final plot-level export above.

## Scope of verification

The package includes a source snapshot manifest with original and packaged hashes. Path adaptations make inventory locations configurable; they do not intentionally change scientific calculations. `provenance/verification.json` records the checks performed when preparing this snapshot.

The pinned requirements capture the six direct Python dependencies from the working analysis environment, not a complete operating-system lockfile. Fonts, GDAL/PROJ, and plotting-library dependencies can affect rendering. The figure generator uses a non-interactive backend when `MPLBACKEND=Agg` is set.

The code and associated software documentation use the MIT licence in `LICENSE`. Derived CSV/JSON data and research figures use CC BY 4.0 in `LICENSE-DATA`. Exact original database hashes, sizes, coverage, and download URLs are recorded in `provenance/inventory_sources.json`. The original raw database bytes remain external inputs.

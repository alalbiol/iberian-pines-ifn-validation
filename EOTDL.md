# Asset 1: Iberian pine inventory-validation package

This is a **dataset entry with accompanying reproduction code**, named
`iberian-pines-ifn-validation`. It contains 364 age-only and 1,393 measured-tree
held-out forecasts, supporting results, fitted model code, and four research figures.
The root README carries EOTDL YAML metadata and the human-readable description.
The dataset licence is CC BY 4.0; the code and software documentation use MIT.

EOTDL ingests a directory of files with README metadata, not a required ZIP
wrapper. The upload folder is the complete small reproduction asset. The MITECO
downloader obtains the larger inventory inputs separately; its source ZIPs and
extracted databases stay outside the asset.

This package supports reproduction from saved forecasts and includes the inventory
analysis source. It does not contain the raw IFN databases or linked-tree caches.
The raw-data route, including `scripts/download_inventories.py`, is described in REPRODUCING.md. It is not a neural-network model
entry, and the held-out forecast tables should not be described as training data.

## Published version

GitHub release `v1.0.0` (commit `dc6d77f407f4e0273e15fac18acf86791cded1a6`) corresponds to EOTDL dataset version **1**. The uploaded manifest was downloaded from EOTDL and all 115 payload hashes matched that commit. See `provenance/releases/v1.0.0.json` for the permanent mapping. Changes on `main` after this tag are not part of the original upload.

## Prepare and check locally

In the Python 3.11 analysis environment:

```bash
python -m pip install -r requirements-lock.txt
python scripts/reproduce.py
python scripts/prepare_eotdl.py
```

The packaging script creates `dist/iberian-pines-ifn-validation/`. It copies only
specified release directories and file types, preserving the existing relative
imports. Git metadata, environments, caches, local outputs, and raw Access
inventories are excluded. `release-manifest.json` records sizes and SHA-256 hashes.
It refuses to overwrite an existing folder; use `--output` for a new snapshot.

Verify the copy independently:

```bash
python dist/iberian-pines-ifn-validation/scripts/reproduce.py
```

This generates an ignored `outputs/` directory inside the upload folder. Remove
that generated directory or rebuild into a fresh folder before ingesting, so
only the manifest-listed payload is uploaded. Plotting can update the packaged
figures, so always build the final upload copy after verification.

## Upload

Use a separate **Python 3.12+** environment for the EOTDL client. Keeping the
uploader separate preserves the tested analysis dependencies.

```bash
python3.12 -m venv /tmp/ifn-eotdl-uploader
source /tmp/ifn-eotdl-uploader/bin/activate
python -m pip install eotdl
eotdl auth login
eotdl datasets ingest -p ./dist/iberian-pines-ifn-validation
```

Run the last command only when ready to publish. It uploads the folder to EOTDL;
local preparation does not. The name must be available to your account. EOTDL
creates a new version when changed contents are ingested again. Review that
version’s description, licence, and files on its dataset page after ingestion.

The GitHub source URL is configured in the metadata; publish the source repository
before the EOTDL release so this link resolves for readers. No thumbnail URL is
specified until a public image URL is available. The README includes a local figure.

## Scope and attribution

- Code and software documentation: `LICENSE` (MIT).
- Derived CSV/JSON results, reference metrics, and research figures:
  `LICENSE-DATA` (CC BY 4.0), attributed to the authors in README.md.
- Natural Earth geographic background: public domain, as described in DATA.md.
- Original MITECO inventory databases: external inputs, retaining source terms.

## Official references

- [Dataset ingestion and README metadata](https://www.eotdl.com/docs/datasets/ingest)
- [Authentication](https://www.eotdl.com/docs/getting-started/authenticate)
- [Client installation](https://www.eotdl.com/docs/getting-started/install)

Preparation does not assign a DOI, upload a dataset, or publish the GitHub repository.

## GitHub About fields

After authenticating with `gh auth login -h github.com`, run
`python scripts/update_github_about.py` to set the repository description,
website, and topics. The website points to the public EOTDL dataset page.

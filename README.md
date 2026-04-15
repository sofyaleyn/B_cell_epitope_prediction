# B-cell Epitope Prediction Pipeline

Predicts B-cell epitopes for a given antigen using three tools:

| Tool | Mode | Output |
|---|---|---|
| **BepiPred 3.0** | web server (manual download) | linear epitope scores |
| **DiscoTope 3.0** | web server (manual download) | conformational epitope scores |
| **GraphBepi** | local, CPU | conformational epitope scores |

Results are aligned to the full antigen sequence, merged into a per-residue table, and exported as plots, an epitope table, and an executed notebook.

---

## Setup

```bash
uv sync
git submodule update --init --recursive   # pulls graphbepi/
brew install brewsci/bio/dssp             # macOS — required by GraphBepi
```

> **Note:** The `graphbepi/` submodule includes local patches for Python 3.11,
> pytorch-lightning 2.x, and macOS ARM64. Do not replace it with upstream.
> See [CLAUDE.md](CLAUDE.md) for details.

---

## Usage

### If you already have PDB files and FASTA for the target

Register the target in `config/targets.yaml` and run:

```bash
# Run all three stages (predict → combine → report)
uv run python scripts/run_pipeline.py --target ERCC1

# Multiple targets at once
uv run python scripts/run_pipeline.py --target ERCC1 --target E3

# All targets registered in config
uv run python scripts/run_pipeline.py

# Options
uv run python scripts/run_pipeline.py --predictors graphbepi --skip-notebook
```

Or run stages individually:

```bash
uv run python scripts/predict.py --target ERCC1 --predictors all   # stage 1
uv run python scripts/combine.py --target ERCC1                    # stage 2
uv run python scripts/report.py  --target ERCC1                    # stage 3
```

---

## Prerequisites — downloading structures from scratch

If you don't have PDB structures yet, use these two steps first.

### 1. Download all PDB structures for a UniProt entry

Fetches all PDB IDs from UniProt cross-references, downloads FASTA + mmCIF (updated)
from PDBe, and writes `chains.csv` mapping each entry to its target chain.

```bash
# default output: data/<uniprot-id>/raw_pdb/
uv run python scripts/fetch_pdbs.py --uniprot-id P07992

# explicit output directory
uv run python scripts/fetch_pdbs.py --uniprot-id P07992 --out-dir data/my_target/raw_pdb
```

### 2. Extract target chains → PDB files

Reads `chains.csv`, extracts the antigen chain from each `*_updated.cif`,
and writes a cleaned single-chain `.pdb` per entry into `data/<uniprot-id>/pdb/`.

```bash
# default: raw_pdb = data/<uniprot-id>/raw_pdb, out = data/<uniprot-id>/pdb
uv run python scripts/split_pdb_target_chain.py --uniprot-id P07992

# explicit directories
uv run python scripts/split_pdb_target_chain.py --uniprot-id P07992 \
    --raw-pdb-dir data/my_target/raw_pdb --out-dir data/my_target/pdb
```

Then register the target in `config/targets.yaml` and run `run_pipeline.py` as above.

> **TODO:** extract the antigen chain sequence from the downloaded per-entry FASTAs
> (currently requires placing a FASTA manually in `data/<target>/`).

---

## Stage outputs

**Stage 1** (`predict.py`) — saves to `outputs/<target>/`:
- `bepipred_raw.csv` — res_id, residue, bepipred_score
- `discotope_raw.csv` — res_id, residue, discotope_score, average_rsa, n_structures
- `graphbepi_raw.csv` — res_id, chain, residue, score, is_epitope

**Stage 2** (`combine.py`) — saves `outputs/<target>/combined_scores.csv`:

| Column | Description |
|---|---|
| `res_id` | PDB residue sequence number |
| `residue` | Single-letter amino acid |
| `bepipred_score` | BepiPred-3.0 score (0–1) |
| `discotope_score` | DiscoTope-3.0 calibrated score |
| `average_rsa` | Mean relative solvent accessibility |
| `graphbepi_score` | GraphBepi score |
| `is_epitope_bepipred` | score ≥ 0.50 |
| `is_epitope_discotope` | score ≥ 0.90 |
| `is_epitope_graphbepi` | score ≥ 0.1763 |
| `is_epitope_AND` | any predictor flags the residue |

Missing predictor results are filled with 0.0 — combine works with whichever outputs exist.

**Stage 3** (`report.py`) — saves to `outputs/<target>/`:
- `B_cell_epitope_combined.png` — 2- or 3-panel score plot
- `B_cell_epitope_combined_with_interface.png` — with interface overlay (`--include-interface`)
- `epitope_table.csv` — residues where `is_epitope_AND == True`
- `summary_report.html` — HTML summary
- `analysis.ipynb` — executed notebook with inline plots

```bash
uv run python scripts/report.py --target ERCC1 --skip-notebook
uv run python scripts/report.py --target E3 --include-interface   # requires antigen_interface in config
```

---

## Utilities

### Split a PDB/mmCIF by chain

General-purpose chain splitter. Use when you have a single structure file and need
to extract specific chains manually (e.g. before registering a new target).

```bash
# split all chains individually
uv run python scripts/split_pdb.py input.pdb

# keep heavy+light chains (H, L) together; split antigen (A) individually
uv run python scripts/split_pdb.py input.pdb --groups H,L

# write to a specific directory
uv run python scripts/split_pdb.py input.cif --groups H,L --out-dir data/mytarget/pdb/
```

### Diagnose BepiPred web submission

Submits a FASTA to the BepiPred-3.0 web server and prints the raw HTTP response.
Use this when debugging the web-scraping predictor.

```bash
uv run python scripts/debug_bepipred.py --target ERCC1
```

---

## Adding a new target

1. Place data under `data/<target>/` (FASTA, predictor result folders, PDB files)
2. Register in [config/targets.yaml](config/targets.yaml):

```yaml
targets:
  my_target:
    chains: [A]           # antigen chain(s) — used by GraphBepi
    bepipred_dir:  data/my_target/<bepipred_folder>
    discotope_dir: data/my_target/<discotope_folder>
    pdb_dir:       data/my_target/pdb
    antigen_interface: data/my_target/antigen_interface.csv  # optional
```

3. Run stages 1–3 above.

---

## Dependencies

- Python 3.11, managed with `uv`
- `mkdssp` 4.x — `brew install brewsci/bio/dssp` (macOS) or `apt install dssp` (Linux)
- CUDA not required

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

The pipeline runs in three stages.

### Stage 1 — run predictors

```bash
uv run python scripts/predict.py --target ERCC1 --predictors all
uv run python scripts/predict.py --target ERCC1 --predictors graphbepi
uv run python scripts/predict.py --target ERCC1 --predictors bepipred discotope
```

Saves raw results to `outputs/<target>/`:
- `bepipred_raw.csv` — res_id, residue, bepipred_score
- `discotope_raw.csv` — res_id, residue, discotope_score, average_rsa, n_structures
- `graphbepi_raw.csv` — res_id, chain, residue, score, is_epitope

### Stage 2 — combine results

```bash
uv run python scripts/combine.py --target ERCC1
```

Merges all available predictor outputs onto the full antigen sequence.
Saves `outputs/<target>/combined_scores.csv` with columns:

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

### Stage 3 — generate report

```bash
uv run python scripts/report.py --target ERCC1
uv run python scripts/report.py --target ERCC1 --skip-notebook
```

Saves to `outputs/<target>/`:
- `B_cell_epitope_combined.png` — 2- or 3-panel score plot
- `bepipred_profile.png`, `discotope_profile.png`, `graphbepi_profile.png`
- `epitope_table.csv` — residues where `is_epitope_AND == True`
- `summary_report.html` — HTML summary
- `analysis.ipynb` — executed notebook with inline plots

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
```

3. Run stages 1–3 above.

---

## Dependencies

- Python 3.11, managed with `uv`
- `mkdssp` 4.x — `brew install brewsci/bio/dssp` (macOS) or `apt install dssp` (Linux)
- CUDA not required

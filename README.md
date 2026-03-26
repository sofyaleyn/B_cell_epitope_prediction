# B-cell Epitope Prediction Pipeline

Predicts B-cell epitopes for a given antigen using three tools:
- **BepiPred** (web server)
- **DiscoTope** (web server)
- **GraphBepi** (local, CPU)

Results are merged into a consensus score and exported as a plot and HTML report.

## Setup

```bash
uv sync
git submodule update --init --recursive   # pulls graphbepi/
```

## Usage

```bash
# Run full pipeline for a target
uv run python scripts/run_pipeline.py --target her2

# Open analysis notebook
uv run jupyter notebook notebooks/epitope_analysis.ipynb
```

## Adding a new target

Edit [config/targets.yaml](config/targets.yaml) — do not hardcode target names anywhere in `src/`.

## Output schema

All predictors return a DataFrame with columns:

| Column | Type | Description |
|---|---|---|
| `residue_id` | int | Residue sequence number |
| `chain` | str | Chain identifier |
| `score` | float | Epitope probability score |
| `is_epitope` | bool | Binary epitope call |

## Dependencies

- Python 3.11, managed with `uv`
- `dssp` system binary (`brew install brewsci/bio/dssp` or `apt install dssp`)
- CUDA not required

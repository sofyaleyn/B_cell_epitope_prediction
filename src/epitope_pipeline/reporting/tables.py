"""Tabular exports for epitope prediction results."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def export_epitope_table(df: pd.DataFrame, output_path: Path) -> None:
    """Save a CSV containing only the predicted epitope residues (is_epitope_AND).

    Args:
        df:          combined_scores DataFrame from load_and_combine().
        output_path: destination CSV path.
    """
    epitope_col = "is_epitope_AND" if "is_epitope_AND" in df.columns else "is_epitope"
    epitopes = df[df[epitope_col]].copy()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    epitopes.to_csv(output_path, index=False)

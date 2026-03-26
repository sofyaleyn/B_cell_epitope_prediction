"""Epitope score visualisation.

Produces a two-panel per-residue score plot (BepiPred + DiscoTope)
saved to the specified output path.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def plot_epitope_scores(
    df: pd.DataFrame,
    target_name: str,
    output_path: Path,
) -> None:
    """Save a two-panel BepiPred + DiscoTope score plot to output_path."""
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(18, 8))
    gs  = fig.add_gridspec(2, 1, hspace=0.4)

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    bp_epitope = df.get("is_epitope_bepipred", df.get("is_epitope", False))
    dt_epitope = df.get("is_epitope_discotope", df.get("is_epitope", False))

    # BepiPred panel
    ax1.plot(df["res_id"], df["bepipred_score"], color="mediumpurple", linewidth=2, label="Score")
    ax1.axhline(y=0.5, color="gray", linestyle="--", linewidth=1, label="Threshold (0.50)")
    ax1.fill_between(df["res_id"], 0, df["bepipred_score"],
                     where=bp_epitope, alpha=0.8, color="lavender", label="Epitope")
    ax1.set_ylabel("BepiPred Score", fontsize=11, fontweight="bold")
    ax1.set_title(f"Linear B-cell Epitopes (BepiPred 3.0) — {target_name}", fontsize=13, fontweight="bold")
    ax1.set_ylim(-0.05, max(0.65, df["bepipred_score"].max() * 1.1))
    ax1.legend()

    # DiscoTope panel
    ax2.plot(df["res_id"], df["discotope_score"], color="green", linewidth=2, label="Score")
    ax2.axhline(y=0.9, color="gray", linestyle="--", linewidth=1, label="Threshold (0.90)")
    ax2.fill_between(df["res_id"], 0, df["discotope_score"],
                     where=dt_epitope, alpha=0.4, color="mediumseagreen", label="Epitope")
    ax2.set_ylabel("DiscoTope Score", fontsize=11, fontweight="bold")
    ax2.set_title(f"Conformational B-cell Epitopes (DiscoTope 3.0) — {target_name}", fontsize=13, fontweight="bold")
    ax2.set_xlabel("Residue Position", fontsize=11, fontweight="bold")
    ax2.legend()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

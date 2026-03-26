"""Summary report export (HTML).

Renders a Jinja2 template with the combined scores DataFrame and
epitope plot, writing outputs/<target>/summary_report.html.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from jinja2 import Environment, BaseLoader


_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>{{ target_name }} — Epitope Prediction Report</title>
  <style>
    body { font-family: sans-serif; max-width: 960px; margin: 2rem auto; }
    table { border-collapse: collapse; width: 100%; font-size: 0.85rem; }
    th, td { border: 1px solid #ccc; padding: 4px 8px; text-align: right; }
    th { background: #f0f0f0; }
    .epitope { background: #ffe0e0; }
    img { max-width: 100%; margin: 1rem 0; }
  </style>
</head>
<body>
  <h1>{{ target_name }} — B-cell Epitope Prediction Report</h1>
  <img src="epitope_plot.png" alt="Epitope score plot">
  <h2>Predicted epitope residues</h2>
  {{ epitope_table }}
  <h2>All residues</h2>
  {{ full_table }}
</body>
</html>
"""


def export_html(
    df: pd.DataFrame,
    target_name: str,
    output_path: Path,
) -> None:
    """Render an HTML report and write it to output_path."""
    env = Environment(loader=BaseLoader())
    template = env.from_string(_TEMPLATE)

    epitope_col = "is_epitope_AND" if "is_epitope_AND" in df.columns else "is_epitope"
    epitope_df = df[df[epitope_col]]

    html = template.render(
        target_name=target_name,
        epitope_table=_to_html(epitope_df),
        full_table=_to_html(df),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")


def _to_html(df: pd.DataFrame) -> str:
    return df.to_html(index=False, classes="data-table", border=0)

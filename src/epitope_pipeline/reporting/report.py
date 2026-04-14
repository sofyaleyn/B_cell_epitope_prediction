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
    body { font-family: sans-serif; max-width: 1100px; margin: 2rem auto; }
    table { border-collapse: collapse; width: 100%; font-size: 0.85rem; }
    th, td { border: 1px solid #ccc; padding: 4px 8px; text-align: right; }
    th { background: #f0f0f0; }
    .epitope { background: #ffe0e0; }
    img { max-width: 100%; margin: 1rem 0; display: block; }
  </style>
</head>
<body>
  <h1>{{ target_name }} — B-cell Epitope Prediction Report</h1>

  <h2>Combined predictor scores</h2>
  <img src="B_cell_epitope_combined.png" alt="Combined epitope score plot">

  {% if discotope_per_struct_img %}
  <h2>DiscoTope 3.0 — per-structure predictions</h2>
  <img src="{{ discotope_per_struct_img }}" alt="DiscoTope per-structure plot">
  {% endif %}

  {% if graphbepi_per_struct_img %}
  <h2>GraphBepi — per-structure predictions</h2>
  <img src="{{ graphbepi_per_struct_img }}" alt="GraphBepi per-structure plot">
  {% endif %}

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
    discotope_per_struct_img: Path | None = None,
    graphbepi_per_struct_img: Path | None = None,
) -> None:
    """Render an HTML report and write it to output_path.

    Image paths are stored relative to output_path so the HTML is
    self-contained within the same directory.
    """
    env = Environment(loader=BaseLoader())
    template = env.from_string(_TEMPLATE)

    epitope_col = "is_epitope_AND" if "is_epitope_AND" in df.columns else "is_epitope"
    epitope_df = df[df[epitope_col]]

    html = template.render(
        target_name=target_name,
        epitope_table=_to_html(epitope_df),
        full_table=_to_html(df),
        discotope_per_struct_img=(
            discotope_per_struct_img.name if discotope_per_struct_img else None
        ),
        graphbepi_per_struct_img=(
            graphbepi_per_struct_img.name if graphbepi_per_struct_img else None
        ),
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")


def _to_html(df: pd.DataFrame) -> str:
    return df.to_html(index=False, classes="data-table", border=0)

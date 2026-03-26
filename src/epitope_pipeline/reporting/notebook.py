"""Programmatic execution of the analysis notebook template via papermill."""

from __future__ import annotations

from pathlib import Path


def execute_template(
    template: Path,
    output: Path,
    target: str,
    **extra_params: object,
) -> None:
    """Execute the analysis notebook template for a given target.

    Injects ``TARGET`` (and any extra keyword arguments as parameters) via
    papermill, writing the executed notebook to ``output``.

    Args:
        template:     path to notebooks/analysis.ipynb
        output:       destination path for the executed notebook
                      (typically outputs/<target>/analysis.ipynb)
        target:       value injected as the ``TARGET`` parameter
        **extra_params: additional papermill parameters to inject
    """
    try:
        import papermill as pm
    except ImportError:
        raise ImportError("papermill is required: uv add papermill")

    output.parent.mkdir(parents=True, exist_ok=True)

    params = {"TARGET": target, **extra_params}
    pm.execute_notebook(
        str(template),
        str(output),
        parameters=params,
    )

#!/usr/bin/env python3
"""Summarize two independent implementation-vs-R comparison artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROWS = (
    ("Median raw theta relative error", "theta_raw_relative_error_median", "low", "%"),
    ("P90 raw theta relative error", "theta_raw_relative_error_p90", "low", "%"),
    ("Intercept RMSE", "intercept_rmse", "low", "number"),
    ("Residual-variance slope error", "residual_variance_regression_slope", "one", "number"),
    ("Top-3,000 feature overlap", "top_feature_overlap", "high", "%"),
    ("Residual RMSE / R residual SD", "residual_rmse_over_oracle_sd", "low", "%"),
    ("Residual slope error", "residual_regression_slope", "one", "number"),
    ("Fit-gene overlap", "fit_gene_overlap_min_set", "high", "%"),
    ("Transform time", "biolang_elapsed_seconds", "low", "seconds"),
    ("Process wall time", "biolang_process_wall_seconds", "low", "seconds"),
    ("Peak working set", "biolang_peak_working_set_gib", "low", "gib"),
)


def score(value: float, direction: str) -> float:
    if direction == "low":
        return value
    if direction == "high":
        return -value
    return abs(value - 1.0)


def display(value: float, unit: str) -> str:
    if unit == "%":
        return f"{value * 100:.3f}%"
    if unit == "seconds":
        return f"{value:.3f} s"
    if unit == "gib":
        return f"{value:.3f} GiB"
    return f"{value:.6f}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("builtin_json", type=Path)
    parser.add_argument("external_json", type=Path)
    parser.add_argument("output_markdown", type=Path)
    args = parser.parse_args()

    builtin = json.loads(args.builtin_json.read_text(encoding="utf-8"))["metrics"]
    external = json.loads(args.external_json.read_text(encoding="utf-8"))["metrics"]
    lines = [
        "# Three-way SCTransform comparison",
        "",
        "Each accuracy value is measured independently against the same R `sctransform` oracle.",
        "Lower is better for errors and resources; higher is better for overlap.",
        "",
        "| Measurement | BioLang built-in | GPL executable | Closer / lower |",
        "|---|---:|---:|---|",
    ]
    for label, key, direction, unit in ROWS:
        if key not in builtin or key not in external:
            continue
        builtin_value = float(builtin[key])
        external_value = float(external[key])
        winner = "built-in" if score(builtin_value, direction) < score(external_value, direction) else "GPL executable"
        if score(builtin_value, direction) == score(external_value, direction):
            winner = "tie"
        lines.append(
            f"| {label} | {display(builtin_value, unit)} | "
            f"{display(external_value, unit)} | {winner} |"
        )

    lines.extend(
        [
            "",
            "This table does not turn correlation into a parity claim. Review the source JSON",
            "for scale-sensitive gates, slopes, offsets, percentile errors, and exact run metadata.",
            "",
        ]
    )
    args.output_markdown.parent.mkdir(parents=True, exist_ok=True)
    args.output_markdown.write_text("\n".join(lines), encoding="utf-8")
    print(args.output_markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

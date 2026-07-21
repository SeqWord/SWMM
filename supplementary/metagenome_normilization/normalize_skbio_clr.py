#!/usr/bin/env python3
"""
normalize_skbio_clr.py

Normalize a taxa count table using scikit-bio CLR
(closure + multiplicative zero replacement + CLR).

Input:
  --input_folder / -I   (default: input)
  -i / --input          filename inside input_folder

Output:
  --output_folder / -O  (default: output)
  -o / --output         filename inside output_folder

Requires:
  pandas, numpy, scikit-bio
"""

from __future__ import annotations
import argparse
import sys
import os
import numpy as np
import pandas as pd

from skbio.stats.composition import closure, clr, multi_replace


def read_table(path: str, sep: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep, dtype=str)
    if df.shape[1] < 2:
        raise ValueError("Input must contain Taxon + ≥1 sample column.")
    df = df.rename(columns={df.columns[0]: "Taxon"})
    df = df.set_index("Taxon")

    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="raise")

    if (df.values < 0).any():
        raise ValueError("Counts must be non-negative.")
    return df


def main() -> int:
    ap = argparse.ArgumentParser(
        description="CLR normalization of metagenomic count tables using scikit-bio."
    )

    ap.add_argument("-i", "--input", required=True,
                    help="Input filename (inside input_folder).")
    ap.add_argument("-o", "--output", required=True,
                    help="Output filename (inside output_folder).")

    ap.add_argument("-I", "--input_folder", default="input",
                    help="Input folder (default: input).")
    ap.add_argument("-O", "--output_folder", default="output",
                    help="Output folder (default: output).")

    ap.add_argument("--sep", default="\t",
                    help="Field separator (default: tab).")
    ap.add_argument("--delta", type=float, default=None,
                    help="Zero-replacement delta for multi_replace.")
    ap.add_argument("--keep_orientation", action="store_true",
                    help="Keep samples as rows (default: taxa as rows).")

    args = ap.parse_args()

    in_path = os.path.join(args.input_folder, args.input)
    out_path = os.path.join(args.output_folder, args.output)

    os.makedirs(args.output_folder, exist_ok=True)

    df_taxa_by_sample = read_table(in_path, args.sep)

    # scikit-bio expects samples as rows
    df_samples_by_taxa = df_taxa_by_sample.T
    mat = df_samples_by_taxa.to_numpy(dtype=float)

    mat_closed = closure(mat)
    mat_repl = multi_replace(mat_closed, delta=args.delta)
    mat_clr = clr(mat_repl)

    df_clr = pd.DataFrame(
        mat_clr,
        index=df_samples_by_taxa.index,
        columns=df_samples_by_taxa.columns
    )

    if not args.keep_orientation:
        df_clr = df_clr.T

    df_clr.to_csv(out_path, sep="\t", index=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        raise SystemExit(2)

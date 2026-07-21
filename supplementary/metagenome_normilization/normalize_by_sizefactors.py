#!/usr/bin/env python3
"""
normalize_pydeseq2_sizefactors.py (PyDESeq2 v0.5.3-compatible)

Input:  TSV with first column = Taxon, other columns = samples (raw integer counts).
Output: normalized counts (taxa x samples) using DESeq2-style size factors.

PyDESeq2 v0.5.3 note:
- Size factors are stored in dds.obs["size_factors"] (AnnData obs). :contentReference[oaicite:1]{index=1}

Requires:
  pip install pandas numpy pydeseq2==0.5.3
"""

from __future__ import annotations
import argparse
import os
import sys
import numpy as np
import pandas as pd

from pydeseq2.dds import DeseqDataSet


def read_table(path: str, sep: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep, dtype=str)
    if df.shape[1] < 2:
        raise ValueError("Input must contain Taxon + at least one sample column.")

    df = df.rename(columns={df.columns[0]: "Taxon"}).set_index("Taxon")

    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="raise")

    if (df.values < 0).any():
        raise ValueError("Counts must be non-negative.")
    if not np.all(np.equal(np.mod(df.values, 1), 0)):
        raise ValueError("Counts must be integers (or integer-like).")

    return df.astype(int)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="DESeq2 size-factor normalization (PyDESeq2 v0.5.3)."
    )

    ap.add_argument("-i", "--input", required=True,
                    help="Input filename (inside input_folder).")
    ap.add_argument("-o", "--output", required=True,
                    help="Output filename (inside output_folder).")
    ap.add_argument("-I", "--input_folder", default="input",
                    help="Input folder (default: input).")
    ap.add_argument("-O", "--output_folder", default="output",
                    help="Output folder (default: output).")
    ap.add_argument("--size_factors_out", default=None,
                    help="Optional filename for size factors (inside output_folder).")
    ap.add_argument("--sep", default="\t",
                    help="Field separator (default: tab). Use ',' for CSV.")
    ap.add_argument("--method", choices=["ratio", "poscounts", "iterative"],
                    default="poscounts",
                    help="Normalization fit_type (default: poscounts).")

    args = ap.parse_args()

    in_path = os.path.join(args.input_folder, args.input)
    out_path = os.path.join(args.output_folder, args.output)
    os.makedirs(args.output_folder, exist_ok=True)

    df_taxa_by_sample = read_table(in_path, args.sep)

    # PyDESeq2 expects samples x features (samples x taxa)
    counts_df = df_taxa_by_sample.T

    # Minimal metadata; design needs a column that exists
    metadata = pd.DataFrame({"condition": ["A"] * counts_df.shape[0]}, index=counts_df.index)

    # Create dataset (do NOT pass size_factors_fit_type here; use fit_size_factors below)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design="~condition",
        quiet=True,
    )

    # Fit size factors (v0.5.3 stores them in dds.obs["size_factors"])
    dds.fit_size_factors(fit_type=args.method)

    if "size_factors" not in dds.obs:
        raise RuntimeError(
            "PyDESeq2 did not store size_factors in dds.obs. "
            "Confirm you are using pydeseq2==0.5.3 and that fit_size_factors ran."
        )

    sf = dds.obs["size_factors"].to_numpy(dtype=float)

    if sf.shape[0] != counts_df.shape[0]:
        raise RuntimeError("Unexpected size_factors length does not match number of samples.")
    if np.any(sf <= 0) or np.any(~np.isfinite(sf)):
        raise RuntimeError("Invalid size factors (non-positive or non-finite).")

    # Normalize counts
    norm = counts_df.to_numpy(dtype=float) / sf[:, None]
    norm_df = pd.DataFrame(norm, index=counts_df.index, columns=counts_df.columns)

    # Write back as taxa x samples
    norm_df.T.to_csv(out_path, sep="\t", index=True)

    if args.size_factors_out:
        sf_path = os.path.join(args.output_folder, args.size_factors_out)
        pd.DataFrame({"sample": counts_df.index, "size_factor": sf}).to_csv(
            sf_path, sep="\t", index=False
        )

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        raise SystemExit(2)

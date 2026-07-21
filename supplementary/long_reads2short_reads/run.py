#!/usr/bin/env python3
"""
pseudo_hiseq_pe_from_longreads.py

Reads long-read FASTQ inputs from --input_folder and generates pseudo Illumina HiSeq paired-end FASTQ files.

Input files supported (mixed in the same folder):
  - plain:  *.fastq, *.fq
  - gzip:   *.fastq.gz, *.fq.gz
  - zip:    *.fastq.zip, *.fq.zip   (ZIP containing one FASTQ-like member)

For each long-read record, the script tiles along the sequence and emits paired-end reads:
  - R1: forward read_length bases from the fragment start
  - R2: reverse-complement read_length bases from the fragment end
Fragment length is sampled around --insert_size (with --insert_stdev); this yields a
gap (positive) or overlap (negative) between R1 and R2.

Output per input file:
  <basename>_1.fastq   and   <basename>_2.fastq   in --output_folder

Example:
  python pseudo_hiseq_pe_from_longreads.py -i input -o output -l 350 --insert_size 900 --insert_stdev 80 --step 200
"""

from __future__ import annotations

import argparse
import gzip
import io
import os
import random
import sys
import zipfile
from pathlib import Path
from typing import Iterator, Tuple, Optional, TextIO


DNA_COMP = str.maketrans({
    "A": "T", "C": "G", "G": "C", "T": "A",
    "a": "t", "c": "g", "g": "c", "t": "a",
    "N": "N", "n": "n"
})


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


def iter_fastq_records(handle: TextIO) -> Iterator[Tuple[str, str, str]]:
    """
    Yield (name, seq, qual) from a FASTQ text stream.
    name excludes the leading '@' and excludes trailing spaces (keeps first token).
    """
    while True:
        h = handle.readline()
        if not h:
            return
        s = handle.readline()
        plus = handle.readline()
        q = handle.readline()
        if not (s and plus and q):
            return

        h = h.rstrip("\n\r")
        s = s.rstrip("\n\r")
        q = q.rstrip("\n\r")

        if not h.startswith("@"):
            # Forgiving: skip until a header-like line
            continue

        name = h[1:].split()[0]
        yield name, s, q


def sample_fragment_len(read_length: int, insert_size: int, insert_stdev: int) -> int:
    """
    Sample a fragment (template) length around insert_size.
    Allows overlap if frag_len < 2*read_length.
    Ensures frag_len >= read_length.
    """
    if insert_stdev <= 0:
        frag = insert_size
    else:
        frag = int(round(random.gauss(insert_size, insert_stdev)))

    if frag < read_length:
        frag = read_length
    return frag


def write_fastq(out_handle: TextIO, name: str, seq: str, qual: str) -> None:
    out_handle.write(f"@{name}\n{seq}\n+\n{qual}\n")


def normalize_basename(p: Path) -> str:
    """
    Strip typical FASTQ-related extensions to produce a clean base name for outputs.
    Handles .fastq/.fq with optional .gz/.zip.
    """
    name = p.name
    lower = name.lower()

    for ext in (".fastq.gz", ".fq.gz", ".fastq.zip", ".fq.zip", ".fastq", ".fq", ".gz", ".zip"):
        if lower.endswith(ext):
            return name[: -len(ext)]
    return p.stem


def open_fastq_text_from_path(p: Path) -> Tuple[TextIO, Optional[zipfile.ZipFile]]:
    """
    Return a text handle for FASTQ content and (optionally) an owning ZipFile that must be closed by caller.
    For .zip: opens first member ending with .fastq/.fq.
    For .gz: opens gzip text stream.
    For plain: opens text stream.

    Caller must close the returned text handle. If ZipFile is not None, caller must close it too.
    """
    lower = p.name.lower()

    if lower.endswith(".zip"):
        zf = zipfile.ZipFile(p, "r")
        members = [m for m in zf.namelist()
                   if m.lower().endswith((".fastq", ".fq")) and not m.endswith("/")]
        if not members:
            zf.close()
            raise RuntimeError(f"No .fastq/.fq file found inside ZIP: {p}")
        member = members[0]
        raw = zf.open(member, "r")
        txt = io.TextIOWrapper(raw, encoding="utf-8", newline="")
        return txt, zf

    if lower.endswith(".gz"):
        # gzip.open supports text mode directly
        txt = gzip.open(p, "rt", encoding="utf-8", newline="")
        return txt, None

    # plain fastq
    txt = open(p, "r", encoding="utf-8", newline="")
    return txt, None


def process_fastq_stream(
    inp: TextIO,
    out1: TextIO,
    out2: TextIO,
    read_length: int,
    insert_size: int,
    insert_stdev: int,
    step: int,
) -> None:
    for rec_name, seq, qual in iter_fastq_records(inp):
        # If quality is missing/mismatched, fabricate constant qualities
        if len(qual) != len(seq) or len(qual) == 0:
            qual = "I" * len(seq)

        L = len(seq)
        if L < read_length:
            continue

        pos = 0
        frag_idx = 0
        while pos + read_length <= L:
            frag_len = sample_fragment_len(read_length, insert_size, insert_stdev)
            frag_end = pos + frag_len
            if frag_end > L:
                break

            r1_seq = seq[pos:pos + read_length]
            r1_qual = qual[pos:pos + read_length]

            r2_start = frag_end - read_length
            if r2_start < pos:
                break

            r2_seq_fwd = seq[r2_start:frag_end]
            r2_qual_fwd = qual[r2_start:frag_end]
            r2_seq = revcomp(r2_seq_fwd)
            r2_qual = r2_qual_fwd[::-1]

            base = f"{rec_name}|pos={pos}|frag={frag_len}|i={frag_idx}"
            write_fastq(out1, base + "/1", r1_seq, r1_qual)
            write_fastq(out2, base + "/2", r2_seq, r2_qual)

            frag_idx += 1
            pos += step

def iter_with_progress(items, desc: str = "Files"):
    """
    Progress iterator over items. Uses tqdm if available, otherwise prints a simple counter.
    """
    try:
        from tqdm import tqdm  # type: ignore
        return tqdm(items, desc=desc, unit="file")
    except Exception:
        total = len(items)
        def _gen():
            for i, x in enumerate(items, 1):
                print(f"[{desc}] {i}/{total}: {getattr(x, 'name', str(x))}")
                yield x
        return _gen()

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Create pseudo Illumina HiSeq paired-end FASTQ files from long-read FASTQ inputs."
    )
    ap.add_argument("-i", "--input_folder", default="input",
                    help="Folder with input FASTQ / FASTQ.GZ / FASTQ.ZIP files (default: input)")
    ap.add_argument("-o", "--output_folder", default="output",
                    help="Output folder for *_1.fastq and *_2.fastq (default: output)")
    ap.add_argument("-l", "--read_length", type=int, default=350,
                    help="PE read length (default: 350)")
    ap.add_argument("--insert_size", type=int, default=700,
                    help="Mean fragment (template) length in bp (default: 700). "
                         "Gap = insert_size - 2*read_length; negative means overlap.")
    ap.add_argument("--insert_stdev", type=int, default=60,
                    help="Stddev for fragment length sampling (default: 60). Use 0 for fixed insert_size.")
    ap.add_argument("--step", type=int, default=200,
                    help="Step (bp) between successive fragments tiled along a long read (default: 200). "
                         "Smaller step => more overlap between consecutive pairs.")
    ap.add_argument("--seed", type=int, default=None,
                    help="Random seed for reproducibility (default: unset)")
    args = ap.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    in_dir = Path(args.input_folder)
    out_dir = Path(args.output_folder)

    if args.read_length <= 0:
        raise SystemExit("--read_length must be > 0")
    if args.step <= 0:
        raise SystemExit("--step must be > 0")

    if not in_dir.exists() or not in_dir.is_dir():
        raise SystemExit(f"Input folder does not exist or is not a directory: {in_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)

    # Accept plain, gz, zip
    patterns = ["*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz", "*.fastq.zip", "*.fq.zip", "*.gz", "*.zip"]
    seen = set()
    inputs: list[Path] = []
    for pat in patterns:
        for p in in_dir.glob(pat):
            if p.is_file() and p not in seen:
                # Filter to likely FASTQ containers
                lower = p.name.lower()
                if lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fastq.zip", ".fq.zip", ".gz", ".zip")):
                    inputs.append(p)
                    seen.add(p)

    inputs = sorted(inputs)
    if not inputs:
        raise SystemExit(f"No input FASTQ files found in: {in_dir} "
                         f"(expected .fastq/.fq, optionally .gz or .zip)")

    for p in iter_with_progress(inputs, desc="Processing"):
        base = normalize_basename(p)
        out1_path = out_dir / f"{base}_1.fastq"
        out2_path = out_dir / f"{base}_2.fastq"

        try:
            inp, zf = open_fastq_text_from_path(p)
            try:
                with open(out1_path, "w", encoding="utf-8", newline="") as out1, \
                     open(out2_path, "w", encoding="utf-8", newline="") as out2:
                    process_fastq_stream(
                        inp=inp,
                        out1=out1,
                        out2=out2,
                        read_length=args.read_length,
                        insert_size=args.insert_size,
                        insert_stdev=args.insert_stdev,
                        step=args.step,
                    )
            finally:
                inp.close()
                if zf is not None:
                    zf.close()

            print(f"[OK] {p.name} -> {out1_path.name}, {out2_path.name}")

        except Exception as e:
            print(f"[ERROR] Failed on {p}: {e}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

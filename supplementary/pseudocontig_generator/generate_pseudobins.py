#!/usr/bin/env python3

import argparse
import random
from pathlib import Path


def read_fasta(path):
    records = {}
    name = None
    seq = []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)

    if name:
        records[name] = "".join(seq)

    return records


def read_gff(path):
    records = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 9:
                continue

            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue

            records.append({
                "seqid": parts[0],
                "source": parts[1],
                "type": parts[2],
                "start": start,
                "end": end,
                "score": parts[5],
                "strand": parts[6],
                "phase": parts[7],
                "attributes": parts[8],
            })

    return records


def write_wrapped_fasta(fh, name, seq, width=80):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def find_input_genomes(input_dir):
    genomes = []

    for folder in sorted(input_dir.iterdir()):
        if not folder.is_dir():
            continue

        fasta_files = list(folder.glob("*.fa"))
        gff_files = list(folder.glob("*.gff"))

        if not fasta_files or not gff_files:
            print(f"Skipping {folder}: missing .fa or .gff file")
            continue

        genomes.append({
            "name": folder.name,
            "folder": folder,
            "fasta": fasta_files[0],
            "gff": gff_files[0],
            "seqs": read_fasta(fasta_files[0]),
            "gff_records": read_gff(gff_files[0]),
        })

    return genomes


def choose_fragment(genome, min_len, max_len):
    candidates = [
        (seqid, seq)
        for seqid, seq in genome["seqs"].items()
        if len(seq) >= min_len
    ]

    if not candidates:
        return None

    seqid, seq = random.choice(candidates)

    length = random.randint(min_len, min(max_len, len(seq)))
    start = random.randint(1, len(seq) - length + 1)
    end = start + length - 1

    fragment = seq[start - 1:end]

    return seqid, start, end, fragment


def get_overlapping_gff_records(genome, original_seqid, start, end, contig_name):
    output_records = []

    for rec in genome["gff_records"]:
        '''
        if rec["seqid"] != original_seqid:
            continue
        '''

        if rec["end"] < start or rec["start"] > end:
            continue

        new_start = max(rec["start"], start) - start + 1
        new_end = min(rec["end"], end) - start + 1

        output_records.append([
            contig_name,
            rec["source"],
            rec["type"],
            str(new_start),
            str(new_end),
            rec["score"],
            rec["strand"],
            rec["phase"],
            rec["attributes"],
        ])

    return output_records


def main():
    parser = argparse.ArgumentParser(
        prog="make_pseudobins.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""
    Generate pseudo-metagenomic bins from complete bacterial genomes.
    
    The program searches the folder 'input/' for genome subfolders. Each genome
    folder must contain:
    
        *.fa   Genome sequence(s) in FASTA format
        *.gff  GFF file containing methylated nucleotide predictions
    
    Each output bin consists of randomly selected genome fragments ("contigs")
    taken from several species. For every output FASTA file a corresponding GFF
    file is generated containing only methylation records overlapping the selected
    fragments. Coordinates in the GFF file are converted to contig-relative
    coordinates.
    
    Output files are written as
    
        output/<output_folder>/1.fa
        output/<output_folder>/1.gff
        output/<output_folder>/2.fa
        output/<output_folder>/2.gff
        ...
    
    Contig names are written as
    
        GenomeFolderName [start..end]
    
    where start and end are coordinates in the original genome.
    """,
        epilog="""
    Example:
    
    python generate_bins.py \\
        -o test_dataset \\
        -n 200 \\
        -m 3000 \\
        -x 15000 \\
        -d 5 \\
        -v 10 \\
        -s 2 \\
        -t 4 \\
        --seed 17
    """
    )
    
    parser.add_argument("--input_path", default="input")
    parser.add_argument("--output_path", default="output")
    parser.add_argument("-o", "--output_folder", default="", help="Output subfolder name, default=''")
    parser.add_argument("-n", "--bin_number", type=int, default=100, help="Number of bins to generate, default=100")
    parser.add_argument("-m", "--min_contig_length", type=int, default=5000, help="Minimum contig length, default=5,000")
    parser.add_argument("-x", "--max_contig_length", type=int, default=20000, help="Maximum contig length, default=20,000")
    parser.add_argument("-d", "--min_contig_bin", type=int, default=4, help="Minimum contig number per bin, default=4")
    parser.add_argument("-v", "--max_contig_bin", type=int, default=10, help="Maximun contig number per bin, default=10")
    parser.add_argument("-s", "--min_species_number", type=int, default=2, help="Minimum taxon number per bin, default=2")
    parser.add_argument("-t", "--max_species_number", type=int, default=5, help="Maximun taxon number per bin, default=2")
    parser.add_argument("--seed", type=int, default=17)

    args = parser.parse_args()
    random.seed(args.seed)

    input_dir = Path(args.input_path)
    if not input_dir.exists():
        parser.error(f"Input path does not exist: {input_dir}")
    
    if not input_dir.is_dir():
        parser.error(f"Input path is not a directory: {input_dir}")
    
    output_dir = Path(args.output_path) / args.output_folder
    output_dir.mkdir(parents=True, exist_ok=True)

    genomes = find_input_genomes(input_dir)

    if not genomes:
        raise RuntimeError("No valid genome folders found in input/")

    if args.min_species_number > len(genomes):
        raise RuntimeError("min_species_number is larger than the number of available genome folders")

    for bin_id in range(1, args.bin_number + 1):
        contig_number = random.randint(args.min_contig_bin, args.max_contig_bin)

        species_number = random.randint(
            args.min_species_number,
            min(args.max_species_number, contig_number, len(genomes))
        )

        selected_genomes = random.sample(genomes, species_number)

        fasta_out = output_dir / f"{bin_id}.fa"
        gff_out = output_dir / f"{bin_id}.gff"

        with open(fasta_out, "w") as fa_fh, open(gff_out, "w") as gff_fh:
            for contig_id in range(1, contig_number + 1):
                genome = random.choice(selected_genomes)

                fragment_data = choose_fragment(
                    genome,
                    args.min_contig_length,
                    args.max_contig_length
                )

                if fragment_data is None:
                    continue

                original_seqid, start, end, fragment = fragment_data

                contig_name = f"{genome['name']} [{start}..{end}]"

                write_wrapped_fasta(fa_fh, contig_name, fragment)

                gff_records = get_overlapping_gff_records(
                    genome,
                    original_seqid,
                    start,
                    end,
                    contig_name
                )
                
                for rec in gff_records:
                    gff_fh.write("\t".join(rec) + "\n")

        print(f"Created bin {bin_id}: {fasta_out}, {gff_out}")


if __name__ == "__main__":
    main()

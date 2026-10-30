import sys, os
import argparse
from metasorter import MetaSorter as MS

__version__ = version = "1.0.0"
date_of_creation = "21/07/2026"

     
###############################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Program description"
    )

    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"%(prog)s version {__version__} ({date_of_creation})",
        help="Show program's version number and exit"
    )

    parser.add_argument(
        "-i", "--input_folder",
        default="input",
        help="Input folder (default: input)"
    )

    parser.add_argument(
        "-o", "--output_folder",
        default="output",
        help="Output folder (default: output)"
    )

    parser.add_argument(
        "-p", "--project_folder",
        default="",
        help="Project folder (default: current directory)"
    )

    parser.add_argument(
        "--bin_path",
        default="",
        help="Path to executable blast files (default: '', standard location)"
    )

    parser.add_argument(
        "--tmp_path",
        default="",
        help="Path to temporary blast files (default: '', standard location)"
    )

    parser.add_argument(
        "-g", "--gff_file",
        default="",
        help="GFF file with epigenetic predictions (default: '', search for bin-specific gff files)"
    )
    
    parser.add_argument(
        "-m", "--motifs",
        default=("GATC,2,-2; AATT,2,-2; GANTC,2,-2"
            ),
        help="Semicolon-separated motif definitions: 'GATC,2,-2; CCWGG,2,-2"
    )
    
    # Grouping parameters
    parser.add_argument(
        "--grouping_p_cutoff",
        type=float,
        default=0.05,
        help="Cutoff for biased distribution within contigs (default = 0.05)"
    )

    parser.add_argument(
        "--min_motifs",
        type=int,
        default=20,
        help="Minimal number of modified sites per contig (default = 20)"
    )

    parser.add_argument(
        "--min_coverage",
        type=int,
        default=20,
        help="Minimal average coverage at modified sites (default = 20)"
    )

    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="alpha-cutoff for modified motif distribution across motifs (default = 0.05)"
    )

    parser.add_argument(
        "--max_group_difference",
        type=float,
        default=0.1,
        help="Maximum difference (0–1) between methylation fractions within a group (default = 0.1)"
    )

    parser.add_argument(
        "--consensus_threshold",
        type=float,
        default=0.65,
        help="Minimum weighted agreement (0–1) for grouping contigs (default = 0.05)"
    )

    parser.add_argument(
        "--min_motif_weight",
        type=float,
        default=0.1,
        help="Minimum weight (0-1) assigned to an informative motif (default = 0.1)"
    )

    parser.add_argument(
        "--max_motif_weight",
        type=float,
        default=1.0,
        help="Maximum weight (0-1) assigned to an informative motif (default = 1.0)"
    )

    parser.add_argument(
        "--uncertain_weight",
        type=float,
        default=0.1,
        metavar="FLOAT",
        help=(
            "Relative weight (0–1) of evidence involving uncertain contigs. "
            "A value of 0.1 assigns 10%% of the normal motif weight; "
            "0 disables uncertain-contig evidence. (default: %(default)s)"
        )
    )

    # Output parameters
    parser.add_argument(
        "--output_graph_format",
        choices=["SVG", "HTML", "PDF", "EPS", "JPG", "JPEG", "TIF", "TIFF", "PNG", "BMP"],
        default="SVG",
        help="Output graph format (default: SVG)"
    )
    
    parser.add_argument(
        "-u", "--unresolved_contigs",
        choices=["keep", "split"],
        default="keep",
        help="Either keep or split unresolved contigs with insignificant methylation patterns (default: keep)"
    )
        
    parser.add_argument(
        "--filter_chimeric_contigs",
        default="No",
        help="Filter contigs with non-random distribution of modified sites, accepts 'Yes', 'No', or 0 < p-value <= 0.05 (default: No)"
    )
        
    parser.add_argument(
        "-c", "--circular_graph",
        choices=["Y", "y", "N", "n"],
        default="Y",
        help="Produce cicular graph in graphical output (default: Y)"
    )
        
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="Output image DPI (default: 600)"
    )

    parser.add_argument(
        "--mqs",
        type=int,
        default=20,
        help="Minimum methylation quality score (default: 20)"
    )

    parser.add_argument(
        "-l", "--sliding_window_length",
        type=int,
        default=1000,
        help="Sliding window length (default: 2000)"
    )

    parser.add_argument(
        "-w", "--sliding_window_step",
        type=int,
        default=300,
        help="Sliding window step (default: 500)"
    )

    parser.add_argument(
        "--file_name_separator",
        type=str,
        default="",
        help="Separator in file name to select basename as the first part of the name (default: '')"
    )

    parser.add_argument(
        "-s", "--save_graphs",
        choices=["Y", "N", "y", "n"],
        default="Y",
        help="Save graphs (Y/N, default: Y)"
    )

    args = parser.parse_args()
    oMS = MS(version, date_of_creation)
    oMS.execute(args)
    

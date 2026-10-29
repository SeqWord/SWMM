#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse

__version__ = "1.0.2"
date_of_creation = "29 October 2026"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

PROGRAMS = {
    "sort-contigs": os.path.join(SCRIPT_DIR, "lib", "sort_contigs.py"),
    "methylation-pattern": os.path.join(SCRIPT_DIR, "lib", "methylation_pattern.py")
}

# Abbreviated program names
ALIASES = {
    "sc": "sort-contigs",
    "mp": "methylation-pattern"
}

DESCRIPTION = """
Wrapper for sort-contigs and methylation-pattern.

Available programs:

  sort-contigs (sc)
      Sorts binned contigs based on patterns of methylated motifs.

  methylation-pattern (mp)
      Visualizes genome methylation patterns.

Examples:
  python run.py sc --help
  python run.py mp --help
  python run.py sort-contigs --help
  python run.py methylation-pattern --help

All additional arguments are passed directly to the selected program.
"""


def main():

    # No positional argument supplied
    if len(sys.argv) == 1:
        print("Usage:")
        print(f"  {sys.argv[0]} [sort-contigs|sc|methylation-pattern|mp] [options]")
        print(f"  {sys.argv[0]} --help")
        print(f"  {sys.argv[0]} --version")
        return

    first_arg = sys.argv[1]

    # Resolve abbreviated program names
    program_name = ALIASES.get(first_arg, first_arg)

    # Run selected subprogram
    if program_name in PROGRAMS:

        program_path = PROGRAMS[program_name]

        if not os.path.isfile(program_path):
            sys.exit(
                f"ERROR: Program file not found: {program_path}"
            )

        cmd = [
            sys.executable,
            program_path
        ] + sys.argv[2:]

        sys.exit(subprocess.call(cmd))

    # Wrapper-level help/version
    parser = argparse.ArgumentParser(
        prog=os.path.basename(sys.argv[0]),
        description=DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"%(prog)s {__version__} ({date_of_creation})"
    )

    parser.parse_args()


if __name__ == "__main__":
    main()



README.docx
Document
convert the attached file into downloadable properly formatted README.md file for GitHub


README.pdf
PDF
Convert this file to README.md

You're very close. The problem is that you're using the git commit command with the wrong syntax.

1. git add .
This worked correctly:

warning: in the working copy of 'README.md', LF will be replaced by CRLF...
This is only the line-ending warning we discussed before. Your README.md has been staged.

2. The correct git commit syntax
You wrote:

git commit -a .
and

git commit -a "Add README.md"
Both are incorrect.

The -a option means:

Automatically stage modified tracked files before committing.

It does not take a commit message. The commit message is specified with -m.

Since you already ran git add ., you do not need -a at all.

Instead, run:

git commit -m "Add README.md"
3. Push to GitHub
After the commit succeeds, upload it:

git push
4. Ignore the __pycache__ directory
Git is also telling you about:

lib/__pycache__/
You normally do not want to commit Python bytecode. Add the following lines to your .gitignore file:

__pycache__/
*.pyc
If your repository doesn't already have a .gitignore, create one in the repository root with at least:

# Python
__pycache__/
*.py[cod]

# Virtual environments
.venv/
venv/

# OS files
.DS_Store
Thumbs.db
Then run:

git add .gitignore
git commit -m "Ignore Python cache files"
git push
Recommended workflow
For future updates, you can simply use:

git status
git add .
git commit -m "Describe your changes"
git push
This is the standard workflow used by most Git users.


how to see the GitHub SWMM page as a user?

There are several ways to view your GitHub repository exactly as other users see it.

1. Open the repository page (recommended)
Go to:

https://github.com/<username>/<repository>
For example, if your repository belongs to the SeqWord account and is named SWMM, it would be:

https://github.com/SeqWord/SWMM
If the repository is public, this is exactly what any visitor will see.

2. View it while logged out
To ensure you're not seeing anything available only to you:

Open a private/incognito browser window.

Visit the repository URL.

If the repository is public, it will display normally. If it's private, you'll see a 404 Not Found page because GitHub hides private repositories from unauthenticated users.

3. Check how the README is rendered
The repository's README.md is rendered automatically on the repository's main page. After you push a new version:

git add .
git commit -m "Update README"
git push
refresh the page in your browser (press Ctrl+F5 to bypass the cache if necessary).

4. Verify the repository is public
On the repository page:

Click Settings.

Scroll to Danger Zone.

Under Repository visibility, it should say Public.

If it says Private, other users won't be able to see it.

If you tell me the repository URL (or your GitHub username and repository name), I can also suggest improvements to make the page look more polished and professional for visitors.


check README at SeqWord/SWMM. Formatting is not correct. Update the readme.md file

Worked for 1m 36s
I corrected the Markdown structure and formatting. The updated version uses one top-level heading, consistent section levels, properly fenced command examples, and GitHub-compatible tables.


Replace your local README.md with this file, then run:

git add README.md
git commit -m "Fix README formatting"
git push
The original repository page could not be fetched directly through the web tool, so the correction was based on the README content previously generated from your uploaded SWMM documentation.


check the problem:
H:\MyPrograms\GitHub\github_clones\SWMM>git status
On branch main
Your branch is ahead of 'origin/main' by 1 commit.
  (use "git push" to publish your local commits)

Untracked files:
  (use "git add <file>..." to include in what will be committed)
        README(1).md

nothing added to commit but untracked files present (use "git add" to track)

H:\MyPrograms\GitHub\github_clones\SWMM>git add .
warning: in the working copy of 'README(1).md', LF will be replaced by CRLF the next time Git touches it

H:\MyPrograms\GitHub\github_clones\SWMM>git commit -m "Fix README formatting"
[main 5446cb1] Fix README formatting
 1 file changed, 189 insertions(+)
 create mode 100644 README(1).md

H:\MyPrograms\GitHub\github_clones\SWMM>git push
To https://github.com/SeqWord/SWMM.git
 ! [rejected]        main -> main (fetch first)
error: failed to push some refs to 'https://github.com/SeqWord/SWMM.git'
hint: Updates were rejected because the remote contains work that you do not
hint: have locally. This is usually caused by another repository pushing to
hint: the same ref. If you want to integrate the remote changes, use
hint: 'git pull' before pushing again.
hint: See the 'Note about fast-forwards' in 'git push --help' for details.

There are two separate problems.

First, you committed the wrong filename:

README(1).md
GitHub only renders a file named exactly:

README.md
Second, the remote repository has changes that are not in your local clone, so Git rejected the push.

Use these commands.

1. Remove the incorrectly named file from Git tracking
git rm "README(1).md"
2. Replace the real README
Copy the corrected file over the existing file so that the filename is exactly:

README.md
Then stage it:

git add README.md
3. Commit the correction
git commit -m "Replace README with corrected formatting"
4. Pull the remote changes safely
Use rebase so your local commits are placed on top of the latest GitHub version:

git pull --rebase origin main
If Git reports a conflict in README.md, keep the corrected version, then run:

git add README.md
git rebase --continue
5. Push again
git push origin main
Afterwards, check:

git status
It should show:

On branch main
Your branch is up to date with 'origin/main'.

nothing to commit, working tree clean
The key mistake was that Windows saved the downloaded file as README(1).md instead of replacing README.md. GitHub will not use README(1).md as the repository landing-page README.


Library
/
README.md


# SWMM — SeqWord Motif Mapper

**SWMM (SeqWord Motif Mapper)** is a Python program for visualizing and statistically evaluating the distribution of epigenetically modified nucleotides and motifs in bacterial genomes and metagenomic contigs.

SWMM provides two main workflows:

- **`sort-contigs`** — sorts previously binned metagenomic contigs according to their methylation-motif patterns.
- **`methylation-pattern`** — visualizes genome methylation patterns and evaluates the distribution of epigenetically modified bases.

## Contents

- [General usage](#general-usage)
- [Sort contigs](#sort-contigs)
- [Methylation pattern](#methylation-pattern)
- [Contact](#contact)

## General usage

Display the main help message:

```bash
python3 swmm.py --help
```

```text
usage: swmm.py [-h] [-v]
```

| Option | Description |
| --- | --- |
| `-h`, `--help` | Show the help message and exit. |
| `-v`, `--version` | Show the program version and exit. |

## Sort contigs

The `sort-contigs` command sorts binned contigs based on patterns of methylated motifs.

### Help

```bash
python3 swmm.py sort-contigs --help
```

### Usage

```text
usage: sort_contigs.py [-h]
                       [-i INPUT_FOLDER]
                       [-o OUTPUT_FOLDER]
                       [-p PROJECT_FOLDER]
                       [-g GFF_FILE]
                       [-m MOTIFS]
                       [--output_graph_format {SVG,HTML,PDF,EPS,JPG,JPEG,TIF,TIFF,PNG,BMP}]
                       [-u {keep,split}]
                       [--filter_chimeric_contigs FILTER_CHIMERIC_CONTIGS]
                       [-c {Y,y,N,n}]
                       [--dpi DPI]
                       [--mqs MQS]
                       [-l SLIDING_WINDOW_LENGTH]
                       [-w SLIDING_WINDOW_STEP]
                       [--file_name_separator FILE_NAME_SEPARATOR]
                       [-s {Y,N,y,n}]
```

### Arguments

| Option | Description | Default |
| --- | --- | --- |
| `-h`, `--help` | Show the help message and exit. | — |
| `-i INPUT_FOLDER`, `--input_folder INPUT_FOLDER` | Input folder. | `input` |
| `-o OUTPUT_FOLDER`, `--output_folder OUTPUT_FOLDER` | Output folder. | `output` |
| `-p PROJECT_FOLDER`, `--project_folder PROJECT_FOLDER` | Project folder. | Current directory |
| `-g GFF_FILE`, `--gff_file GFF_FILE` | GFF file containing epigenetic predictions. If omitted, the program searches for bin-specific GFF files. | Empty |
| `-m MOTIFS`, `--motifs MOTIFS` | Semicolon-separated motif definitions, for example `GATC,2,-2; CCWGG,2,-2`. | — |
| `--output_graph_format FORMAT` | Output format: `SVG`, `HTML`, `PDF`, `EPS`, `JPG`, `JPEG`, `TIF`, `TIFF`, `PNG`, or `BMP`. | `SVG` |
| `-u {keep,split}`, `--unresolved_contigs {keep,split}` | Keep or split unresolved contigs with statistically insignificant methylation patterns. | `keep` |
| `--filter_chimeric_contigs VALUE` | Filter contigs with a non-random distribution of modified sites. Accepts `Yes`, `No`, or a p-value satisfying `0 < p <= 0.05`. | `No` |
| `-c {Y,y,N,n}`, `--circular_graph {Y,y,N,n}` | Produce a circular graph. | `Y` |
| `--dpi DPI` | Output image resolution in dots per inch. | `600` |
| `--mqs MQS` | Minimum methylation quality score. | `20` |
| `-l LENGTH`, `--sliding_window_length LENGTH` | Sliding-window length. | `2000` |
| `-w STEP`, `--sliding_window_step STEP` | Sliding-window step. | `500` |
| `--file_name_separator SEPARATOR` | Separator used to select the basename as the first part of a filename. | Empty |
| `-s {Y,N,y,n}`, `--save_graphs {Y,N,y,n}` | Save generated graphs. | `Y` |

## Methylation pattern

The `methylation-pattern` command processes GFF files generated by `ipdSummary`, produces graphical representations, and performs statistical analysis of the distribution of epigenetically modified bases.

### Help

```bash
python3 swmm.py methylation-pattern --help
```

### General usage

```bash
python run.py [arguments]
python run.py -h
python run.py --help
python run.py -v
python run.py --version
```

### Example

```bash
python program.py \
    -i input.gff \
    -g genome.gbk \
    -d example \
    -mm Y \
    -w GATC,2,-2 \
    -p 100 \
    -sp Y
```

This example:

1. Reads `input.gff` and `genome.gbk` from `./input/example`.
2. Writes graphical and text output to `./output/example`.
3. Searches for the motif `GATC,2,-2`.
4. Treats the second nucleotide from the left on the direct strand and the second nucleotide from the right on the reverse-complement strand as modified positions.
5. Sets the promoter length to `100` bp.
6. Enables circular-map visualization and motif-distribution statistics.

### General settings

| Option | Description | Default |
| --- | --- | --- |
| `-d FOLDER`, `--project_directory FOLDER` | Project subfolder. | Empty |
| `-i FILE`, `--input_GFF FILE` | Input GFF filename. Required. | Empty |
| `-g FILE`, `--input_GBK FILE` | Input GenBank filename. Required. | Empty |
| `-m FILE`, `--filter_file FILE` | File defining regions to filter. | Empty |
| `-ft NAME`, `--generic_file_name NAME` | Generic output filename. | Empty |
| `-ogf FORMAT`, `--output_graph_format FORMAT` | Output format: `SVG`, `HTML`, `PDF`, `EPS`, `JPG`, `JPEG`, `TIF`, `TIFF`, `PNG`, or `BMP`. | `HTML` |
| `-dpi DPI`, `--dpi DPI` | Raster-image resolution from 100 to 1200 dpi. | `300` |
| `-p LENGTH`, `--promoter_length LENGTH` | Promoter-region length. | `75` |
| `-r NUMBER`, `--maximum_sites NUMBER` | Maximum number of sites to verify. Use `0` to skip checking. | `10000` |
| `-n NUMBER`, `--blast_context_mismatch NUMBER` | Allowed number of context mismatches. | `2` |
| `-z Yes/No`, `--blast_motif_mismatch Yes/No` | Allow motif mismatches. | `Yes` |
| `-u FOLDER`, `--input_folder FOLDER` | Input folder. | `input` |
| `-o FOLDER`, `--output_folder FOLDER` | Output folder. | `output` |
| `-x FOLDER` | Executable folder. | `./lib/bin` |
| `-tmp FOLDER` | Temporary folder. | `./lib/bin/tmp` |

### Circular-map settings

| Option | Description | Default |
| --- | --- | --- |
| `-mm Y/N`, `--circular_map Y/N` | Generate a circular plot. | `Y` |
| `-w MOTIF`, `--cmap_motif MOTIF` | Motif definition, for example `GATC,2,-2`. | — |
| `-s MODE`, `--sites_or_motifs MODE` | Search for sites or motifs. Accepted values include `sites`, `motifs`, `S`, and `M`. | `sites` |
| `-f M/U`, `--modified_or_unmodified M/U` | Display modified or unmodified motifs. | `M` |
| `-wl LENGTH`, `--window_length LENGTH` | Sliding-window length. | `8000` |
| `-ws STEP`, `--window_step STEP` | Sliding-window step. | `2000` |
| `-c SCORE`, `--cmap_score_cutoff SCORE` | Circular-plot score cutoff. | `21` |
| `-cmt TITLE`, `--cmap_graph_title TITLE` | Circular-map graph title. | Empty |

### Dot-plot settings

| Option | Description | Default |
| --- | --- | --- |
| `-dp Y/N`, `--dotplot Y/N` | Generate a dot plot. | `N` |
| `-dpn BASES`, `--nucleotides BASES` | Nucleotides to display: `A`, `C`, `G`, `T`, or an empty value. | `A,C` |
| `-dpm TYPES`, `--methylation_types TYPES` | Methylation types to display, such as `m6A` or `m4C`. | Empty |
| `-dpf MOTIFS`, `--dotplot_motifs MOTIFS` | Motifs to include or exclude, for example `GATC,2,-2; -CRGKGATC,1,6,-2`. | — |
| `-dpc SCORE`, `--dotplot_score_cutoff SCORE` | Dot-plot score cutoff. | `21` |
| `-dpw VALUE`, `--maximum_coverage VALUE` | Maximum coverage on the x-axis. Use `0` for automatic scaling. | `0` |
| `-dps VALUE`, `--maximum_score VALUE` | Maximum score on the y-axis. Use `0` for automatic scaling. | `0` |
| `-dpt TITLE`, `--dp_graph_title TITLE` | Dot-plot graph title. | Empty |

### Statistics settings

| Option | Description | Default |
| --- | --- | --- |
| `-sp Y/N`, `--statplot Y/N` | Generate a statistics panel. | `Y` |
| `-tsk TASKS`, `--tasks TASKS` | Statistical tasks, such as `gc`, `gcs`, and `mge`. | `gc,gcs` |
| `-std MODE`, `--strand MODE` | Exclude the `leading` or `lagging` strand, or use `off` to disable strand exclusion. | `off` |
| `-spt TITLE`, `--sp_graph_title TITLE` | Statistical-panel graph title. | Empty |

## Contact

For further information, refer to the project documentation or contact:

**Oleg Reva**  
University of Pretoria  
Email: `oleg.reva@up.ac.za`
Library
/
README.md


# SWMM — SeqWord Motif Mapper

**SWMM (SeqWord Motif Mapper)** is a Python program for visualizing and statistically evaluating the distribution of epigenetically modified nucleotides and motifs in bacterial genomes and metagenomic contigs.

SWMM provides two main workflows:

- **`sort-contigs`** — sorts previously binned metagenomic contigs according to their methylation-motif patterns.
- **`methylation-pattern`** — visualizes genome methylation patterns and evaluates the distribution of epigenetically modified bases.

## Contents

- [General usage](#general-usage)
- [Sort contigs](#sort-contigs)
- [Methylation pattern](#methylation-pattern)
- [Contact](#contact)

## General usage

Display the main help message:

```bash
python3 swmm.py --help
```

```text
usage: swmm.py [-h] [-v]
```

| Option | Description |
| --- | --- |
| `-h`, `--help` | Show the help message and exit. |
| `-v`, `--version` | Show the program version and exit. |

## Sort contigs

The `sort-contigs` command sorts binned contigs based on patterns of methylated motifs.

### Help

```bash
python3 swmm.py sort-contigs --help
```

### Usage

```text
usage: sort_contigs.py [-h]
                       [-i INPUT_FOLDER]
                       [-o OUTPUT_FOLDER]
                       [-p PROJECT_FOLDER]
                       [-g GFF_FILE]
                       [-m MOTIFS]
                       [--output_graph_format {SVG,HTML,PDF,EPS,JPG,JPEG,TIF,TIFF,PNG,BMP}]
                       [-u {keep,split}]
                       [--filter_chimeric_contigs FILTER_CHIMERIC_CONTIGS]
                       [-c {Y,y,N,n}]
                       [--dpi DPI]
                       [--mqs MQS]
                       [-l SLIDING_WINDOW_LENGTH]
                       [-w SLIDING_WINDOW_STEP]
                       [--file_name_separator FILE_NAME_SEPARATOR]
                       [-s {Y,N,y,n}]
```

### Arguments

| Option | Description | Default |
| --- | --- | --- |
| `-h`, `--help` | Show the help message and exit. | — |
| `-i INPUT_FOLDER`, `--input_folder INPUT_FOLDER` | Input folder. | `input` |
| `-o OUTPUT_FOLDER`, `--output_folder OUTPUT_FOLDER` | Output folder. | `output` |
| `-p PROJECT_FOLDER`, `--project_folder PROJECT_FOLDER` | Project folder. | Current directory |
| `-g GFF_FILE`, `--gff_file GFF_FILE` | GFF file containing epigenetic predictions. If omitted, the program searches for bin-specific GFF files. | Empty |
| `-m MOTIFS`, `--motifs MOTIFS` | Semicolon-separated motif definitions, for example `GATC,2,-2; CCWGG,2,-2`. | — |
| `--output_graph_format FORMAT` | Output format: `SVG`, `HTML`, `PDF`, `EPS`, `JPG`, `JPEG`, `TIF`, `TIFF`, `PNG`, or `BMP`. | `SVG` |
| `-u {keep,split}`, `--unresolved_contigs {keep,split}` | Keep or split unresolved contigs with statistically insignificant methylation patterns. | `keep` |
| `--filter_chimeric_contigs VALUE` | Filter contigs with a non-random distribution of modified sites. Accepts `Yes`, `No`, or a p-value satisfying `0 < p <= 0.05`. | `No` |
| `-c {Y,y,N,n}`, `--circular_graph {Y,y,N,n}` | Produce a circular graph. | `Y` |
| `--dpi DPI` | Output image resolution in dots per inch. | `600` |
| `--mqs MQS` | Minimum methylation quality score. | `20` |
| `-l LENGTH`, `--sliding_window_length LENGTH` | Sliding-window length. | `2000` |
| `-w STEP`, `--sliding_window_step STEP` | Sliding-window step. | `500` |
| `--file_name_separator SEPARATOR` | Separator used to select the basename as the first part of a filename. | Empty |
| `-s {Y,N,y,n}`, `--save_graphs {Y,N,y,n}` | Save generated graphs. | `Y` |

## Methylation pattern

The `methylation-pattern` command processes GFF files generated by `ipdSummary`, produces graphical representations, and performs statistical analysis of the distribution of epigenetically modified bases.

### Help

```bash
python3 swmm.py methylation-pattern --help
```

### General usage

```bash
python run.py [arguments]
python run.py -h
python run.py --help
python run.py -v
python run.py --version
```

### Example

```bash
python program.py \
    -i input.gff \
    -g genome.gbk \
    -d example \
    -mm Y \
    -w GATC,2,-2 \
    -p 100 \
    -sp Y
```

This example:

1. Reads `input.gff` and `genome.gbk` from `./input/example`.
2. Writes graphical and text output to `./output/example`.
3. Searches for the motif `GATC,2,-2`.
4. Treats the second nucleotide from the left on the direct strand and the second nucleotide from the right on the reverse-complement strand as modified positions.
5. Sets the promoter length to `100` bp.
6. Enables circular-map visualization and motif-distribution statistics.

### General settings

| Option | Description | Default |
| --- | --- | --- |
| `-d FOLDER`, `--project_directory FOLDER` | Project subfolder. | Empty |
| `-i FILE`, `--input_GFF FILE` | Input GFF filename. Required. | Empty |
| `-g FILE`, `--input_GBK FILE` | Input GenBank filename. Required. | Empty |
| `-m FILE`, `--filter_file FILE` | File defining regions to filter. | Empty |
| `-ft NAME`, `--generic_file_name NAME` | Generic output filename. | Empty |
| `-ogf FORMAT`, `--output_graph_format FORMAT` | Output format: `SVG`, `HTML`, `PDF`, `EPS`, `JPG`, `JPEG`, `TIF`, `TIFF`, `PNG`, or `BMP`. | `HTML` |
| `-dpi DPI`, `--dpi DPI` | Raster-image resolution from 100 to 1200 dpi. | `300` |
| `-p LENGTH`, `--promoter_length LENGTH` | Promoter-region length. | `75` |
| `-r NUMBER`, `--maximum_sites NUMBER` | Maximum number of sites to verify. Use `0` to skip checking. | `10000` |
| `-n NUMBER`, `--blast_context_mismatch NUMBER` | Allowed number of context mismatches. | `2` |
| `-z Yes/No`, `--blast_motif_mismatch Yes/No` | Allow motif mismatches. | `Yes` |
| `-u FOLDER`, `--input_folder FOLDER` | Input folder. | `input` |
| `-o FOLDER`, `--output_folder FOLDER` | Output folder. | `output` |
| `-x FOLDER` | Executable folder. | `./lib/bin` |
| `-tmp FOLDER` | Temporary folder. | `./lib/bin/tmp` |

### Circular-map settings

| Option | Description | Default |
| --- | --- | --- |
| `-mm Y/N`, `--circular_map Y/N` | Generate a circular plot. | `Y` |
| `-w MOTIF`, `--cmap_motif MOTIF` | Motif definition, for example `GATC,2,-2`. | — |
| `-s MODE`, `--sites_or_motifs MODE` | Search for sites or motifs. Accepted values include `sites`, `motifs`, `S`, and `M`. | `sites` |
| `-f M/U`, `--modified_or_unmodified M/U` | Display modified or unmodified motifs. | `M` |
| `-wl LENGTH`, `--window_length LENGTH` | Sliding-window length. | `8000` |
| `-ws STEP`, `--window_step STEP` | Sliding-window step. | `2000` |
| `-c SCORE`, `--cmap_score_cutoff SCORE` | Circular-plot score cutoff. | `21` |
| `-cmt TITLE`, `--cmap_graph_title TITLE` | Circular-map graph title. | Empty |

### Dot-plot settings

| Option | Description | Default |
| --- | --- | --- |
| `-dp Y/N`, `--dotplot Y/N` | Generate a dot plot. | `N` |
| `-dpn BASES`, `--nucleotides BASES` | Nucleotides to display: `A`, `C`, `G`, `T`, or an empty value. | `A,C` |
| `-dpm TYPES`, `--methylation_types TYPES` | Methylation types to display, such as `m6A` or `m4C`. | Empty |
| `-dpf MOTIFS`, `--dotplot_motifs MOTIFS` | Motifs to include or exclude, for example `GATC,2,-2; -CRGKGATC,1,6,-2`. | — |
| `-dpc SCORE`, `--dotplot_score_cutoff SCORE` | Dot-plot score cutoff. | `21` |
| `-dpw VALUE`, `--maximum_coverage VALUE` | Maximum coverage on the x-axis. Use `0` for automatic scaling. | `0` |
| `-dps VALUE`, `--maximum_score VALUE` | Maximum score on the y-axis. Use `0` for automatic scaling. | `0` |
| `-dpt TITLE`, `--dp_graph_title TITLE` | Dot-plot graph title. | Empty |

### Statistics settings

| Option | Description | Default |
| --- | --- | --- |
| `-sp Y/N`, `--statplot Y/N` | Generate a statistics panel. | `Y` |
| `-tsk TASKS`, `--tasks TASKS` | Statistical tasks, such as `gc`, `gcs`, and `mge`. | `gc,gcs` |
| `-std MODE`, `--strand MODE` | Exclude the `leading` or `lagging` strand, or use `off` to disable strand exclusion. | `off` |
| `-spt TITLE`, `--sp_graph_title TITLE` | Statistical-panel graph title. | Empty |

## Contact

For further information, refer to the project documentation or contact:

**Oleg Reva**  
University of Pretoria  
Email: `oleg.reva@up.ac.za`

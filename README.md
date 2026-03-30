# SimpleMetaPipeline v3.0

[![Tests](https://github.com/J-Cos/SimpleMetaPipeline/actions/workflows/test.yml/badge.svg)](https://github.com/J-Cos/SimpleMetaPipeline/actions/workflows/test.yml)

A reproducible bioinformatic pipeline for metabarcoding data — integrating
Cutadapt, DADA2, LULU, VSEARCH/Swarm, IDTAXA and BLAST into a single
workflow controlled via a validated configuration object.

Published in *Methods in Ecology and Evolution* (2024):
[DOI: 10.1111/2041-210X.14434](https://doi.org/10.1111/2041-210X.14434)

## What's New in v3.0

> **Note**: The v3.0 refactoring was developed with the support of
> [Claude Opus 4.6](https://www.anthropic.com/claude) (Anthropic), which
> assisted with architecture design, code generation, bug fixing, testing,
> and documentation.

- **R package architecture**: `smp_config()` + `smp_run()` replace global variables and sourced scripts
- **Backwards compatible**: loads intermediate outputs from v2 runs automatically
- **Bug fixes**: single-read `inflatedErr` crash, unsorted BLAST top-hit, hardcoded 250 cores, `set.seed(0.1)`
- **Error handling**: all system calls (vsearch, swarm, cutadapt, blastn) now check exit codes
- **Input validation**: `smp_validate()` checks directories, tools, and parameter consistency upfront
- **Test suite**: 26 unit tests + 19 integration tests across 5 datasets
- **Migration tool**: `smp_migrate()` converts old control scripts to v3 format

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/J-Cos/SimpleMetaPipeline.git
cd SimpleMetaPipeline

# 2. Create the conda environment
conda env create -f environment.yml -p ./env

# 3. Activate the environment
conda activate ./env

# 4. Install the R package
Rscript -e 'devtools::install(".")'

# 5. Install LULU (not on CRAN/Bioconductor)
Rscript -e 'devtools::install_github("tobiasgf/lulu")'
```

## Quick Start

```r
library(SimpleMetaPipeline)

config <- smp_config(
    path     = "../BioinformaticPipeline_Env",
    dataname = "Example",
    n_runs   = 2L,
    trunc_len = list(c(150, 150), c(150, 150)),
    trim_left = list(c(0, 0), c(0, 0)),
    training_set = "GTDB_r207-mod_April2022.RData"
)

smp_run(config)
```

Results are saved to `BioinformaticPipeline_Env/Results/`.

## Migrating from v2

If you have an existing control script:

```r
library(SimpleMetaPipeline)
smp_migrate("path/to/your/old/ControlScript.r", output_path = "my_config.R")
```

This reads the old global variables and generates the equivalent `smp_config()` call.
Existing intermediate outputs are loaded automatically — no re-processing needed.

## Directory Structure

```
BioinformaticPipeline_Env/
├── BioinformaticPipeline/    # this repository (R package)
│   ├── R/                    # package source code
│   ├── tests/                # unit + integration tests
│   ├── inst/templates/       # control script template
│   └── env/                  # conda environment
├── FASTQs/
│   └── <dataname>/           # one subdirectory per sequencing run
│       ├── Run1/
│       └── Run2/
├── Data/
│   ├── Classifiers/          # IDTAXA training sets (.RData)
│   └── BlastDBs/             # BLAST databases
├── IntermediateOutputs/      # auto-populated, enables resume
└── Results/                  # final outputs
```

## Configuration

All parameters are passed via `smp_config()`:

| Parameter      | Default            | Description                                       |
| -------------- | ------------------ | ------------------------------------------------- |
| `path`         | —                  | Path to `BioinformaticPipeline_Env`               |
| `dataname`     | —                  | Name matching your FASTQs subdirectory            |
| `read_type`    | `"paired-end"`     | `"paired-end"`, `"single-read"`, or `"premerged"` |
| `n_runs`       | `1`                | Number of sequencing runs                         |
| `trunc_len`    | `list(c(150,150))` | Truncation lengths per run                        |
| `linkage`      | `"complete"`       | `"complete"` (VSEARCH) or `"single"` (Swarm)      |
| `use_cutadapt` | `FALSE`            | Whether to run primer trimming                    |
| `use_idtaxa`   | `TRUE`             | Whether to run IDTAXA classification              |
| `use_blast`    | `FALSE`            | Whether to run BLAST                              |
| `multithread`  | `4`                | Threads (or `TRUE` for auto-detect)               |

See `?smp_config` for the full list.

## Resuming Interrupted Runs

The pipeline caches each step to `IntermediateOutputs/`. Re-run
`smp_run(config)` to resume from the last completed step.

To start fresh: `smp_clear_cache(config)`

## Running Tests

```bash
conda activate ./env
cd BioinformaticPipeline

# Unit tests
Rscript -e 'devtools::test()'

# Integration tests (requires test FASTQs)
Rscript tests/integration/test_all_datasets.R

# End-to-end test with comparison
Rscript tests/integration/test_end_to_end.R
```

## Citation

Williams et al. (2024). A simple bioinformatic pipeline for metabarcoding.
*Methods in Ecology and Evolution*, 15, 2225-2232.
DOI: [10.1111/2041-210X.14434](https://doi.org/10.1111/2041-210X.14434)

## License

GPL-3.0 — see [LICENSE](LICENSE) for details.

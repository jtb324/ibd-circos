# Circos Plot Generator

A command-line tool for visualizing Identity By Descent (IBD) networks using Circos plots. This tool allows users to map pairwise IBD segments shared between individuals and overlay phenotype or variant carrier status for visual analysis.

## Usage

```bash
./circos_maker.R --network <network_file> --ibd <ibd_file> [options]
```

## Arguments

| Flag | Argument | Description | Required |
| :--- | :--- | :--- | :--- |
| `-n`, `--network` | `character` | Path to the input network file (e.g., DRIVE networks). | **Yes** |
| `-i`, `--ibd` | `character` | Path to the pairwise IBD segments file (hap-IBD format). | **Yes** |
| `-p`, `--phenotype` | `character` | Tab-separated file with two columns: `GRID` and `Status`. Overrides network file status. | No |
| `--id` | `character` | Specific Network ID to plot. If provided, only this network is processed. | No |
| `--min-network-size`| `integer` | Minimum number of individuals in a network to generate a plot (Default: 2). | No |
| `--pheno-column` | `character` | Specific column name in the network file to use for phenotype status. | No |
| `-o`, `--output` | `character` | Path for the output image (Default: `test.png`). | No |
| `-h`, `--help` | N/A | Show the help message and exit. | No |

## Requirements

The script requires the following R libraries:
- `ComplexHeatmap`
- `circlize`
- `stringr`
- `optparse`
- `data.table`
- `tidyverse`

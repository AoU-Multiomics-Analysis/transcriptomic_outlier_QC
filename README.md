# Transcriptomic outlier QC

This pipeline is used to identify outlier samples for eQTL calling using WGCNA. Most approaches for identifying outlier samples tend to rely on using PCA, in which case outliers may only be detected by standardizing first one or two principal components. This may not accurately capture all samples that are outliers in the data.

Using WGCNA, we first calculate a correlation matrix (biweight midcorrelation) across all samples after filtering out lowly expressed genes and then performing a rank normal transformation. From this correlation matrix a network is then derived where nodes are samples and the edges describe the weighted correlation between samples. Finally, using the network, outlier samples can be identified by standardizing the connectivity score of each sample. Outlier samples are identified as samples with a connectivity Z score < -3.

---

## Repository Contents

| File | Description |
|------|-------------|
| `transcriptomic_outliers_QC.wdl` | WDL workflow that orchestrates the outlier detection pipeline |
| `identify_sample_outliers.R` | R script that performs gene filtering, normalization, and outlier detection |
| `Dockerfile` | Docker image definition containing all R dependencies |
| `.dockstore.yml` | Dockstore registration configuration for the WDL workflow |

---

## Workflow: `transcriptomic_outliers_QC.wdl`

The WDL workflow wraps the R script into a single task (`call_outliers`) and runs it in a Docker container. It is designed to run on cloud platforms such as Terra or any Cromwell-compatible executor.

### Runtime Configuration

| Parameter | Value |
|-----------|-------|
| Docker image | `evinpadhi/transcriptomic_outlier_qc:latest` |
| CPU | 1 |
| Disk | 500 GB SSD |
| Boot disk | 25 GB |
| Zone | `us-central1-c` |
| Memory | Configurable (see Inputs) |

---

## Script: `identify_sample_outliers.R`

This R script performs the following steps:

1. **Load expression data** – Reads TPM and raw count matrices in RSEM output format (skipping the first 2 header lines).
2. **Filter lowly expressed genes** – Retains only genes where at least 20% of samples have a raw count greater than 6.
3. **Rank-normalize TPMs** – Applies a rank normal transformation (via `RNOmni::RankNorm`) to the TPM values of the filtered genes.
4. **Compute correlation matrix** – Calculates a biweight midcorrelation (bicor) matrix across all samples using WGCNA.
5. **Build network and compute connectivity** – Derives a weighted network where nodes are samples and edge weights are soft-thresholded correlations. Computes the connectivity (sum of edge weights) for each sample.
6. **Identify outliers** – Z-scores the connectivity values; samples with Z score < -3 are flagged as outliers.
7. **Write outputs** – Saves the full Z-score table and a filtered outlier table to TSV files.

### Dependencies

All R package dependencies are installed in the Docker image. The required packages are:

- [`WGCNA`](https://cran.r-project.org/package=WGCNA) – network construction and connectivity calculations
- [`RNOmni`](https://cran.r-project.org/package=RNOmni) – rank normal transformation
- [`tidyverse`](https://www.tidyverse.org/) – data manipulation and I/O
- [`data.table`](https://cran.r-project.org/package=data.table) – fast file reading
- [`optparse`](https://cran.r-project.org/package=optparse) – command-line argument parsing
- [`R.utils`](https://cran.r-project.org/package=R.utils) – general R utilities

---

## Inputs

| Input | Type | Description |
|-------|------|-------------|
| `TPM_path` | `File` | Path to a gene-by-sample TPM expression matrix in RSEM format (tab-separated, first two lines are header rows) |
| `Count_path` | `File` | Path to a gene-by-sample raw count matrix in RSEM format (tab-separated, first two lines are header rows) |
| `OutputPrefix` | `String` | Prefix string used to name the output files (e.g. `my_cohort`) |
| `Memory` | `Int` | Memory to allocate for the job in GB (e.g. `64`) |

### Input File Format

Both the TPM and count files are expected to match the RSEM aggregated output format (e.g., from `rsem-merge-expr-matrices` or equivalent tools). The first two lines are skipped by the script, and the remaining tab-separated table must contain:

- Column `Name` – Ensembl gene IDs (or equivalent gene identifiers)
- Column `Description` – Gene descriptions (dropped during processing)
- One column per sample with expression values

---

## Outputs

| Output | Description |
|--------|-------------|
| `{OutputPrefix}_connectivity_outliers.tsv` | Tab-separated file containing the `ResearchID` and `Z_score` for all samples flagged as outliers (Z score < -3) |
| `{OutputPrefix}_connectivity_scores.tsv` | Tab-separated file containing the `ResearchID` and `Z_score` for **all** samples |

---

## Running the Workflow

### On Terra (or any Cromwell executor)

1. Import `transcriptomic_outliers_QC.wdl` into your Terra workspace (or submit via Cromwell directly).
2. Provide the required inputs in a JSON file:

```json
{
  "transcriptomic_outliers_QC.call_outliers.TPM_path": "gs://your-bucket/path/to/tpm_matrix.tsv",
  "transcriptomic_outliers_QC.call_outliers.Count_path": "gs://your-bucket/path/to/count_matrix.tsv",
  "transcriptomic_outliers_QC.call_outliers.OutputPrefix": "my_cohort",
  "transcriptomic_outliers_QC.call_outliers.Memory": 64
}
```

3. Submit the workflow and retrieve the output TSV files from the execution bucket.

### Running the R Script Directly

The script can also be run outside of WDL using Docker or a local R environment with all dependencies installed:

```bash
Rscript identify_sample_outliers.R \
    --TPM_file /path/to/tpm_matrix.tsv \
    --count_file /path/to/count_matrix.tsv \
    --prefix my_cohort
```

### Using the Docker Image

```bash
docker pull evinpadhi/transcriptomic_outlier_qc:latest

docker run --rm \
    -v /path/to/data:/data \
    evinpadhi/transcriptomic_outlier_qc:latest \
    Rscript /tmp/identify_sample_outliers.R \
        --TPM_file /data/tpm_matrix.tsv \
        --count_file /data/count_matrix.tsv \
        --prefix /data/my_cohort
```

---

## Building the Docker Image

The `Dockerfile` uses [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) to install all R package dependencies into a base conda environment.

```bash
docker build -t transcriptomic_outlier_qc:latest .
```
 

# Pathogenicity Benchmark and Visualization

This repository provides tools for analyzing and visualizing pathogenicity predictions from variant classification
platforms. It compares predictions from Genomize-Seq and Competitor tools against ClinGen's curated truthset, utilizing
metrics like F1, precision, and recall. The results are presented in Sankey diagrams and radar charts.

## Project Overview

Variant classification in Mendelian diseases is a crucial yet challenging task. This project extends the research
presented
in ["Comparative Interpretation of Possible Consequences of Variant Classification Performed on Different Platforms" by İşlek et al.](https://drive.google.com/file/d/1U7TBaK-z9GdFcJzHU92DPZIw4NdJLhL2/view?usp=sharing),
analyzing discrepancies in classifications and providing insights through statistical evaluations and visualizations.

Key goals include:

- Evaluating variant classification accuracy based on ClinGen's five-tier (P, LP, VUS, LB, B) and consolidated
  three-tier models.
- Identifying patterns in overclassification and underclassification across platforms.
- Visualizing comparisons with intuitive plots.

## Key Features

### Comparison Analysis

- `Pathogenicity Classification`: Compares predictions for each variant from Genomize-Seq and Competitor platforms
  against ClinGen's curated truthset.
- `Evidence Code Evaluation`: Evaluates the supporting evidence for each variant classification.
- `Tier Merging`: Supports merging of tiers (e.g., VUS++ and VUS+) for more flexible analysis.
- `Metrics Calculation`: Computes evaluation metrics including F1-score, precision, and recall for each platform's
  predictions.

### Visualization

- `Sankey Diagrams`: Visualizes transitions between the truthset's pathogenicity categories and platform predictions,
  illustrating the degree of classification accuracy or misclassification.
- `Radar Charts`: Visualizes the performance metrics (F1, recall, and precision) for the different platforms in a
  comparative and intuitive manner.

### Data Export

- `Export to TSV`: All comparison results are saved in tab-separated values (TSV) format for easy export and further
  analysis.

## Repository Structure

- `pathogenicity_benchmark.py`: Core logic for comparison and metrics calculation.

- `sankey_diagram.py`: Visualization scripts for Sankey diagrams and radar charts.

- `lib/`: Helper modules for JSON handling, file parsing, and static variables.

- `data/`: Directory for input datasets and generated outputs.

- `Dockerfile and docker-compose.yml`: Setup files for containerized execution.

- `requirements.txt`: List of required Python dependencies.

## Prerequisites

### Using pip

To use this repository with pip, you need Python 3.10 or later and the dependencies listed in `requirements.txt`.
Install the dependencies as follows:

1. Ensure Python 3.10 is installed on your system.
2. Install the required packages using:

    ```bash
    pip install -r requirements.txt
    ```

### Using Conda

To use this repository with Conda:

1. Install [Miniconda or Anaconda](https://docs.anaconda.com/).
2. Create a new Conda environment with Python 3.10:
    ```bash
    conda create -n pathogenicity-benchmark python=3.10
    conda activate pathogenicity-benchmark
    ```
3. Install the required packages using:
    ```bash
    conda install --file requirements.txt
    ```

### Using Docker

To use this repository with Docker:

1. Install [Docker](https://www.docker.com/) on your system.
2. Build the Docker image using the provided `Dockerfile`:

    ```bash
    docker build -t pathogenicity-benchmark .
    ```

## Input Data

1. ClinGen Truthset: JSON file (e.g., `clingen_variant_hg38.json.gz`).

2. Platform Annotations:

- Genomize-Seq Ensembl (e.g., `clingen_annotation_2024feb_ensembl.json.gz`).

- Genomize-Seq RefSeq (e.g., `clingen_annotation_2024feb_refseq.json.gz`).

- Competitor (e.g., `competitor_clingen2024feb_hg38_variant_annotation.tsv.gz`).

## Usage

### Running locally

1. Clone the repository:

    ```bash
    git clone
   
    cd pathogenicity-benchmark
    ```

    2. Install the required dependencies as described in the [Prerequisites](#prerequisites) section.
    3. Run the benchmark script with the required input files:

    ```bash
    python sankey_diagram.py
    ```

### Running with Docker

1. Build the Docker image as described in the [Prerequisites](#using-docker) section.
2. Run the Docker container:

    ```bash
    docker run -it pathogenicity-benchmark
    ```

## Output

The output files are saved in the `data/output` directory:

### Sankey Diagrams

Illustrates transitions between truthset pathogenicity categories and platform predictions.
Example: [sankey_diagram.png](https://drive.google.com/file/d/1rSRh65TstpZkiQBM9LlE-zVP-RS0ROT1/view?usp=sharing).

### Radar Charts

Visualizes F1, precision, and recall metrics for the five-tier or three-tier models.
Example: [radar_chart.pdf](https://drive.google.com/file/d/1gHpmMrpZ8ktnUwKOfeQCe6pxMLTMZ1pZ/view?usp=sharing).

## Statistical Analysis

Metrics such as F1, precision, and recall are calculated for each pathogenicity tier. The results are stored in TSV
files for easy access and interpretation.
Example: [cmpetitor.tsv](https://drive.google.com/file/d/1zPpcYsv3A_1QacJv9mUNM8vdQT4ZiptN/view?usp=sharing).

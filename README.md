# RNAi Target Discovery Pipeline

![Python](https://img.shields.io/badge/python-3.10-blue.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)
![Status](https://img.shields.io/badge/status-active-success.svg)

## 📌 Overview

This bioinformatics pipeline automates the discovery of **double-stranded RNA (dsRNA)** targets for the control of fungal pathogens (specifically designed for *Colletotrichum* spp., but adaptable to other fungi).

The pipeline addresses the three critical challenges in RNAi fungicide design:
1.  **Conservation:** Identifies genomic regions conserved across multiple fungal species/strains to prevent resistance.
2.  **Biosafety (Off-Target Analysis):** Rigorously scans targets against non-target organisms (plants, pollinators, humans) to ensure environmental safety.
3.  **Functional Specificity:** Filters out "housekeeping genes" (e.g., Ribosomes, Actin) to select high-value, pathogen-specific targets.

## ⚙️ Workflow

The pipeline executes the following steps automatically:

1.  **Genome Acquisition:** Downloads reference genomes (CDS) from NCBI for the pathogen and off-target species.
2.  **Ortholog Finding:** Uses **BLASTn** to find genes present in all target fungal strains.
3.  **Multiple Sequence Alignment (MSA):** Aligns orthologs using **Clustal Omega** to identify conserved blocks.
4.  **Conserved Region Scan:** Identifies continuous blocks of DNA (e.g., >60bp) with 95% identity across fungal strains.
5.  **Biosafety Check (In Silico):** BLASTs candidate regions against non-target genomes (e.g., *Citrus sinensis*, *Apis mellifera*, *Homo sapiens*).
    * *Rejection Criteria:* Any match >=21bp with >91% identity is flagged as a risk.
6.  **Functional Annotation:** Retrieves gene descriptions and applies a **Blacklist Filter** to remove generic housekeeping genes.

---

## 🚀 Installation (Conda/Mamba)

This pipeline relies on external binaries (**BLAST+**, **Clustal Omega**, **NCBI Datasets**). The recommended way to install all dependencies correctly is using **Conda** or **Mamba**.

### 1. Clone the Repository
```bash
git clone https://github.com/denilsonfbar/RNAi-Target-Discovery.git
cd RNAi-Target-Discovery
```

### 2. Create the Environment
We provide an `environment.yml` file that installs Python, libraries, and all necessary tools.

```bash
# Create the environment
conda env create -f environment.yml

# Activate the environment
conda activate rnai_pipeline
```

> **Note:** If you don't have Conda installed, we recommend [Miniforge](https://github.com/conda-forge/miniforge).

---

## 📝 Configuration

All pipeline parameters are controlled via the `config.yaml` file. You do not need to edit the Python code.

Open `config.yaml` and adjust the following:

### 1. Define Species
```yaml
# The main pathogen to control
main_species: "Colletotrichum abscissum"

# Related species to ensure the dsRNA works against (Conservation)
conservation_targets:
  - "Colletotrichum gloeosporioides"
  - "Colletotrichum nymphaeae"

# Non-target organisms to protect (Biosafety Check)
off_target_species:
  - "Citrus sinensis"       # Host plant
  - "Apis mellifera"        # Pollinator
  - "Homo sapiens"          # Human safety
```

### 2. Adjust Biosafety Thresholds
These parameters define what constitutes a "dangerous" match based on RNAi biology.
```yaml
biosafety_risk:
  risk_min_length: 21      # Minimum siRNA length (bp)
  risk_min_identity: 91.0  # Max tolerated identity before rejection
```

---

## 🏃 Usage

Once configured, run the pipeline with a single command:

```bash
python main.py
```

The pipeline will print real-time logs to the console, detailing every step of the process.

---

## 📂 Output Structure

All results are saved in the `results/` folder (or whatever you defined in `config.yaml`).

```text
results/
├── 01_genomes/             # Downloaded FASTA files
├── 02_orthologs/           # CSV files with orthologous gene matches
├── 03_alignments/          # Clustal Omega alignments (.fasta)
├── 04_conserved_regions/   # Raw list of conserved sequences
├── 05_biosafety_check/     # BLAST results against non-targets
└── 06_final_report/        # FINAL RESULTS
    ├── FINAL_CANDIDATES_ANNOTATED.csv   # Raw list of molecularly safe IDs
    ├── TARGETS_GOLD_PREMIUM.csv         # <--- START HERE (Approved candidates)
    └── TARGETS_REJECTED_FUNCTIONAL.csv  # Candidates rejected due to function (e.g., Actin)
```

### Interpreting the Results

* **`TARGETS_GOLD_PREMIUM.csv`**: Contains the best candidates. These genes are:
    1.  **Conserved** across all target fungi.
    2.  **Safe** (No significant homology to plants/bees/humans).
    3.  **Specific** (Not generic housekeeping genes like Ribosomes).
* **`TARGETS_REJECTED_FUNCTIONAL.csv`**: Genes that are molecularly safe but biologically risky (e.g., *Actin*, *Tubulin*). Use with caution.

---

## 🛠 Dependencies

* [Python 3.10](https://www.python.org/)
* [Biopython](https://biopython.org/)
* [Pandas](https://pandas.pydata.org/)
* [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi) (via Conda)
* [Clustal Omega](http://www.clustal.org/omega/) (via Conda)
* [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/) (via Conda)

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**Developed for scientific research on RNAi-based crop protection.**
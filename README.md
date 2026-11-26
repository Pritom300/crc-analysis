# crc-analysis

This repository contains the data analysis pipeline and replication package for the research paper on **Colorectal Cancer (CRC)**.

---

##  Project Directory Structure
All necessary dataset components and code structures are organized as follows:
* **`Data/input/`**: 
  * `GPL570-55999(2).txt` — Microarray annotation file.
  * `GSE142987_sample_count_matrix.txt.gz` 
  * `GSE164191_series_matrix.txt.gz`
  * `GSE110223_series_matrix.txt.gz`
  * `GSE41258_series_matrix.txt.gz`
  * `GSE10715_series_matrix.txt.gz`
  * `GSE14333_clinical_data.csv`
  * `GSE14333_series_matrix.txt.gz`  
* **`scripts/`**: Sequential R pipelines named numerically (e.g., `1.Data_Classification.R`, `2...`, `3...`).
* **`output/`**: Execution results including `Stage_Classification_Importance.csv`.
* **`notebooks/`**: `blood_panel.ipynb` (Google Colab deployment ready).

---

## Step-by-Step Execution Workflow

To ensure correct reproducible results, the workflow must be run in the exact order below:

### Phase 1 (R Scripts, RStudio)
You must execute the scripts inside the `scripts/` directory in numerical sequential order (**1 ➔ 2 ➔ 3**):


### Phase 2: (Python, Google Colab)
If you intend to run the Python validation pipeline (`notebooks/blood_panel.ipynb`), follow these configuration steps:

1. **Clone/Load Project to Google Colab**: Load the `blood_panel.ipynb` file from this repository directly inside Google Colab.
2. **Environment Setup**: Ensure you replicate the environment mapping. You must load the following files inside your Colab environment matching the repository tree structure:
   * **`Data/input/`**:
     * `GPL570-55999(2).txt`
     * `GSE142987_sample_count_matrix.txt.gz`
     * `GSE164191_series_matrix.txt.gz`
     * `GSE110223_series_matrix.txt.gz`
     * `GSE41258_series_matrix.txt.gz`
     * `GSE10715_series_matrix.txt.gz`
   * **`output/`**:
     * `Stage_Classification_Importance.csv` *(Generated from the R scripts in Phase 1)*
3. **Run Notebook Cells**: Execute all blocks inside `blood_panel.ipynb` sequentially to reproduce.

---
> **Note:** To avoid any directory errors when saving results, please create the `output` folder along with its subfolder `LIMMA_Results` (`output/LIMMA_Results`) in the root directory before running the code.

## Open Source License
Distributed under the **MIT License**. 


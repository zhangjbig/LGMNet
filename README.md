# LGMNet
## Project Structure

### `R_analysis/` – Macrophage Subtype Identification and Analysis (R)

This folder contains R scripts used for scRNA-seq preprocessing, subtype annotation, functional characterization, and survival analysis of tumor-associated macrophages (TAMs).

#### 1. `seurat/`

Scripts to load scRNA-seq gene expression data, convert them into **Seurat** objects, and identify macrophage populations via clustering and annotation.

#### 2. `metaneighbor.R`

Performs **MetaNeighbor** similarity analysis to assess cross-sample consistency of identified macrophage subtypes.

#### 3. `cluster_sight.R`

Visualizes the distribution of each macrophage cluster and subtype across individual samples, aiding spatial and subtype pattern interpretation.

#### 4. `function.R`

Extracts **shared upregulated DEGs** from macrophage subtypes and performs **GO** and **KEGG** enrichment analysis to infer biological functions.

#### 5. `survival.R`

Assesses the survival relevance of macrophage subtypes using:

* Kaplan–Meier analysis
* Cox proportional hazards regression (uni-/multivariate)
* Outputs subtype scores and labels required by the deep learning model in `DL_model/`.

---

### `DL_model/` – Deep Learning Prediction of Angio-TAMs (Python)

This folder contains Python scripts to build and evaluate deep learning models for **Angio-TAMs** enrichment prediction based on radiomics features extracted from T2-weighted MRI.

#### `LR_select_feature.py`

Uses **logistic regression (LR)** to rank and select the top 20 radiomic features most associated with Angio-TAMs enrichment scores.

#### `LGMNet.py`

Main deep learning model combining:

* Feature selection
* Modified **Inception V2** backbone
* Customizable **MLP classifier** head
  for noninvasive Angio-TAMs prediction.

#### `model_compare.py`

Benchmarks **21 machine learning and deep learning models**, comparing AUC, precision, recall, and F1-score to evaluate model performance against `LGMNet`.


### Dependencies

* R = 4.3.2 (with Seurat=4.4.0, MetaNeighbor=1.22.0, clusterProfiler=4.10.1, survival=3.7.0, etc.)
* Python ≥ 3.10 
* All required Python packages are listed in the `requirements.txt` file.
To install them, run the following command in your terminal:
```bash
pip install -r requirements.txt
```


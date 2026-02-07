# kaggle-melting-point-prediction
Classical ML solution for the Kaggle “Thermophysical Property: Melting Point” competition. End-to-end pipeline with RDKit preprocessing, scaffold-aware CV, XGBoost ensembles, and external data integration. Ranked 479/1177 as a learning project.

## TL;DR

- Classical ML solution for Kaggle melting point prediction
- Rank **479 / 1177**
- End-to-end pipeline with RDKit + XGBoost
- Scaffold-aware GroupKFold validation
- Offline preprocessing + reproducible experiments
- Notebook → script refactor with minimal score drift (~0.3 °C MAE)

**Tech:** Python 3.11, RDKit, scikit-learn, XGBoost  
**Focus:** chemistry-aware validation, clean engineering, reproducibility

# Thermophysical Property Prediction: Melting Point  
**Kaggle Competition – Classical Machine Learning Pipeline**

---

## 1. Introduction

This repository contains a classical machine learning solution for the Kaggle community competition **“Thermophysical Property: Melting Point”**.  
I participated as a beginner for educational purposes and achieved **rank 479 out of 1177** successful competitors.

The task is to predict the **melting temperature (Tm)** of organic molecules based on:
- molecular structure encoded as **SMILES strings**
- **424 molecular group features** provided by the competition

The dataset consists of **3328 compounds**, split by Kaggle into:
- **Training set (≈80%)** with known melting points
- **Test set (≈20%)** without targets

External datasets were permitted under Kaggle’s rules and were used in the final solution.

This repository represents a **fully refactored, script-based version** of the original Kaggle notebook, with explicit preprocessing, modular modeling code, and reproducible experiments.

---

## 2. Project Goals

The main goals of this project were:

- Build an **end-to-end ML pipeline** for molecular property prediction
- Practice **feature engineering with RDKit**
- Apply **scaffold-aware validation** for chemical datasets
- Transition from a **Jupyter notebook** to a maintainable GitHub repository
- Document modeling and engineering decisions transparently

---

## 3. Repository Structure

```text
melting-point-prediction/
│
├── data/
│   ├── 01_raw/              # Original Kaggle data
│   ├── 02_external/         # External dataset (optional)
│   └── 03_processed/        # Fully preprocessed CSV files
│
├── preprocessing/
│   ├── 01_compute_morgan_fp.py
│   ├── 02_compute_descriptors.py
│   ├── 03_merge_external_data.py
│   └── 04_compute_scaffolds_merge.py
│
├── features/
│   ├── feature_config.py    # Feature column definitions
│   ├── scaffold_features.py # Scaffold encoding logic
│
├── modeling/
│   ├── preprocessing.py     # sklearn preprocessing pipeline
│   ├── xgb_models.py        # XGBoost model factory
│   └── ensemble.py          # Seed ensemble + CV logic
│
├── experiments/
│   └── exp_ensemble.py      # Main experiment script
│
├── kaggle_original_setup/
│   └── notebook.ipynb       # Archived Kaggle notebook
│
├── requirements.txt
├── README.md
└── Makefile

```

## 4. Data Preprocessing

All preprocessing is performed offline and orchestrated via a Makefile.

Key steps include:
- Morgan fingerprint computation (RDKit)
- Physicochemical descriptor calculation
- Murcko scaffold generation
- External dataset merging and deduplication
- Consistent preprocessing for train and test sets

Scaffold generation is applied to both train and test datasets, ensuring compatibility with scaffold-based validation and encoding.

## 5. Feature Engineering

The final model uses a combination of:
- Binary Morgan fingerprints
- Continuous RDKit descriptors
- Murcko scaffold identifiers (used for grouping, not as predictors)

All transformations are handled through a unified sklearn preprocessing pipeline.

## 6. Model and Validation Strategy

Model
- XGBoost regressor
- Fixed hyperparameters tuned empirically

Validation
- GroupKFold cross-validation
- Grouping by Murcko scaffold
- Prevents scaffold leakage between folds

Ensembling
- Final predictions are averaged over 5 random seeds
- Identical preprocessing and model configuration per seed

Sample weighting
- External dataset samples are down-weighted to reduce distribution shift

## 7. Final Results

- Public and private Kaggle leaderboard scores for all tried setups are documented in:

documentation/experiment_results.xlsx

- The refactored pipeline reproduces the original Kaggle score within ≈0.3 °C MAE
- This difference is expected due to:
  - seed averaging
  - pipeline refactoring
  - minor numerical differences between environments

## 8. Requirements and Compatibility Notes

- Python ≤ 3.11.9 is required
- RDKit is mandatory for preprocessing
- RDKit depends on NumPy/Pandas wheels unavailable for Python ≥ 3.12
- Kaggle notebooks run in a frozen environment, which differs slightly from local setups
- This repository consolidates all steps into a single virtual environment for simplicity and reproducibility.

## 9. Discussion and Lessons Learned

Several alternative setups were explored, including:
- baseline models without scaffold grouping
- descriptor-only and fingerprint-only models
- different external data weighting schemes

Key takeaways include:
- Scaffold-aware validation is critical in chemistry ML
- Data preprocessing and validation strategy often matter more than model choice
- Refactoring notebooks into scripts improves reproducibility and clarity

## 10. Reproducibility

To reproduce the full pipeline:

```bash
make all
python -m experiments.exp_ensemble
```

This generates a Kaggle-compatible submission file.

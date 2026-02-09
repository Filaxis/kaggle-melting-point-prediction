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

Several modeling setups were explored during the course of this project:

- **Descriptor + NMF baseline**, using 8 RDKit molecular descriptors combined with 50 NMF-reduced molecular group features.
- **Extended descriptor + NMF model**, using 18 molecular descriptors with the same NMF representation. This configuration highlighted the sensitivity of unsupervised feature extraction to preprocessing order and reinforced the need for strict separation between training and validation data.
- **Fingerprint-based model**, replacing NMF-reduced group features with 512-bit Morgan fingerprints (radius = 2) and a reduced set of molecular descriptors.
- **Fingerprint + external data model (final)**, where the Kaggle training data was merged with a deduplicated external dataset to increase chemical diversity, with external samples down-weighted during training.

A detailed comparison of public and private leaderboard MAE values for all tested configurations is provided in `documentation/experiment_results.xlsx`.

When transitioning from NMF-reduced group features to Morgan fingerprints, the molecular descriptor set was intentionally reduced from 18 to 8 features. Many of the removed descriptors (e.g. ExactMolWt, HeavyAtomCount, MolMR, LabuteASA, BertzCT, BalabanJ, and partial charge extrema) are strongly correlated with molecular size or connectivity patterns already captured by the fingerprint representation. Retaining only a compact, chemically interpretable descriptor set helped reduce redundancy and model variance without sacrificing predictive signal.

Key lessons learned include:

- **Feature representation dominates model choice**: replacing NMF-reduced group features with Morgan fingerprints had a larger impact than subsequent hyperparameter tuning.
- **Descriptor redundancy matters**: combining high-dimensional fingerprints with a large number of correlated global descriptors can degrade performance; chemistry-driven pruning improves robustness.
- **Validation strategy is critical in chemistry ML**: scaffold-aware cross-validation provides a more realistic performance estimate, even when it does not immediately translate into public leaderboard gains.
- **External data can improve generalization**, but only when combined with deduplication, scaffold awareness, and appropriate sample weighting.
- **Refactoring notebooks into modular scripts improves reliability**: migrating from a monolithic Kaggle notebook to a structured pipeline made hidden assumptions explicit and improved reproducibility.

Leaderboard selection note

- The code in this repository corresponds to the **best-performing configuration on the private leaderboard**, which used an extended, deduplicated external dataset and scaffold-aware validation.
- Due to Kaggle’s automatic submission selection based on the public leaderboard, a different submission (trained on Kaggle-only data) was selected, reaching the official competition rank **479 / 1177**.
- Had the external-data submission been selected, the final rank would have been approximately **320 / 1177**.

## 10. Reproducibility

To reproduce the full pipeline:

```bash
make all
python -m experiments.exp_ensemble
```

This generates a Kaggle-compatible submission file.

PYTHON := rdkit/Scripts/python.exe

RAW_DIR := data/01_raw
EXT_DIR := data/02_external
PROC_DIR := data/03_processed
PREP := preprocessing

TRAIN := $(RAW_DIR)/melt_train.csv
TEST := $(RAW_DIR)/melt_test.csv
EXT := $(EXT_DIR)/Bradley_Melting_Point_Dataset.csv

all: train test external merge

train:
	$(PYTHON) $(PREP)/01_compute_morgan_fp.py \
		--input $(TRAIN) \
		--output $(PROC_DIR)/train_fp.csv
	$(PYTHON) $(PREP)/02_compute_descriptors.py \
		--input $(PROC_DIR)/train_fp.csv \
		--output $(PROC_DIR)/train_fp_desc.csv
	$(PYTHON) $(PREP)/04_compute_scaffolds.py \
		--input $(PROC_DIR)/train_fp_desc.csv \
		--output $(PROC_DIR)/train_fp_desc_scaffold.csv

test:
	$(PYTHON) $(PREP)/01_compute_morgan_fp.py \
		--input $(TEST) \
		--output $(PROC_DIR)/test_fp.csv
	$(PYTHON) $(PREP)/02_compute_descriptors.py \
		--input $(PROC_DIR)/test_fp.csv \
		--output $(PROC_DIR)/test_fp_desc.csv
	$(PYTHON) $(PREP)/04_compute_scaffolds.py \
		--input $(PROC_DIR)/test_fp_desc.csv \
		--output $(PROC_DIR)/test_fp_desc_scaffold.csv

external:
	$(PYTHON) $(PREP)/01_compute_morgan_fp.py \
		--input $(EXT) \
		--output $(PROC_DIR)/external_fp.csv \
		--add-id \
		--start-id 3330
	$(PYTHON) $(PREP)/02_compute_descriptors.py \
		--input $(PROC_DIR)/external_fp.csv \
		--output $(PROC_DIR)/external_fp_desc.csv
	$(PYTHON) $(PREP)/03_deduplicate_external.py \
		--train $(PROC_DIR)/train_fp_desc.csv \
		--test $(PROC_DIR)/test_fp_desc.csv \
		--external $(PROC_DIR)/external_fp_desc.csv \
		--output $(PROC_DIR)/external_fp_desc_ded.csv
	$(PYTHON) $(PREP)/04_compute_scaffolds.py \
		--input $(PROC_DIR)/external_fp_desc_ded.csv \
		--output $(PROC_DIR)/external_fp_desc_scaffold.csv

merge:
	$(PYTHON) $(PREP)/05_merge_train_external.py \
		--train $(PROC_DIR)/train_fp_desc_scaffold.csv \
		--external $(PROC_DIR)/external_fp_desc_scaffold.csv \
		--output $(PROC_DIR)/train_merged_fp_desc.csv
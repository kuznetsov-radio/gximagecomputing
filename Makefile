.PHONY: parity-roundtrip parity-roundtrip-strict docs-html docs-clean

PYTHON ?= python
SAV_PATH ?= test_data/test.chr.sav
H5_PATH ?= test_data/test.chr.h5
ATOL ?= 0
RTOL ?= 0
SPHINXBUILD ?= sphinx-build
DOCS_DIR ?= docs
DOCS_BUILD_DIR ?= $(DOCS_DIR)/_build

parity-roundtrip:
	HDF5_USE_FILE_LOCKING=FALSE PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src $(PYTHON) tests/RegressionRoundTripSavH5.py \
		--sav-path "$(SAV_PATH)" \
		--h5-path "$(H5_PATH)" \
		--rebuild-h5 \
		--atol $(ATOL) \
		--rtol $(RTOL)

# Alias for CI/local scripts that expect a "strict" regression name.
parity-roundtrip-strict: parity-roundtrip

docs-html:
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src $(SPHINXBUILD) -b html $(DOCS_DIR) $(DOCS_BUILD_DIR)/html

docs-clean:
	rm -rf $(DOCS_BUILD_DIR)

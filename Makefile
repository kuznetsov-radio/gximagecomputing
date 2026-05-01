.PHONY: parity-roundtrip parity-roundtrip-strict docs-html docs-clean

PYTHON ?= python
SAV_PATH ?= $(shell PYTHONPATH=src $(PYTHON) -m gxrender.utils.test_data default-model --suffix .sav 2>/dev/null)
H5_PATH ?= /tmp/gximagecomputing_roundtrip_from_sav.h5
ATOL ?= 0
RTOL ?= 0
SPHINXBUILD ?= sphinx-build
DOCS_DIR ?= docs
DOCS_BUILD_DIR ?= $(DOCS_DIR)/_build

parity-roundtrip:
	@if [ -z "$(SAV_PATH)" ] || [ -z "$(H5_PATH)" ]; then \
		echo "Missing external fixture dataset. Clone ../pyGXrender-test-data next to gximagecomputing and run ../pyGXrender-test-data/scripts/install_dataset.sh, or set GXRENDER_TEST_DATA_ROOT to the extracted raw dataset directory."; \
		exit 1; \
	fi
	HDF5_USE_FILE_LOCKING=FALSE PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src $(PYTHON) scripts/python/RegressionRoundTripSavH5.py \
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

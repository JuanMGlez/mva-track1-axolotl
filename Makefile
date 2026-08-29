VCF ?= data/WGS_EX2312012_HGWCNDSX7.vcf.gz
OUT ?= results

.PHONY: install test run reproduce clean

install:
	python -m pip install -e ".[dev]"

test:
	python -m pytest

run:
	python -m mva_track1 --vcf $(VCF) --out $(OUT)

# What the unattended pipeline would have submitted, with no expert review.
heuristic:
	python -m mva_track1 --vcf $(VCF) --out $(OUT)-heuristic --heuristic-epcr

reproduce: test run

clean:
	rm -rf $(OUT) $(OUT)-heuristic .pytest_cache

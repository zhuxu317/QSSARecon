PYTHON ?= python3
CONFIG ?= config/sim_props.yaml
MECHANISM ?=
MECHANISM_FLAG = $(if $(MECHANISM),--mechanism $(MECHANISM),)

.PHONY: run simulate reconstruct compare

run: simulate reconstruct compare

simulate:
	$(PYTHON) src/simulate_counterflow.py --config $(CONFIG) $(MECHANISM_FLAG)

reconstruct: simulate
	$(PYTHON) src/reconstruct_nh2.py --config $(CONFIG) $(MECHANISM_FLAG)

compare: reconstruct
	$(PYTHON) src/plot_nh2_comparison.py --config $(CONFIG)

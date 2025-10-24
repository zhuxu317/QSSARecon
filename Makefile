# Define the Python interpreter
PYTHON ?= python


.PHONY: all execute 
all:
#----------RUN-1D-overall----------
CASES_DIR ?= cases
# --- NP ---
CASE_NAME ?= NH3_KAUST_DEF_1bar #NH3_PCI_Matrix, NH3_PCI_KAUST_5bar, NH3_NP_exp, NH3_NP_sim, NH3_DNS, NH3_KAUST_NP_1bar
# --- PM ---
PM_CASE_NAME ?= NH3_KAUST_DEF_1bar
PM_CASES_DIR ?= cases

INERT     ?= AR
SORET ?= 1
LOGLEVEL ?= 0
ELEMS ?= 5
PROPS_FILE ?= exp_props.yaml
# PROPS_FILE ?= sim_props.yaml   # default if you keep a single props.yaml
EXTRA_FLAGS ?=



# ---------Single Cases-----
MECHS ?= mechanism/KAUST_32s243r/KAUST_32s243r.yaml 
# MECHS ?= mechanism/H2_pNOX_15_94_TC.yaml
# #--------Multiple Cases-----
# MECHS  =mechanism/Han_35s177r/Han_35s177r.yaml \
# 		mechanism/Mei_39s256r/Mei_39s256r.yaml \
# 		mechanism/Otomo_32s213r/Otomo_32s213r.yaml \
# 		mechanism/Marshall_33s221r/Marshall_33s221r.yaml \
# 		mechanism/NUIG_39s306r/NUIG_39s306r.yaml \
# 		mechanism/KAUST_32s243r/KAUST_32s243r.yaml \
# ELEMS  :=  5 7 8 7 8 7 

# NUIG_39s306r 8 
# Otomo_32s213r 5
# Marshall_33s221r 8
# Han_35s177r.  7
# Mei_39s256r.cti 7


#----------run----------# 

.PHONY: run-1D-counter-flow-NP-CSP
run-1D-counter-flow-NP-CSP:
	@echo "Launching all counterflow cases in parallel…"
	@i=1; \
	for mech in $(MECHS); do \
	  echo "  [$$i] mech=$$mech, props=$(PROPS_FILE)"; \
	  $(PYTHON) src/run_1D_counter_flow_NP_CSP.py \
	    --cases_dir $(CASES_DIR) \
	    --case_name $(CASE_NAME) \
	    --mech_file $$mech \
	    --inert_specie $(INERT) \
	    --loglevel $(LOGLEVEL) \
	    --props_file $(PROPS_FILE) \
	    $(if $(OUT_ROOT),--out_root $(OUT_ROOT),) & \
	  i=$$((i+1)); \
	done; \
	wait; \
	echo "All simulations complete."


.PHONY: run_1D_counter_flow_PPM_CSP
run_1D_counter_flow_PPM_CSP:
	@echo "Launching PREMIXED counterflow cases in parallel…"
	@i=1; \
	for mech in $(MECHS); do \
	  echo "  [$$i] mech=$$mech"; \
        /home/zhuxu21/.conda/envs/viz311/bin/python3 src/run_1D_counter_flow_PPM_CSP.py \
	    --cases_dir $(PM_CASES_DIR) \
	    --case_name $(PM_CASE_NAME) \
	    --mech_file $$mech \
	    --loglevel $(LOGLEVEL) \
	    --props_file $(PROPS_FILE) \
		$(if $(SORET),--soret,) \
	    $(if $(OUT_ROOT),--out_root $(OUT_ROOT),) & \
	  i=$$((i+1)); \
	done; \
	wait; \
	echo "All PREMIXED simulations complete."


.PHONY: run_1D_counter_flow_PM_CSP
run_1D_counter_flow_PM_CSP:
	@echo "Launching PREMIXED counterflow cases in parallel…"
	@i=1; \
	for mech in $(MECHS); do \
	  echo "  [$$i] mech=$$mech, props=$(PROPS_FILE)"; \
	  /home/zhuxu21/.conda/envs/viz311/bin/python3 src/run_1D_counter_flow_PM_CSP.py \
	    --cases_dir $(PM_CASES_DIR) \
	    --case_name $(PM_CASE_NAME) \
	    --mech_file $$mech \
	    --loglevel $(LOGLEVEL) \
	    --props_file $(PROPS_FILE) \
	    $(if $(OUT_ROOT),--out_root $(OUT_ROOT),) & \
	  i=$$((i+1)); \
	done; \
	wait; \
	echo "All PREMIXED simulations complete."



# ---------- POST: QSSA/CSP (generic) ----------
.PHONY: recon-QSSA-CSP
recon-QSSA-CSP:
	@echo "Launching recon (QSSA/CSP) cases in parallel…"
	@i=1; \
	for mech in $(MECHS); do \
	  elem=$$(echo $(ELEMS) | cut -d' ' -f$$i); \
	  echo "  [$$i] mech=$$mech, element_num=$$elem, props=$(PROPS_FILE) EXTRA_SPECIES='$(EXTRA_SPECIES)'"; \
	  /home/zhuxu21/.conda/envs/viz311/bin/python3 src/recon_QSSA_CSP.py \
	    --cases_dir $(CASES_DIR) \
	    --case_name $(CASE_NAME) \
	    --mech_file $$mech \
	    --element_num $$elem \
	    --inert_specie $(INERT) \
	    --loglevel $(LOGLEVEL) \
	    --props_file $(PROPS_FILE) \
	    $(if $(EXTRA_FLAGS),$(EXTRA_FLAGS),) \
	    $(if $(EXTRA_SPECIES),--extra_species $(EXTRA_SPECIES),) \
	    & \
	  i=$$((i+1)); \
	done; \
	wait; \
	echo "All recon runs complete."
	
#----------RUN----------# 
.PHONY:post-PaSR-NH3
post-PaSR-NH3:
	python3 src/post_PaSR_CEQ.py


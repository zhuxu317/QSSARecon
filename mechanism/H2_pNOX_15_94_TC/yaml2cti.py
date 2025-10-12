from pathlib import Path
from cantera.yaml2ck import convert

# Only ask for mechanism & thermo
mech, thermo, _ = convert(
    solution=Path("H2_pNOX_15_94_TC.yaml"),
    phase_name="gas",
    mechanism_path="chem.inp",
    thermo_path="therm.dat",
    transport_path="tran.dat",
    # no transport_path
    sort_elements="alphabetical",
    sort_species="alphabetical",
    overwrite=True
)

print("Wrote:", mech, thermo)
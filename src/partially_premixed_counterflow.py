# Cantera does not provide a direct function for solving partially premixed counterflow flames. 
# It can only handle premixed or diffusion counterflow flames. 
# When air is added to the fuel side, Cantera throws an error when calculating the stoichiometric mixture fraction on the fuel side. 
# To proceed with the calculation, this function must be skipped. 
# This case presents a method for computing partially premixed counterflow flames based on existing functions.
# Developed by Jianyi Jiang on September 2, 2025
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os
# -----------------------
# problem parameter
# H2 example, left:H2/air(phi=1);right:air
# -----------------------
p = ct.one_atm            # Pressure
tin_f = 290            # fuel temperature
tin_o = 290            # oxidizer temperature
width  = 0.02             # width 2 cm
loglevel = 1              # output level
comp_o = 'O2:0.21, N2:0.79'  # air composition
comp_f = "H2:2.0,O2:1.0, N2:3.76" # fuel composition
#choose mass flow rate or velocity inlet
gas = ct.Solution('gri30.yaml')
gas.TPX = tin_o, p, comp_o
rho_o = gas.density    # unit：kg/m³
gas.TPX = tin_f, p, comp_f
rho_f = gas.density    # unit：kg/m³
u_o=0.9 #unit m/s
u_f=0.6 #unit m/s
mdot_o=rho_o * u_o
mdot_f=rho_f * u_f

# ------------------------------------------------
# First step: calculate diffusion counterflow as the initial values
# ------------------------------------------------
# Below is the solution of diffusion counterflow flame 
comp_o_air = 'O2:0.21, N2:0.79'
comp_o_f='H2:1.0'
gas.TP = tin_o, p
# -----------------------
# construct 1D diffusion counterflow
# -----------------------
f = ct.CounterflowDiffusionFlame(gas, width=width)
# Transport model
f.transport_model = 'mixture-averaged'
# solve Criterion
f.set_refine_criteria(ratio=3.0, slope=0.06, curve=0.062)
f.fuel_inlet.mdot = mdot_f
f.fuel_inlet.T = tin_f
f.fuel_inlet.X = comp_o_f
f.oxidizer_inlet.mdot = mdot_o
f.oxidizer_inlet.T = tin_o
f.oxidizer_inlet.X = comp_o_air
# -----------------------
# Robust solve: first disable energy equation, then enable it
# -----------------------
def robust_solve(flame, _log=1):
    try:
        flame.solve(_log, refine_grid=True, auto=True)
        return
    except ct.CanteraError:
        pass  
    flame.energy_enabled = False
    flame.solve(_log, refine_grid=True, auto=False)
    flame.energy_enabled = True
    flame.solve(_log, refine_grid=True, auto=False)

robust_solve(f, loglevel)

# ----------------------------------------------
# Gradual transition: change the fuel side from pure H2 to a mixture of H2:2, O2:1, N2:3.76
# ----------------------------------------------
target_fuel_premix = {'O2':1, 'N2':3.76, 'H2':2}  # target fuel composition
f.set_max_grid_points("flame", 20000)   # mesh limit, change yourself
n_steps=5  # transition step number
alphas = np.linspace(0.2, 1.0, n_steps)  # 0.2, 0.4, ..., 1.0
# ------------------------------
# transition strategy: keep H2 constant, increase air gradually
# --------------------------------
for a in alphas:
    comp_f = f"H2:{target_fuel_premix['H2']}, O2:{a*target_fuel_premix['O2']}, N2:{a*target_fuel_premix['N2']}"
    f.fuel_inlet.X = comp_f
   # Continue using the existing solution as the initial guess, disable auto-initialization (avoid triggering stoich calc)
    f.solve(loglevel, refine_grid=True, auto=False)
    print(f"[step] fuel side now: {comp_f}")
# -------------------------
# Output
# -------------------------
## Simple diagnostics: peak temperature and flame position (coordinate of max temperature)
x = f.grid
T = f.T
Y_H2  = f.Y[f.gas.species_index("H2"), :]
Y_O2  = f.Y[f.gas.species_index("O2"), :]
Y_H2O = f.Y[f.gas.species_index("H2O"), :]
i_max = T.argmax()
print(f"Peak T = {T[i_max]:.2f} K at x = {x[i_max]*1000:.2f} mm")
# -------------------------
# Plot
# -------------------------
plt.figure(figsize=(6,4))
plt.plot(x,T)
plt.xlabel("Position (m)")
plt.ylabel("Temperature")
plt.title("T profile")
plt.grid(True)
plt.tight_layout()
plt.show()
plt.figure(figsize=(6,4))
plt.plot(x, Y_H2,  label="H2")
plt.plot(x, Y_O2,  label="O2")
plt.plot(x, Y_H2O, label="H2O")
plt.xlabel("Position (m)")
plt.ylabel("Mass Fraction")
plt.title("Species Profiles")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
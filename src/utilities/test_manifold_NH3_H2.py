import os, sys, time
import scipy
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from scipy.optimize import linprog
from scipy.linalg import lstsq
from copy import deepcopy
import matplotlib.colors as mcolors

# mech = "mechs/NH3_otomo.cti"
mech = "gri30.yaml"
# mech = "mechs/DRM19_s21r84.yaml"

T0, P0 = 1200, ct.one_atm
phi, fuel, oxid = 1.0, "NH3:0.8, H2:0.2", "O2:0.21, N2:0.79"

def simulate(gas, tot=1e-3):
    r = ct.IdealGasConstPressureReactor(gas)
    sim = ct.ReactorNet([r])
    states = ct.SolutionArray(gas, extra=['t'])
    try:
        while sim.time < tot:
            sim.step()
            states.append(r.thermo.state, t=sim.time)
    except:
        cprint("Error Occured In Cantera Solver!")
        return 0
    return states

gas = ct.Solution(mech)
Ns = gas.n_species
Ne = gas.n_elements

# element matrix B, B[i,j] is the atom number of the jth element for the ith species
B = np.zeros((Ns,Ne))
for i,sp in enumerate(gas.species()):
    for (key,val) in sp.composition.items():
        B[i,gas.element_index(key)] = val

plt.ion()
plt.figure()

# case 1: stoichmetric initial point
# gas.TP = T0, P0
# gas.set_equivalence_ratio(phi, fuel, oxid)
# states = simulate(gas)
# Xeq = states.X[-1,:]
# plt.plot(states.X[:,gas.species_index("CO2")], states.X[:,gas.species_index("H2O")])

# # case 2: highly random start
# NSamples = 10
# X0 = np.exp(np.random.rand(NSamples, Ns)*10-5) # """ DevSkim: ignore DS148264 """
# X0 = (X0.T / np.sum(X0, axis=1)).T
# for i in range(NSamples):
#     gas.TPX = T0, P0, X0[i,:]
#     states = simulate(gas)
#     plt.plot(states.X[:,gas.species_index("CO2")], states.X[:,gas.species_index("H2O")])

# case 3: similar to the paper, start from the equilibrium and perturbate it
gas.TP = T0, P0
gas.set_equivalence_ratio(phi, fuel, oxid)
states = simulate(gas, tot=100)
Xeq = states.X[-1,:]
H0,P0 = gas.HP

c = B.T @ Xeq # element fractions
noise_directions = scipy.linalg.null_space(B.T) # ns*(ns-ne)

idx, idy = gas.species_index("N2"), gas.species_index("H2O")


# idx, idy = gas.species_index("CO2"), gas.species_index("H2O")
# idx, idy = gas.species_index("CO2"), gas.species_index("OH")
plt.scatter(Xeq[idx],Xeq[idy],marker="o",color="k")

NSamples = 25
for i in range(NSamples):
  
    N = deepcopy(Xeq)
    # a, b = np.random.rand(2) * 0.95 #""" DevSkim: ignore DS148264 """
    a, b = (i%5)/5, int(i/5)/5 # grided points
    
    # 减少的N2和H2O 给到NH3和O2和OH
    N_ratio = N[gas.species_index("N2")] * a/100
    H_ratio = N[gas.species_index("H2O")] *  b 
    N[gas.species_index("N2")] -= N_ratio
    N[gas.species_index("N")] += N_ratio*2
    
    N[gas.species_index("H2O")] -= H_ratio
    N[gas.species_index("O2")] += H_ratio/2
    N[gas.species_index("H2")] += H_ratio
    gas.HPX = H0, P0, N/N.sum()

    states = simulate(gas)
    
    # Log-transform the time
    log_time = np.log10(states.t)
    # plt.scatter(states.X[:,idx], states.X[:,idy], marker=".", c=states.t, cmap="jet")
    # Normalize the colormap based on the range of log_time
    norm = mcolors.Normalize(vmin=np.min(log_time), vmax=np.max(log_time))


    idNO = gas.species_index("NO")
    # Scatter plot with log time as color
    scatter = plt.scatter(
        states.X[:, idx],
        states.X[:, idy],
        marker=".",
        c=states.T,
        cmap="viridis",
        s=10  # Marker size in points²
    )
    
    # Scatter id NO
    # scatter = plt.scatter(
    #     states.X[:, idx],
    #     states.X[:, idNO],
    #     marker=".",
    #     c=states.T,
    #     cmap="viridis",
    #     s=10  # Marker size in points²
    # )
        
        
    # scatter = plt.scatter(states.X[:, idy], states.T, marker=".", c=log_time, cmap="viridis", norm=norm)

    # plt.xlim(0, 0.3)
    # plt.ylim(1000, 2800)
    
    plt.draw()
    plt.pause(1e-6)

plt.xlabel("X(H2)")
plt.ylabel("X(H2O)")
# plt.ioff()
plt.savefig("/data/ZhuXu/Cantera/CEMA/figs/test_manifold_new.png")
plt.show(block=True)

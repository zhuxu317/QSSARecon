import os, sys, time
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from scipy.optimize import linprog
from scipy.linalg import lstsq
from functools import lru_cache
import pandas as pd
epsilon = 1e-32
R = 8314.462175
from run_1D_CEQ import run_CEQ_core



class RCCE:
    def __init__(self, mech, constrained_species, use_element_constraints=True):
        """
        Initialize the RCCE class.

        Parameters:
        - mech: Path to the Cantera mechanism file.
        - constrained_species: List of species to constrain.
        - use_element_constraints: If True, uses element constraints (Method 1).
                                    If False, uses Method 2 where Nc = len(constrained_species) + 1.
        """
        self.gas = ct.Solution(mech)
        self.specified_species = constrained_species
        self.use_element_constraints = use_element_constraints

        self.Ns = self.gas.n_species
        self.Ne = self.gas.n_elements
        self.Nr = self.gas.n_reactions
        
        # Select Nc based on the chosen method
        if use_element_constraints:
            self.Nc = self.Ne + len(constrained_species)  # Method 1
        else:
            self.Nc = len(constrained_species) + 1        # Method 2
        
        self.sp_MW = self.gas.molecular_weights
        self.sp_names = np.array(self.gas.species_names)
        self.sp_coeffs = np.zeros((self.gas.n_species, 1 + 7 + 7))
        
        for i, sp in enumerate(self.gas.species()):
            self.sp_coeffs[i, :] = sp.thermo.coeffs
        
        self.nu_f = np.zeros((self.gas.n_reactions, self.gas.n_species))
        self.nu_r = np.zeros((self.gas.n_reactions, self.gas.n_species))

    def set_conditions(self, T, P, X):
        """
        Set the conditions for temperature, pressure, and mole fractions.

        Parameters:
        - T: Temperature.
        - P: Pressure.
        - X: Mole fractions.
        """
        # Set gas state
        self.gas.TPX = T, P, X
        
        # Prepare species dictionary for specified species
        species_dict = {}
        for species in self.specified_species:
            if species in self.gas.species_names:
                species_dict[species] = self.gas[species].X[0]
            else:
                species_dict[species] = None

        # Initialize mole fractions
        mole_fractions = np.zeros(self.gas.n_species)
        for sp in species_dict:
            index = self.gas.species_index(sp)
            mole_fractions[index] = species_dict[sp] if species_dict[sp] is not None else 0.0

        self.X0 = self.gas.X
        N_tot = 1.0
        N0 = self.X0 * N_tot

        # Compute stoichiometric coefficients for forward and reverse reactions
        for i, r in enumerate(self.gas.reactions()):
            for (key, val) in r.reactants.items():
                self.nu_f[i, self.gas.species_index(key)] = val
            for (key, val) in r.products.items():
                self.nu_r[i, self.gas.species_index(key)] = val

        self.el_names = self.gas.element_names
        self.el_numbers = np.zeros((self.gas.n_species, self.Nc))
        c = np.zeros((self.Nc))

        if self.use_element_constraints:
            # Method 1: Use element constraints
            for i, sp in enumerate(self.gas.species()):
                for (key, val) in sp.composition.items():
                    self.el_numbers[i, self.gas.element_index(key)] = val
            B = self.el_numbers
            c = B.T @ N0
            for i, sp in enumerate(self.specified_species):
                self.el_numbers[self.gas.species_index(sp), i + self.Ne] = 1
                c[i + self.Ne] = self.gas[sp].X
        else:
            # Method 2: Nc = len(constrained_species) + 1
            for i, sp in enumerate(self.specified_species):
                self.el_numbers[self.gas.species_index(sp), i + 1] = 1
                c[i + 1] = self.gas[sp].X
            self.el_numbers[:, 0] = 1
            c[0] = self.X0.sum()

        self.gas.set_unnormalized_mole_fractions(mole_fractions)
        self.N0 = self.gas.X
        self.B = self.el_numbers
        
        rank = np.linalg.matrix_rank(self.B)
        full_rank = rank == min(self.B.shape)

        # if full_rank:
            # print("The matrix self.B is full rank.")
        # else:
        #     print("The matrix self.B is not full rank.")
    
        self.c = c

    ##Get properties
    def get_N_ave(self):
        x = self.gas.X
        product = self.Be.T @ x
        N_ave =  self.ce.sum()/product.sum()
        return N_ave
    
    def get_N(self):
        N_ave = self.get_N_ave()
        x = self.gas.X 
        N = x * N_ave
        return N 
    
    def get_H_RT(self, T, Acoeffs):
        T1, T2, T3, T4 = T, T**2, T**3, T**4
        return Acoeffs @ np.array([1, T1/2, T2/3, T3/4, T4/5, 1/T, 0.])

    def get_S0_R(self, T, Acoeffs):
        logT, T1, T2, T3, T4 = np.log(T), T, T**2, T**3, T**4
        return Acoeffs @ np.array([logT, T1, T2/2, T3/3, T4/4, 0., 1])

    def get_U(self,T,P,N): # get total U, mass weighted
        X = N/N.sum()
        sp_coeffs = self.sp_coeffs
        Acoeffs = np.array([sp_coeffs[i,1:8] if T>=sp_coeffs[i,0] else sp_coeffs[i,8:] for i in range(self.Ns)])
        U_k = self.get_H_RT(T, Acoeffs)*R*T - R*T
        return np.dot(U_k, X) / np.dot(self.sp_MW, X)

    def get_H(self,T,P,N): # get total H, mass weighted
        X = N/N.sum()
        sp_coeffs = self.sp_coeffs
        Acoeffs = np.array([sp_coeffs[i,1:8] if T>=sp_coeffs[i,0] else sp_coeffs[i,8:] for i in range(self.Ns)])
        H_RT = self.get_H_RT(T, Acoeffs)
        return np.dot(H_RT*R*T, X) / np.dot(self.sp_MW, X)

    def get_S(self,T,P,N): # get total S, mass weighted
        X = N/N.sum()
        sp_coeffs = self.sp_coeffs
        Acoeffs = np.array([sp_coeffs[i,1:8] if T>=sp_coeffs[i,0] else sp_coeffs[i,8:] for i in range(self.Ns)])
        S_R = self.get_S0_R(T, Acoeffs) - np.log(X+epsilon) - np.log(P/ct.one_atm)
        return np.dot(S_R*R, X) / np.dot(self.sp_MW, X)

    def get_V(self,T,P,N): # get the specific volume
        return np.sum(N)*R*T/P

    def get_T(self,T,P,N):
        return T

    def get_P(self,T,P,N):
        return P

    def get_g0_RT(self,T):
        sp_coeffs = self.sp_coeffs
        Acoeffs = np.array([sp_coeffs[i,1:8] if T>=sp_coeffs[i,0] else sp_coeffs[i,8:] for i in range(self.Ns)])
        Hi0_RT = self.get_H_RT(T, Acoeffs)
        Si0_R  = self.get_S0_R(T, Acoeffs)
        return Hi0_RT - Si0_R

    def get_g_RT(self,T,P): # get G/RT for all species
        return self.get_g0_RT(T) + np.log(P/ct.one_atm)

    def get_mu_RT(self,T, P, N):
        X = N/N.sum()
        return self.get_g0_RT(T) + np.log(P/ct.one_atm) + safe_log(X)

    def get_XP_from_muT(self,mu, T):
        # print(mu, T)
        XP = np.exp(mu - self.get_g0_RT(T)) * ct.one_atm
        P = np.sum(XP)
        X = XP / P
        return X,P

    def get_Gibbs_potential(self,T, P, N):
        X = N/N.sum()
        g_RT = self.get_g_RT(T,P)  # molar specific Gibbs function
        G_tot = np.dot(X, (g_RT + safe_log(X))*R*T) #/ np.dot(sp_MW, X)
        dGdN = (g_RT + safe_log(X) + 1) * R*T #/ np.dot(sp_MW, X)
        return G_tot, dGdN
    
    
    def print_names(self):
        print_idxs = range(self.Ns)
        print(("%-9s "*len(print_idxs))%tuple(self.sp_names[print_idxs]), end=''); print("%-6s"%"Xsum")
        pass

    def print_X(self,X):
        print_idxs = range(self.Ns)
        # print(("%-9.5f "*len(print_idxs))%tuple(X[print_idxs]), end=''); print("%-6.3f"%np.sum(X))
        print(("%-9.3e "*len(print_idxs))%tuple(X[print_idxs]), end=''); print("%-6.3f"%np.sum(X))
        pass


    def ChemEquilibrium(self, T, P, method="CEQ", history = {'step':[], 'X':[]}):
        # Method 3: CEQ, Pope's method Gibbs Function Continuation
        Ns = self.Ns
        if method == "CEQ":
            # global prameters 
            frac_zm = 0.1   # fraction of zm used in inital guess
            ds_inc = 1.4    # factor by which ds is increased after success
            ds_dec = 0.25   # factor by which ds is decreased after failure
            ds_min = 1e-15  # smallest allowed time step
            logy_lim = 120. # upper limit on log(y)  (to prevent overflow
            # err_tol = 1e-6  # convergence tolerance for residual
            err_tol = 1e-9  # convergence tolerance for residual
            
            err_huge = 1e6  # upper limit on error
            # N_tot = np.sum(N0)
            # X0 = self.N0 / N_tot
            X0 = self.N0
            X = X0
            # step 0: construct the constraints matrix B
            #         by default, only element constraints are considered, Nc = Ne
            #         The main species are denoted as determined (d)
            Nc = self.Nc
            
            # TODO: modify the c here, initialize c with more constraints.
            B = self.B # (Ns x Ne)
            # Nb = np.linalg.matrix_rank(B) # number of independent constraints
            # step 1: solve the constraints B.T @ N = c to get the max-min values
            c = self.c
            
            Nx = Ns+1
            Bsum = np.sum(B, axis=0)
            f = np.zeros(Nx)
            f[Ns] = -1.0
            # step 1: solve the constraints B.T @ N = c to get the max-min values
            # the progress described in the FDA report
            AA = np.zeros((Ns,Nx))
            AA[:,:Ns] = -np.diag(np.ones(Ns))
            AA[:,Ns] = 1
            b = np.zeros(Ns)
            Ae = np.zeros((Nc,Nx))
            Ae[:,:Ns] = B.T
            opt = linprog(f,A_ub=AA,b_ub=b,A_eq=Ae,b_eq=c,method='simplex') # default bounds (0,None)
            if opt.status == 0: # success
                Nmin = opt.x[Ns]
                Nmax = opt.x[:Ns]+Nmin
            # get inital guess of zu0
            zm = Nmax
            gu = self.get_g_RT(T,P) # Gibbs function (at s=1) (len(gu)=ns)
            opt = linprog(gu,A_eq=B.T,b_eq=c,method='simplex')
            zg = opt.x
            zg[zg<0] = 0
            # gmin = np.dot(gu,zg)
            zu0 = zg + frac_zm * (zm-zg) # initial guess

            # step 2: set initial conditions for lam0 and g0
            #         lam0 is solved by: log(xu) = -gu + B * lam 
            xu = zu0 / np.sum(zu0) # N_tot
            # rhs = gu + np.log(xu+epsilon)
            # lam0, res, rnk, s = lstsq(B, rhs)
            Binv = np.linalg.pinv(B)
            lam0 = Binv @ (gu + np.log(xu+epsilon))
            gu0 = B @ lam0 - np.log(xu+epsilon)
            
            # step 3: Newton-iteration for g(s) = g(0) + s [g(1) - g(0)]
            s = 0
            send = 1.0
            lam = lam0
            dguds = gu - gu0

            Q = 1.0 + 0.0 # should vary with undetermined species Ns*1
            while (s < 1.0):
                ds = send - s
                if ds<ds_min:
                    sys.exit("Error: integration failure")
                    return
                
                # ceq_rate: get dlam/ds
                gus = gu - (1-s) * dguds
                y = np.exp((-gus + B@lam)/2)
                if np.max(y)>logy_lim:
                    send = s + ds * ds_dec
                    print("[Warning] reducing ds, new ds=", send-s)
                    continue
                # ds = ds/10
                ydgds = y * dguds # Ns
                RR = y * Q
                HH = np.diag(y) @ B # (Ns x Ns) @ (Ns x Nc)
                HHinv = np.linalg.pinv(HH) # (Nc x Ns) #ERROR?
                lamdot_g = HHinv @ ydgds # Nc 
                lamdot_y = HHinv @ y     # Nc
                # print((HH@lamdot_g).shape, (ydgds-HH@lamdot_g).shape, (RR * (HH@lamdot_y)).shape)
                alpha = np.sum(RR * (ydgds-HH@lamdot_g)) / np.sum(RR * (HH@lamdot_y))
                dlamds = lamdot_g + alpha * lamdot_y
                lam1 = lam + ds * dlamds
                
                # ceq_newt: Newton Iteration
                gus = gu - (1-send) * dguds # gu(send)
                lam2 = lam1
                err_best = err_huge
                cnorm = np.linalg.norm(c)

                fail = False
                for i in range(100):
                    y = np.exp((-gus + B@lam2)/2.)
                    if np.max(y)>logy_lim:
                        err = err_huge
                        fail = True
                        break
                    RR = y * Q
                    q  = 1.0 - np.sum( RR * y )
                    HH = np.diag(y) @ B # size(Ns x Nc)
                    HHinv = np.linalg.pinv(HH)
                    v  = y @ HH # (1 x Ns) @ (Ns x Nc), = H'*y
                    vnorm = np.linalg.norm(v)
                    if i>0:
                        err_old = err
                    err = np.linalg.norm((v/vnorm - c/cnorm))  #  residual
                    if err < err_best : #  !  record best solution
                        err_best = err
                        lam_best = lam2
                        y_best   = y
                    if( err + abs(q) < err_tol):
                        break
                    # additional requirements for the Newton iteration
                    # if no more improvement, goback reduce the timestep
                    if( i == 99 or (i>0 and (err-err_old)/(err_old) > 1e-3 )):
                        fail = True
                        break
                    w = vnorm * (c/cnorm - v/vnorm) # size(Nc)
                    lam_v = HHinv @ y  # size(Nc), (Nc x Ns) @ (Ns x 1)
                    lam_w = HHinv @ HHinv.T @ w # size(Nc) (Nc x Ns) @ (Ns x Nc) @ (Nc x 1)
                    P = RR @ HH
                    v_star = np.dot(P, lam_v) #np.sum( P @ lam_v )
                    w_star = np.dot(P, lam_w) # np.sum( P @ lam_w )
                    # print(P.shape, lam_v.shape, lam_w.shape)
                    denom = max( v_star+w_star, 0.5 )
                    d_al  = (q-w_star)/denom # \delta\alpha 
                    d_alp = (q+v_star)/denom # \delta\alpha + 1
                    dlam  = d_alp * lam_w + d_al * lam_v
                    lam2  = lam2 + dlam
                if fail:
                    print("[Warning] reducing ds, new ds=", send-s)
                    send = s + ds * ds_dec
                    continue
                if( err_best < err ):
                    err = err_best
                    lam = lam_best
                    y   = y_best

                lam = lam2
                s    = send  # increment time
                send = min( send + ds * ds_inc,  1.0 )  # increase time step
                X = np.exp(-gus+B@lam)
                if 'out_step' not in history.keys():
                    history['out_step'] = 0
                history['step'].append(s+history['out_step'])
                history['X'].append(X)
                # print_X(X)
                # print("[Info] steping to s =", s)

            zu = y * y
            v  = zu @ B
            zu = zu * cnorm / np.linalg.norm(v)
            N  = zu
        else:
            sys.exit("[ERROR] Method %s is not defined."%method)
        return N, history

    def equilibrium(self, method="TP"):
        # methods: TP; HP, SP, UP; TV; UV, SV
        MaxIter = 100
        History = {}
        # N0 = X0 * N_tot
        N0 = self.N0
        History['step'] = [0]
        History['X'] = [N0]
        gas = self.gas
        T, P = gas.TP
        if method in ["TP", "PT"]:
            nfunc = 0
            N,History = self.ChemEquilibrium(T, P, method="CEQ", history=History)
        elif method in ["HP", "PH"]:
            fixed = "P"
            nfunc = 1
            yfunc = self.get_H
        elif method in ["SP", "PS"]:
            fixed = "P"
            nfunc = 1
            yfunc = self.get_S
        elif method in ["UP", "PU"]:
            fixed = "P"
            nfunc = 1
            yfunc = self.get_U
        elif method in ["TV", "VT"]:
            fixed = "T"
            nfunc = 1
            yfunc = self.get_V # PV=RT, V=RT/P
        elif method in ["UV", "VU"]:
            nfunc = 2
            yfuncs = [self.get_U, self.get_V]
        elif method in ["SV", "VS"]:
            nfunc = 2
            yfuncs = [self.get_S, self.get_V]
        else:
            sys.exit("Not implemented for", method)
        
        if nfunc==1:
            yini = yfunc(T,P,N0)
            print("Initial yval,", yini, " x is", ("T" if fixed=="P" else "P"))
            N, History= self.ChemEquilibrium(T, P, method="CEQ", history=History)
            yold = yfunc(T,P,N)
            dx = (T if fixed=="P" else P) * 0.05
            for i in range(MaxIter):
                History['out_step'] = i+1
                if fixed=="P":
                    T = T + dx
                else:
                    P = P + dx
                N, History= self.ChemEquilibrium(T, P, method="CEQ", history=History)
                ynew = yfunc(T,P,N)
                print("Iter %02d ynew = %15.6f, xnew = %9.4f"%(i,ynew,(T if fixed=="P" else P)))
                if(abs(ynew-yini)<1e-8):
                    break
                dydx = (ynew - yold) / dx
                dx = (yini-ynew) / dydx
                yold = ynew
            print("Finally y,", ynew)
    
        elif nfunc==2:
            def get_y(T,P,N):
                return np.array([yf(T,P,N) for yf in yfuncs])
            yini = get_y(T,P,N0)
            print("Initial yval,", yini)
            N, History = self.ChemEquilibrium(T, P, method="CEQ", history=History)
            yold = get_y(T,P,N)
            jac = np.zeros((nfunc,nfunc))
            res = yold - yini
            for i in range(MaxIter):
                History['out_step'] = i+1

                x = np.array([T,P])
                dx = x*1e-4
                jac[:,0] = (get_y(x[0]+dx[0],x[1],      N) - yold) / dx[0]
                jac[:,1] = (get_y(x[0],      x[1]+dx[1],N) - yold) / dx[1]

                T,P = np.array([T,P]) - 0.5 * np.linalg.inv(jac) @ res
                N, History= self.ChemEquilibrium(T, P, method="CEQ", history=History)
                ynew = get_y(T,P,N)

                print("Iter %02d ynew ="%i, ynew, ", xnew =", np.array([T,P]))
                if np.linalg.norm(ynew/yini-1) < 1e-8:
                    break
                res = ynew - yini
                yold = ynew
            print("Finally y,", ynew)

        X = N / N.sum()
        return (T,P,X), History

def safe_log(x):
    return np.log(np.abs(x)+epsilon)



def test_rcce():
    mech = "src/NH3/NH3_otomo.cti" #"mechs/grimech30_53s325r.cti"
    # mech = "mechs/grimech30_53s325r.cti"
    # T, P, X = 1800, 2*ct.one_atm, "CH4:0.1, O2:0.2, N2:0.7"
    main_species_names = ["NH3", "O2", "H2", "H2O", "N2", "AR"]
    rcce = RCCE(mech,main_species_names, use_element_constraints= True)
    method = "TP" # "TP";  "TV"; "HP", "SP", "UP"; "SV", "UV"
    print("The given condition is:", method)
    print()
    gas = ct.Solution(mech)
    print(f"Number of species: {gas.n_species}")
    print(f"Number of reactions: {gas.n_reactions}")
    print("species_names:", gas.species_names)
    original_data_path = "data/case_NH3_counterflow/N_CF_7.csv"
    df = pd.read_csv(original_data_path)
    row = df.iloc[78]
    Tini = row['T']
    Pini = row['P']
    # input_species_values = [row[name] for name in main_species_names]
    # Xini = dict(zip(main_species_names, input_species_values))
    input_species_values = [row[name] for name in gas.species_names]
    Xini = dict(zip(gas.species_names, input_species_values))
    
    t0 = time.time()
    rcce.set_conditions(Tini, Pini, Xini)
    rcce.gas.equilibrate(method, log_level=0)
    t_opt = time.time() - t0
    Topt, Popt, Xopt = rcce.gas.TPX
    Gopt = rcce.gas.gibbs_mole

    gas.TPX = Tini, Pini, Xini
    Tini, Pini, Xini = gas.TPX
    
    t0 = time.time()
    rcce.set_conditions(Tini, Pini, Xini)
    (Tsim, Psim, Xsim), History = rcce.equilibrium()
    t_sim = time.time() - t0
    Gsim, _ = rcce.get_Gibbs_potential(Tsim, Psim, Xsim)
    
    Tini, Pini, Xini = gas.TPX
    print_idxs = list(set((Xini > 1e-8).nonzero()[0].tolist() + (Xopt > 1e-8).nonzero()[0].tolist()))
    print("Results: Initial, Cantera, Optimal, Newton-")
    rcce.print_X(Xini)
    rcce.print_X(Xopt)
    rcce.print_X(Xsim)
    rcce.print_names()
    
    print("Results: Initial, Cantera, Optimal")
    # print("Initial Gibbs Free Energy = %16.6f, (T,P)=(%.3f, %.1f)" % (Gini, Tini, Pini))
    print("Cantera Gibbs Free Energy = %16.6f, (T,P)=(%.3f, %.1f)" % (Gopt, Topt, Popt))
    print("Optimal Gibbs Free Energy = %16.6f, (T,P)=(%.3f, %.1f)" % (Gsim, Tsim, Psim))
    print("Optimal Maximum error in X = %.4e" % np.max(np.abs(Xopt - Xsim)))
    print("Cantera cost time =", t_opt)
    print("Optimal cost time =", t_sim)
    plt.plot(History['step'], np.array(History['X'])[:, print_idxs])
    if History['out_step'] == 0 and len(History['step']) > 1 and np.max(History['step']) > 1:
        plt.ylim([1e-8, 1])
        plt.xlim([1, np.max(History['step'])])
        plt.xscale('log')
        plt.yscale('log')
    else:
        plt.ylim([1e-12, 1])
        plt.yscale('log')
    plt.xlabel("Step")
    plt.ylabel("Mole Fraction")
    plt.legend(rcce.sp_names[print_idxs], ncol=5)
    plt.savefig("figs/rcce_test.png")
        
    # COMPARE QSSA and RCCE with or without element conservation
    input_species_values = [row[name] for name in main_species_names]
    Xini_CEQ = dict(zip(main_species_names, input_species_values))
    gas = run_CEQ_core(gas, main_species_names, Xini_CEQ, Tini, Pini, dt=1e-1, time_end=1e0)
    X_qssa = gas.X  # Mole fractions from QSSA

    # Calculate relative errors for Xsim and X_qssa with respect to Xini
    epsilon = 1e-12  # Small value to avoid division by zero
    relative_error_sim = np.abs((Xsim - Xini) / (Xini + epsilon) * 100)
    relative_error_qssa = np.abs((X_qssa - Xini) / (Xini + epsilon) * 100)

    # Create the relative error plot
    x = np.arange(len(gas.species_names))  # X positions for the species

    plt.figure(figsize=(14, 6))
    plt.plot(x, relative_error_sim, label='Relative Error (Xsim vs Xini)', marker='o', linestyle='-', color='blue', linewidth=2)
    plt.plot(x, relative_error_qssa, label='Relative Error (X_qssa vs Xini)', marker='s', linestyle='--', color='green', linewidth=2)

    # Set the y-axis to logarithmic scale (if desired, uncomment the next line)
    # plt.yscale('log')

    # Add labels, title, legend, and format the plot
    plt.xlabel('Species', fontsize=14)
    plt.ylabel('Relative Error (%)', fontsize=14)
    plt.title('Relative Error Plot (Xsim and X_qssa vs Xini)', fontsize=16)
    plt.xticks(x, gas.species_names, rotation=90, fontsize=10)  # Rotate x-axis labels for clarity
    plt.legend(fontsize=12)
    plt.grid(axis='y', linestyle='--', linewidth=0.5)  # Add gridlines for better readability

    # Set y-axis limits (adjust as needed)
    plt.ylim([0, 500])  # Set upper limit or adjust based on the data

    # Save the plot
    plt.tight_layout()  # Adjust layout to prevent label overlap
    plt.savefig("figs/rcce_relative_error_plot.png")

def main():
    mech = "src/NH3/NH3_otomo.cti" #"mechs/grimech30_53s325r.cti"
    
    ############cal_init#################
    T, P, X = 1800, 2*ct.one_atm, "H2:0.2, O2:0.2, N2:0.7"
    gas = ct.Solution(mech)
    gas.TPX = T, P, X
    gas.equilibrate('TP')
    Tini, Pini, Xini = gas.TPX
    print_idxs = list(set((Xini > 1e-8).nonzero()[0].tolist()))
    ##########################################
    main_species_names = ["NH3", "O2", "H2", "H2O", "N2", "AR"]
    rcce = RCCE(mech, main_species_names)
    #  T, P, Xini,
    Gini, _ = rcce.get_Gibbs_potential(Tini, Pini, Xini)


############cal_rcce####################
    t0 = time.time()
    rcce.set_conditions(T, P, Xini)
    (Tsim, Psim, Xsim), History = rcce.equilibrium(method="TP")
    t_sim = time.time() - t0
    Gsim, _ = rcce.get_Gibbs_potential(Tsim, Psim, Xsim)
    
    print("Results: Initial, Cantera, Optimal, Newton-")
    rcce.print_X(Xini)
    rcce.print_X(Xsim)
    rcce.print_names()
    
    print("Results: Initial, Cantera, Optimal")
    print("Initial Gibbs Free Energy = %16.6f, (T,P)=(%.3f, %.1f)" % (Gini, Tini, Pini))
    print("Optimal Gibbs Free Energy = %16.6f, (T,P)=(%.3f, %.1f)" % (Gsim, Tsim, Psim))
    print("Optimal cost time =", t_sim)
    plt.plot(History['step'], np.array(History['X'])[:, print_idxs])
    if History['out_step'] == 0 and len(History['step']) > 1 and np.max(History['step']) > 1:
        plt.ylim([1e-8, 1])
        plt.xlim([1, np.max(History['step'])])
        plt.xscale('log')
        plt.yscale('log')
    else:
        plt.ylim([1e-12, 1])
        plt.yscale('log')
    plt.xlabel("Step")
    plt.ylabel("Mole Fraction")
    plt.legend(rcce.sp_names[print_idxs], ncol=5)
    plt.savefig("figs/rcce_test.png")
    
if __name__ == '__main__':
    test_rcce()
    # main()
    
    
    
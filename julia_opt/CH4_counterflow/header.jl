using Arrhenius
using ForwardDiff
using LinearAlgebra
using DiffEqSensitivity
# using DiffEqSensitivity: ForwardDiffSensitivity,ZygoteVJP,InterpolatingAdjoint
using DiffEqSensitivity: ForwardDiffSensitivity
using DifferentialEquations

using Sundials

using Random
using ProgressBars
using DelimitedFiles
using Plots, Colors, Printf, Profile

using Flux
using Flux: crossentropy
using Flux.Losses: mae
using Flux.Optimise: update!
using LatinHypercubeSampling
using Statistics
using DiffEqFlux
using Revise
using CSV
using DataFrames,FilePathsBase  

sensealg = ForwardDiffSensitivity()
# solver = CVODE_BDF()
# solver = KenCarp4()
solver = TRBDF2()
# solver   = Rodas5()    # try 5th-order implicit Runge-Kutta


BLAS.set_num_threads(20)
Threads.nthreads() = 20
# Global settings
abstol    = 1e-9
reltol    = 1e-12

global npr = 1;

# Settings
P = one_atm
# gas = CreateSolution("../../mechanism/mechanism_julia/gri30_217.yaml")
gas = CreateSolution("../../mechanism/mechanism_julia/gri12.yaml")


ns = gas.n_species
nr = gas.n_reactions
# p_true = zeros(nr)
p_true = zeros(nr * npr);

@show gas.species_names

rng = MersenneTwister(0x7777777)


function dudtp!(du, u, p, t)
    T = u[end]
    Y = @view(u[1:ns])
    mean_MW = 1. / dot(Y, 1 ./ gas.MW)
    ρ_mass  = P / R / T * mean_MW
    X       = Y2X(gas, Y, mean_MW)
    C       = Y2C(gas, Y, ρ_mass)

    cp_mole, cp_mass = get_cp(gas, T, X, mean_MW)
    h_mole           = get_H(gas, T, Y, X)
    S0               = get_S(gas, T, P, X)

    # """ for npr = 3
    #     The original Arrhenius equation:     ω =  A *      T^b  * exp(-Ea/RT)
    # kp =  exp( p1 + p2 * log(T) - p3 / RT) = exp(p1) * T^p2 * exp(-p3/RT)
    #     so ω .* kp =  A * exp(p1) * T^(b+p2) *  exp(-(Ea+p3)/RT)
    #     which modifies the parameters to `A * exp(p1)`, `b+p2`, `Ea+p3`
    # """
    if npr==3
        # print("in dudt! npr=$npr (3)\r\n")
        _p = reshape(p, nr, 3)
        kp = @. @views(exp(_p[:, 1] + _p[:, 2] * log(T) - _p[:, 3] * 4184.0 / R / T))
    end
    # """ for npr = 1
    #     The original Arrhenius equation:     ω =  A *      T^b  * exp(-Ea/RT)
    if npr==1
        kp = exp.(p)
    end
    #     so ω .* kp =  A * exp(p) * T^b *  exp(-Ea/RT)
    #     which modifies the parameters to `A * exp(p)`, (`b` and `Ea` not changed)
    # """

    qdot = wdot_func(gas.reaction, T, C, S0, h_mole; get_qdot=true)
    wdot = gas.reaction.vk * (qdot .* kp)

    Ydot = wdot ./ ρ_mass .* gas.MW
    # mask out gradients for main species
    mask = ones(ns)
    mask[main_species_indices] .= 0.0
    Ydot = Ydot .* mask
    Tdot = 0.0
    du .= vcat(Ydot, Tdot)
end


# Function to predict ODE solutions
function predict_ode(u0, p)
    _prob = ODEProblem(dudtp!, u0,tspan, p)
    sol = solve(_prob, solver, saveat=tspan[2],
                reltol=reltol, abstol=abstol,  sensealg=sensealg, verbose=true)
    if sol.retcode != ReturnCode.Success
        println("ODE solver failed during reconstruction")
    end
    return sol[end]
end

# Function to read the species data from the CSV
function read_species_data(file_name)
    # Read the CSV into a DataFrame
    df = CSV.File(file_name) |> DataFrame
    
    # Extract species data columns that start with "Y_" and also include "T" and "P"
    species_columns = filter(name -> startswith(name, "Y_") || name == "T" || name == "P", names(df))
    species_data = select(df, species_columns)
    
    return species_data
end


function get_initial_conditions(species_data::DataFrame,
                                main_species_names::Vector{String})
    ns = gas.n_species
    initial_conditions = Vector{Vector{Float64}}()
    cols = names(species_data)

    for (i,row) in enumerate(eachrow(species_data))
        # start with all zeros
        Y = zeros(ns) 

        # fill in only your main species
        for sp in main_species_names
            col = "Y_" * sp
            if col in cols
                j = species_index(gas, sp)      # species’ index in gas.species_names
                Y[j] = row[col]                 # pull from DataFrameRow by String key
            end
        end
        T0 = row[:T]
        push!(initial_conditions, vcat(Y, T0))
    end

    return initial_conditions
end


function export_predictions(case_name, u0_list, p; suffix="")
    # 1) gather end-state mass fractions
    preds = [ begin
        sol = predict_ode(u0, p)
        sol[1:end-1]                # drop temperature, keep Y₁…Y_ns
    end for u0 in u0_list ]

    # 2) assemble into matrix & DataFrame
    mat = reduce(hcat, preds)'      # (n_grid × n_species)
    df  = DataFrame(mat, header, copycols=true)

    # 3) ensure your targeted recon_dir exists
    mkpath(recon_dir)
    # 4) write out into recon/<species_str>/
    fname = joinpath(recon_dir, "$(case_name)_recon$(suffix).csv")
    CSV.write(fname, df)
end

all_cases = ["CH4_NP_0","CH4_NP_1","CH4_NP_2","CH4_NP_3","CH4_NP_4"]

# Precompute any globals once
all_species_names  = gas.species_names
header = ["Y_" * sp for sp in all_species_names]
main_species_names = ["H2", "O2", "H2O", "CH4", "CO", "CO2", "N2"]
main_species_indices = indexin(main_species_names, gas.species_names)
println("main_species_indices = ", main_species_indices)

tspan = (0.0, 1e0)   # if you need this for export_predictions internally
recon_dir = "recon"

# for case_name in all_cases
#     println("┌─ Processing case: $case_name")
#     # 1) read data
#     case_file = "../../SIM_results/case_CH4_counterflow/$(case_name).csv"
#     local species_data = read_species_data(case_file)
#     # 2) build initial conditions
#     initial_conditions = get_initial_conditions(species_data, main_species_names)
#     # 3) time the prediction export
#     elapsed = @elapsed export_predictions(case_name, initial_conditions, p_true; suffix="_ori")
#     println("└─  export_predictions for $case_name took $(round(elapsed, digits=4)) seconds\n")
# end
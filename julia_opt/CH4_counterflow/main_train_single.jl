########################################################################
# train.jl – fine‑tune reaction log‑factors on Y_H and Y_CH with NeuralODE
#######################################################################
include("header.jl")                     
# sensealg = InterpolatingAdjoint(autojacvec=ZygoteVJP())   # adjoint mode
# ───────────────────────────────────────────────────────────────────
# 0.  USER OPTIONS – adjust in one place
# ───────────────────────────────────────────────────────────────────
train_case_name      = "CH4_NP_0"
train_case_file      = "../../SIM_results/case_CH4_counterflow/$(train_case_name).csv"
all_cases      = ["CH4_NP_0","CH4_NP_1","CH4_NP_2","CH4_NP_3","CH4_NP_4"]
target_species = ["OH",  "H", "O", "CH2"]         # Y_H  &  Y_CH
save_dir       = "finetune/$(train_case_name)"
main_species_names = ["H2","O2","H2O","CH4","CO","CO2","N2"]
# main_species_indices = indexin(main_species_names, gas.species_names)
continue_train  = false
finetune_file = joinpath(save_dir, "p_finetuned.csv")
n_epochs_adam  = 20                 # quick noisy pass
batchsize      = 32                  # grid points / batch
grad_clip      = 1e2
lr_adam        = 1e-2
eps_reg = 1e-5 # 1e-5
tspan = [0,5e-3]
mkpath(save_dir)
# ───────────────────────────────────────────────────────────────────
# 1.  LOAD GRID‑POINT DATA FROM YOUR CSV
# ───────────────────────────────────────────────────────────────────
species_df  = CSV.read(train_case_file, DataFrame)
grid_pts    = size(species_df,1)
OH_init      = species_df[!,:Y_OH]
mask          = (OH_init .>= 3e-4)   .& (OH_init .<= 1)
train_indices = findall(mask)
println("🔍  Training on ", length(train_indices),
        " points where 1e-4 ≤ Y_OH ≤ 1") 

# -------------------------------------------------------------------
# Build u0 with *only* the main‑species mass fractions + temperature
# -------------------------------------------------------------------
"""
    build_u0(row::DataFrameRow) -> Vector{Float64}

Create the initial‑state vector for one grid point.

* Fills `Y[i]` **only** for the species listed in `main_species_names`;
  every other entry stays zero.
* Appends the temperature `row.T` as the last element.
"""
function build_u0(row::DataFrameRow)
    Y = zeros(ns)                                   # all species → 0
    for sp in main_species_names
        col = Symbol("Y_" * sp)                     # e.g. :Y_H2
        if hasproperty(row, col)
            Y[species_index(gas, sp)] = row[col]    # insert the value
        end
    end
    return vcat(Y, row.T)                           # Y₁…Y_ns,  T
end

u0_list = [build_u0(species_df[i,:]) for i=1:grid_pts]

# Pull the experimental target curves (same shape as predictions)

Yexp = zeros(length(target_species), grid_pts)
for (k,sp) in enumerate(target_species)
    col = "Y_"*sp
    Yexp[k,:] .= species_df[!, col]
end

# ───────────────────────────────────────────────────────────────────
# 1b) Compute per‐species scales (Y_max–Y_min) and show them
# ───────────────────────────────────────────────────────────────────
yscale = vec(maximum(Yexp, dims=2) .- minimum(Yexp, dims=2))
# guard against zero range
yscale .= max.(yscale, 1e-16)
println("▶️  Scaling factors (Y_max–Y_min) for each target species:")
for (i,sp) in enumerate(target_species)
    println("   ", sp, ": ", yscale[i])
end
# ───────────────────────────────────────────────────────────────────
# 2.  PREPARE A PROBLEM *PER* GRID‑POINT (saveat = t=1 only)
# ───────────────────────────────────────────────────────────────────
main_species_indices = indexin(main_species_names, gas.species_names)
function model(u0, p)
    # _prob = remake(prob, u0=u0, p=p, tspan=(0.0,1.0))
    _prob = ODEProblem(dudtp!, u0,tspan, p)
    sol   = solve(_prob, solver; saveat=tspan[2],
                  reltol=reltol, abstol=abstol,
                  sensealg = sensealg,          # ←  important
                  verbose=false)
    return sol[end]
end

# Fast wrapper: return *only* Y_H & Y_CH
tgt_idx = indexin(target_species, gas.species_names)
@inline select_targets(u) = @view(u[tgt_idx])
  
function predict_grid(i, p)
    full = model(u0_list[i], p)
    return select_targets(full)
end

# ───────────────────────────────────────────────────────────────────
# 3.  LOSS  (MAE on Y_H & Y_CH, with tiny L2‑reg on p)
# ───────────────────────────────────────────────────────────────────
function loss_ode(p, i_exp)
    # raw and predicted
    ŷ      = predict_grid(i_exp, p)
    y       = Yexp[:, i_exp]

    # normalize each component by its swing
    y_norm  = y  ./ yscale
    ŷ_norm  = ŷ ./ yscale

    return mae(y_norm, ŷ_norm) +
           eps_reg * sum(abs2, p) / length(p)
end

# ───────────────────────────────────────────────────────────────────
# 4.  TRAINING DRIVERS
# ───────────────────────────────────────────────────────────────────
p_pred  = copy(p_true)           # start from reconstructed vector
ps      = Flux.Params([p_pred])  # tell Zygote what to track
loss_hist = Float64[]

function train!(opt; 
                n_epoch=10, 
                batchsize=32, 
                reltol=1e-6, 
                abstol=1e-9)

    println("Start training process!")
    p_pred  = deepcopy(p_true)
    epochs  = ProgressBar(1:n_epoch)
    N       = length(train_indices)

    # buffers indexed over the *subset* 1:N
    grad_norm = zeros(Float64, N)
    losses    = Float64[]
    deltas_p  = Float64[]

    for epoch in epochs
        # shuffle only the selected indices
        perm = randperm(N)
        for j in perm
            i_exp = train_indices[j]

            # forward‑mode gradient of loss_ode(p, i_exp)
            grad = ForwardDiff.gradient(p -> loss_ode(p, i_exp), p_pred)

            # record & clip
            grad_norm[j] = norm(grad)
            if grad_norm[j] > grad_clip
                grad .= (grad ./ grad_norm[j]) .* grad_clip
            end

            # update parameters
            update!(opt, p_pred, grad)
        end

        # compute mean loss & Δp over the *subset*
        ℓ  = mean(loss_ode(p_pred, i) for i in train_indices)
        Δp = maximum(abs.(p_pred .- p_true))

        push!(losses,  ℓ)
        push!(deltas_p, Δp)

        set_description(epochs,
            @sprintf("Epoch %3d | loss %.2e | Δp %.2e | gnorm %.2e",
                     epoch, ℓ, Δp, mean(grad_norm)))

        # snapshot every 5 epochs
        if epoch % 5 == 0
            fname = joinpath(save_dir, "p_epoch_$(epoch).csv")
            CSV.write(fname, DataFrame(p = p_pred))
        end
    end

    # final save & report
    CSV.write(joinpath(save_dir,"p_finetuned.csv"),
              DataFrame(p = p_pred))
    println("✔ Training complete: final loss = $(losses[end]), Δp = $(deltas_p[end])")

    return p_pred
end

#───────────────────────────────────────────────────────────────────
# 5a) Either load existing p_pred, or run train!
#
println(">>> ADAM pass")
if continue_train && isfile(finetune_file)
  println("⚡️  Found existing fine‑tuned parameters at\n    $finetune_file\n   → loading and skipping training")
  df_p   = CSV.read(finetune_file, DataFrame)
  # assume your CSV has a column named `p`
  p_pred = Vector(df_p.p)
else
  p_pred = train!(ADAM(lr_adam); n_epoch=n_epochs_adam, batchsize=batchsize)
end

# 5.1  Save final parameters for the trained case
CSV.write(joinpath(save_dir,"p_finetuned.csv"),
          DataFrame(p = p_pred))

# 5.3  Now loop over all_cases and do the same
for cname in all_cases
    # read that case's CSV
    fname = "../../SIM_results/case_CH4_counterflow/$(cname).csv"
    df    = read_species_data(fname)

    # rebuild its full u0_list
    u0_list_case = get_initial_conditions(df, main_species_names)

    # original vs finetuned
    export_predictions(cname,u0_list_case, p_true;
                       suffix = "_ori")
    export_predictions(cname,u0_list_case, p_pred;
                       suffix = "_finetuned")
    println("✔ Export complete – see $(cname) and finetuned.csv")

end

# 5.4  Finally, plot your loss curve as before
using Plots
pngfile = joinpath(save_dir,"loss_curve.png")
plot(1:length(loss_hist), loss_hist;
     yscale=:log10,
     xlabel="Epoch", ylabel="Loss",
     label="MAE")
savefig(pngfile)
println("✔ Fine‑tuning complete – see $(pngfile) and p_finetuned.csv")



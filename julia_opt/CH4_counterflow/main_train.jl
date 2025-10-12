########################################################################
# main_train.jl – multi‑case fine‑tuning with user‑options up top
########################################################################
include("header.jl")
########################################################################
# Patch Base.log and Base.log10 to return NaN on invalid inputs
########################################################################
# ───────────────────────────────────────────────────────────────────
# 0.  USER OPTIONS – adjust everything here
# ───────────────────────────────────────────────────────────────────
all_cases        = ["CH4_NP_0","CH4_NP_1","CH4_NP_2",
                    "CH4_NP_3","CH4_NP_4"]
train_cases      = all_cases[1:2]            # use first 3 for training
target_species   = ["CH3","H"]                # optimize Y_H & Y_CH
n_epochs_adam    = 10                        # ADAM iters
batchsize        = 2                        # mini‑batch size
grad_clip        = 1e2                       # gradient clipping
lr_adam          = 1e-2                      # ADAM lr
eps_reg          = 1e-6                       # L²‑reg weight
save_dir         = "finetune/multi_case"     # where to dump p's
main_species_names = ["H2","O2","H2O","CH4","CO","CO2","N2"]

mkpath(save_dir)                                    # make sure it exists

# ───────────────────────────────────────────────────────────────────
# 1.  READ & PREPARE DATA FOR *ALL* CASES
# ───────────────────────────────────────────────────────────────────
# We'll store each case's DataFrame and its u0_list
case_to_df = Dict{String,DataFrame}()
case_to_u0 = Dict{String,Vector{Vector{Float64}}}()

for cname in all_cases
    file = "../../SIM_results/case_CH4_counterflow/$(cname).csv"
    df   = read_species_data(file)
    case_to_df[cname] = df
    case_to_u0[cname] = get_initial_conditions(df, all_species_names)
end

ntargets = length(target_species)

u0_list   = Vector{Vector{Float64}}()
Yexp_list = Vector{Vector{Float64}}()

for cname in train_cases
    df   = case_to_df[cname]
    u0_c = case_to_u0[cname]
    for (i,row) in enumerate(eachrow(df))
        # only include if OH mass fraction > 1e-4
        if row[:Y_OH] > 3e-3
            push!(u0_list, u0_c[i])
            # grab each target by symbol
            push!(Yexp_list,
                 [row[ Symbol("Y_" * sp) ] for sp in target_species])
        end
    end
end

# now finalize the filtered dataset
grid_pts = length(u0_list)
Yexp     = reduce(hcat, Yexp_list)
println("🔍  Training on ", grid_pts,
        " grid‑points (Y_OH > 1e‑4) from cases ",
        join(train_cases, ", "))

# ───────────────────────────────────────────────────────────────────
# 1b) Compute per‑species scales (max–min) to normalize the MAE
# ───────────────────────────────────────────────────────────────────
yscale = vec(maximum(Yexp, dims=2) .- minimum(Yexp, dims=2))
# guard against zero ranges
yscale .= max.(yscale, 1e-16)


# ───────────────────────────────────────────────────────────────────
# 2.  DEFINE THE ODE + PREDICTOR
# ───────────────────────────────────────────────────────────────────
main_species_indices = indexin(main_species_names, gas.species_names)

function model(u0, p)
    _prob = remake(prob, u0=u0, p=p, tspan=(0.0,1.0))
    sol   = solve(_prob, solver; saveat=[1.0],
                  reltol=reltol, abstol=abstol,
                #   sensealg = sensealg,          # ←  important
                  verbose=false)
    return sol[end]
end
tgt_idx = indexin(target_species, gas.species_names)
@inline select_targets(u) = @view(u[tgt_idx])

function predict_grid(i, p)
    full = model(u0_list[i], p)
    return select_targets(full)
end

# ───────────────────────────────────────────────────────────────────
# 3.  LOSS (MAE on Y_H & Y_CH + tiny L²‑reg on p)
# ───────────────────────────────────────────────────────────────────

function loss_ode(p, i_exp)
    ŷ = predict_grid(i_exp, p)
    y  = Yexp[:, i_exp]
    return mae(y, ŷ) + eps_reg*sum(abs2, p)/length(p)
end

function loss_ode(p, i_exp)
    ŷ     = predict_grid(i_exp, p)
    y      = Yexp[:, i_exp]
    # normalize each component by its precomputed scale
    y_norm  = y  ./ yscale
    ŷ_norm  = ŷ ./ yscale
    return mae(y_norm, ŷ_norm) +
           eps_reg * sum(abs2, p) / length(p)
end

# ───────────────────────────────────────────────────────────────────
# 4.  TRAINING LOOP (ForwardDiff + Flux.Optimise.update!)
# ───────────────────────────────────────────────────────────────────
function train!(opt; n_epoch=50, batchsize=32)
    println("▶️  Start training process!")
    p_pred    = deepcopy(p_true)
    ps = Flux.params(p_pred)  # 定义参数集合
    epochs    = ProgressBar(1:n_epoch)
    N         = grid_pts
    losses    = Float64[]
    deltas_p  = Float64[]
    grad_norm = zeros(Float64, N)
    for epoch in epochs
        for j in randperm(N)
            grad = ForwardDiff.gradient(p-> loss_ode(p,j), p_pred)
            grad_norm[j] = norm(grad)
            if grad_norm[j] > grad_clip
                grad .= (grad ./ grad_norm[j]) .* grad_clip
            end
            update!(opt, p_pred, grad)
            println("iteration finished for j", j)
        end

        # epoch metrics
        ℓ  = mean(loss_ode(p_pred,i) for i in 1:N)
        Δp = maximum(abs.(p_pred .- p_true))
        push!(losses, ℓ);  push!(deltas_p, Δp)

        # update bar description
        set_description(epochs,
            @sprintf("Epoch %3d | loss %.2e | Δp %.2e",
                     epoch, ℓ, Δp))

        # snapshot every 10 epochs
        if epoch % 10 == 0
            CSV.write(joinpath(save_dir,
                               "p_epoch_$(epoch).csv"),
                      DataFrame(p=p_pred))
        end
    end

    # final save
    CSV.write(joinpath(save_dir,"p_finetuned.csv"),
              DataFrame(p=p_pred))
    println("✔ Training done: loss=$(losses[end]) Δp=$(deltas_p[end])")
    return p_pred, losses
end

# ───────────────────────────────────────────────────────────────────
# 5.  RUN TRAINING
# ───────────────────────────────────────────────────────────────────
println("=== ADAM PASS ===")
p1, losses = train!(ADAM(lr_adam);
                    n_epoch=n_epochs_adam,
                    batchsize=batchsize)

# ───────────────────────────────────────────────────────────────────
# 6.  EXPORT PREDICTIONS FOR *ALL* CASES
# ───────────────────────────────────────────────────────────────────
for cname in all_cases
    export_predictions(case_to_u0[cname],
                       p_true;
                       suffix="_ori_$(cname)")
    export_predictions(case_to_u0[cname],
                       p1;
                       suffix="_finetuned_$(cname)")
end

# ───────────────────────────────────────────────────────────────────
# 7.  PLOT LOSS ‑‑ optional
# ───────────────────────────────────────────────────────────────────
using Plots
pngfile = joinpath(save_dir,"loss_curve.png")
plot(1:length(losses), losses;
     yscale=:log10,
     xlabel="Epoch", ylabel="Loss",
     label="MAE")
savefig(pngfile)
println("✔ Saved loss curve: ", pngfile)

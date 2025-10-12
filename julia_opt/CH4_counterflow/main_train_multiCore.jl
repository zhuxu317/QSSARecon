########################################################################
# train.jl – fine-tune reaction log-factors on selected species with NeuralODE
########################################################################
include("header.jl")

using Base.Threads                    # for @threads, nthreads()
println("▶️  Running with $(nthreads()) threads")

# ───────────────────────────────────────────────────────────────────
# 0.  USER OPTIONS
# ───────────────────────────────────────────────────────────────────
train_case_name    = "CH4_NP_0"
train_case_file    = "../../SIM_results/case_CH4_counterflow/$(train_case_name).csv"
all_cases          = ["CH4_NP_0","CH4_NP_1","CH4_NP_2","CH4_NP_3","CH4_NP_4"]

target_species     = ["CH2"]
species_str        = join(target_species, "_")

save_dir           = "finetune/$(train_case_name)"
main_species_names = ["H2","O2","H2O","CH4","CO","CO2","N2"]

continue_train     = false
finetune_file      = joinpath(save_dir, "p_finetuned_$(species_str).csv")
eps_reg            = 0
tspan              = [0.0, 1e1]
grad_clip          = 1e2
mkpath(save_dir)
recon_dir    = joinpath("recon", species_str)
mkpath(recon_dir)

# ───────────────────────────────────────────────────────────────────
# 1.  LOAD DATA & INITIAL CONDITIONS
# ───────────────────────────────────────────────────────────────────
species_df   = CSV.read(train_case_file, DataFrame)
grid_pts     = size(species_df,1)
# CH3_init     = species_df[!,:Y_CH3]
# mask         = (CH3_init .>= 2.7e-4) .& (CH3_init .<= 1)
CH2_init     = species_df[!,:Y_CH2]
mask         = (CH2_init .>= 1e-6) .& (CH2_init .<= 1)

train_indices = findall(mask)
println("🔍  Training on ", length(train_indices),
        " points where 2.5e-4 ≤ Y_CH3 ≤ 1")
batchsize = length(train_indices)

function build_u0(row::DataFrameRow)
    Y = zeros(ns)
    for sp in main_species_names
        col = Symbol("Y_" * sp)
        if hasproperty(row, col)
            Y[species_index(gas, sp)] = row[col]
        end
    end
    return vcat(Y, row.T)
end

u0_list = [build_u0(species_df[i,:]) for i in 1:grid_pts]

Yexp = zeros(length(target_species), grid_pts)
for (k,sp) in enumerate(target_species)
    Yexp[k, :] .= species_df[!, Symbol("Y_"*sp)]
end

yscale = vec(maximum(Yexp, dims=2) .- minimum(Yexp, dims=2))
yscale .= max.(yscale, 1e-10)
println("▶️  Scaling factors (Y_max–Y_min):")
for (i,sp) in enumerate(target_species)
    println("   ", sp, ": ", yscale[i])
end

# ───────────────────────────────────────────────────────────────────
# 2.  ODE PROBLEM SETUP
# ───────────────────────────────────────────────────────────────────
main_species_indices = indexin(main_species_names, gas.species_names)

function model(u0, p)
    prob = ODEProblem(dudtp!, u0, tspan, p)
    sol  = solve(prob, solver;
                 saveat = tspan[2],
                 reltol = reltol, abstol = abstol,
                 sensealg = sensealg,
                 verbose = false)
    return sol[end]
end

tgt_idx = indexin(target_species, gas.species_names)
@inline select_targets(u) = @view(u[tgt_idx])
predict_grid(i,p) = select_targets(model(u0_list[i],p))

# ───────────────────────────────────────────────────────────────────
# 3.  LOSS PIECES
# ───────────────────────────────────────────────────────────────────
mae_ode(p,i_exp) = begin
    ŷ     = predict_grid(i_exp, p)
    y      = Yexp[:, i_exp]
    y_norm = y  ./ yscale
    ŷ_norm = ŷ ./ yscale
    mae(y_norm, ŷ_norm)
end

reg_ode(p) = eps_reg * sum(abs2, p) / length(p)
loss_ode(p,i_exp) = mae_ode(p,i_exp) + reg_ode(p)

# ───────────────────────────────────────────────────────────────────
# 4.  TRAINING DRIVER
# ───────────────────────────────────────────────────────────────────
function train!(p0, opt; n_epoch=10, batchsize=32, reltol=1e-6, abstol=1e-9)
    println("Start training process!")
    p_pred = deepcopy(p0)
    loss_hist = Float64[]
    N = length(train_indices)

    # storage for per-epoch mean‐gradients
    gbar_hist = Array{Float64,2}(undef, length(p0), n_epoch)



    for epoch in 1:n_epoch
        perm = randperm(N)

        # batched, threaded gradient updates
        for bstart in 1:batchsize:N
            batch = perm[bstart:min(bstart+batchsize-1, N)]
            grads = Vector{Vector{Float64}}(undef, length(batch))

            @threads for t in eachindex(batch)
                j       = batch[t]
                grads[t]= ForwardDiff.gradient(p->loss_ode(p,j), p_pred)
            end

            Gmat   = reduce(hcat, grads)
            ḡ      = vec(mean(Gmat, dims=2))
            norm_g = norm(ḡ)

            if norm_g > grad_clip
                @printf("⚠️  Gradient norm %.2e exceeds clip %.2e\n", norm_g, grad_clip)
                println("Gradient (first 10): ", ḡ[1:min(end,10)])
                ḡ .= (ḡ ./ norm_g) .* grad_clip
            end
            p_pred = update!(opt, p_pred, ḡ)
            # update!(opt, p_pred, ḡ)
        end

        mae_term   = mean(mae_ode(p_pred, i) for i in train_indices)
        reg_term   = reg_ode(p_pred)
        total_loss = mae_term + reg_term

        @printf("Epoch %3d | MAE = %.2e | Reg = %.2e | Total = %.2e\n",
                epoch, mae_term, reg_term, total_loss)
        push!(loss_hist, total_loss)

        if epoch % 5 == 0
            CSV.write(joinpath(save_dir, "p_epoch_$(epoch).csv"),
                      DataFrame(p = p_pred))
        end
    end
    df_changes = DataFrame(
        index    = 1:length(p0),
        p_init   = p0,
        p_pred   = p_pred,
        delta    = p_pred .- p0,
        gbar_end = grads[:, end]
    )
   CSV.write(joinpath(save_dir, "p_changes.csv"), df_changes)
   println("✔  Saved parameter changes + final ḡ → $(joinpath(save_dir, "p_changes.csv"))")


    CSV.write(finetune_file, DataFrame(p = p_pred))
    println("✔ Training complete: final loss = $(loss_hist[end])")
    return p_pred, loss_hist
end

# ───────────────────────────────────────────────────────────────────
# 5a) LOAD / CONTINUE / OR START FRESH
# ───────────────────────────────────────────────────────────────────
println(">>> ADAM pass")

# determine initial parameters
p_init = p_true
if continue_train && isfile(finetune_file)
    println("⚡️  Loading existing parameters from $finetune_file")
    p_init = Vector(CSV.read(finetune_file, DataFrame).p)
end

# now chain your training stages
opt1 = ADAMW(1e-1, (0.9,0.999), 1f-6)
p1, hist1 = train!(p_init, opt1; n_epoch=5, batchsize=batchsize)

opt2 = ADAMW(5e-2, (0.9,0.999), 1f-6)
p2, hist2 = train!(p1,    opt2; n_epoch=5, batchsize=batchsize)

opt3 = ADAMW(5e-2, (0.9,0.999), 1f-6)
p3, hist3 = train!(p1,    opt3; n_epoch=5, batchsize=batchsize)
# final result
p_pred = p3
loss_hist = vcat(hist1, hist2, hist3)

# ───────────────────────────────────────────────────────────────────
# 5b) EXPORT PREDICTIONS
# ───────────────────────────────────────────────────────────────────
for cname in all_cases
    df  = read_species_data("../../SIM_results/case_CH4_counterflow/$(cname).csv")
    u0c = get_initial_conditions(df, main_species_names)
    export_predictions(cname, u0c, p_true;  suffix="_ori")
    export_predictions(cname, u0c, p_pred;  suffix="_finetuned")
    println("✔ Exported for $cname → $recon_dir")
end

# ───────────────────────────────────────────────────────────────────
# 5c) PLOT LOSS HISTORY
# ───────────────────────────────────────────────────────────────────
pngfile = joinpath(save_dir, "loss_curve_$(species_str).png")
plot(1:length(loss_hist), loss_hist;
     yscale = :log10,
     xlabel = "Epoch", ylabel = "Loss",
     label   = "Total")
savefig(pngfile)
println("✔ Fine-tuning complete – see $pngfile and $finetune_file")

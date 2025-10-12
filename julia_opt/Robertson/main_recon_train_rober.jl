###########################################################################
# main_recon_train_rober.jl
#
# A professional, modular version for training a neural ODE by 
# reconstructing state trajectories and fine-tuning parameters.
###########################################################################

###############################
# 1. Load Dependencies & Header
###############################
include("header.jl")

# sensealg = InterpolatingAdjoint(autojacvec=ReverseDiffVJP(false))
# sensealg = BacksolveAdjoint(autojacvec=ReverseDiffVJP(false))
# sensealg = ForwardSensitivity(autojacvec=true)

# Custom vector norm (shallow version for regularization)
vecnorm(x) = sum(abs2, x) / length(x)

###############################
# 2. Define ODE Functions
###############################

# Main ODE function
function f(y, p, t)
    k = [0.04, 3e7, 1e4] .* p
    dydt1 = -k[1] * y[1] + k[3] * y[2] * y[3]
    dydt2 =  k[1] * y[1] - k[2] * y[2]^2 - k[3] * y[2] * y[3]
    dydt3 =  k[2] * y[2]^2
    return [dydt1, dydt2, dydt3]
end

# Reconstruction ODE function: Only update y₂; keep y₁ and y₃ fixed.
function recon_f(y, p, t)
    k = [0.04, 3e7, 1e4] .* p
    dydt1 = 0
    dydt2 = k[1] * y[1] - k[2] * y[2]^2 - k[3] * y[2] * y[3]
    dydt3 = 0
    return [dydt1, dydt2, dydt3]
end

###############################
# 3. ODE Solution & Reconstruction Helpers
###############################

# Predict full ODE trajectory given initial condition `u0` and parameters `p`
function predict_ode(u0, p; sample = samplesteps)
    _prob = remake(prob, u0=u0, p=p, tspan=[0, tsteps[sample]])
    _sol = solve(_prob, solver, saveat=tsteps,
                 reltol=reltol, abstol=abstol, verbose=false)
    pred = Array(_sol)
    if _sol.retcode == ReturnCode.Success
        # Optionally do nothing
    else
        println("ODE solver failed")
    end
    return pred
end

# Run a reconstruction starting from a single state vector, forcing y₂ = 0
function recon_ode(y_vec, p; sample = tsteps_numbers)
    # Enforce reconstruction initial condition: use y₁ and y₃ from y_vec, zero out y₂.
    y0 = [y_vec[1], 0.0, y_vec[3]]
    
    # Create a local time grid for reconstruction
    local tsteps_local = 10 .^ range(log10(1e-4), log10(1e4), length=sample)
    local tspan_local = (0.0, tsteps_local[end] + 1e-3)
    
    # Define and solve the ODE problem using recon_f
    local prob_local = ODEProblem(recon_f, y0, tspan_local, p)
    local sol = solve(prob_local, solver, saveat=tsteps_local[end],
                      reltol=reltol, abstol=abstol, verbose=false)
    
    if sol.retcode != ReturnCode.Success
        println("ODE solver failed during reconstruction")
    end
    return sol[end]  # Return final state (3-element vector)
end

# Apply recon_ode columnwise to reconstruct the entire data set.
function rober_reconstruct_solution(y_true, datasize, recon_f, p_true)
    n = size(y_true, 2)  # number of columns (samples)
    y_recon = []         # store reconstructed columns
    for i in 1:n
        y_vec = y_true[:, i]
        y_pred = recon_ode(y_vec, p_true; sample=datasize)
        push!(y_recon, y_pred)
    end
    y_recon_matrix = hcat(y_recon...)
    println("size of y_recon = ", size(y_recon_matrix))
    return y_recon_matrix
end

###############################
# 4. Loss Function & Training
###############################
# Function to save DataFrame to CSV
function save_to_csv(filename::String, tsteps::Vector, y_true::Array, y_recon::Array)
    # Create a DataFrame for saving results
    N = length(tsteps)
    df = DataFrame(
        t = tsteps,
        y_true_1 = y_true[1, :],
        y_true_2 = y_true[2, :],
        y_true_3 = y_true[3, :],
        y_recon_1 = y_recon[1, :],
        y_recon_2 = y_recon[2, :],
        y_recon_3 = y_recon[3, :]
    )
    
    # Write the DataFrame to CSV
    CSV.write(filename, df)
end

# Loss function: Only compare the second state since y₁ and y₃ are identical.
function loss_ode(p, i_exp; abstol=1e-12, sample=tsteps_numbers)
    # Use the i_exp-th column from u0_list (each column is a training example)
    y_recon = recon_ode(u0_list[:, i_exp], p; sample=sample)
    loss = mae(u0_list[2, i_exp], y_recon[2])  # Compare only y₂
    return loss
end

function train(opt; n_epoch=10, batchsize=50, reltol=1e-6, abstol=1e-9)
    println("Start training process!")
    
    # Compute initial reconstruction on noisy data using the initial parameters
    y_recon = rober_reconstruct_solution(y_noise, tsteps_numbers, recon_f, p_init)
    
    # Save initial state to CSV (epoch 0)
    save_to_csv("figures/train_interval_epoch_0.csv", tsteps, y_true, y_recon)
    
    p_pred = deepcopy(p_init)
    epochs = ProgressBar(1:n_epoch)
    loss_epoch = zeros(Float64, n_exp)
    grad_norm = zeros(Float64, n_exp_train)
    
    # Create containers to store losses for each epoch
    losses_y_train = Float64[]
    losses_y_valid = Float64[]
    losses_p = Float64[]

    for epoch in epochs
        # Update parameters for training examples (random permutation)
        for i_exp in randperm(n_exp_train)
            grad = ForwardDiff.gradient(
                        x -> loss_ode(x, i_exp; abstol=abstol, sample=tsteps_numbers),
                        p_pred)
            grad_norm[i_exp] = norm(grad, 2)
            if grad_norm[i_exp] > grad_max
                grad = grad ./ grad_norm[i_exp] .* grad_max
            end
            update!(opt, p_pred, grad)
        end

        # Compute loss for all examples
        for i_exp in 1:n_exp
            loss_epoch[i_exp] = loss_ode(p_pred, i_exp; abstol=abstol, sample=tsteps_numbers)
        end
        loss_y_train = mean(loss_epoch[1:n_exp_train])
        loss_y_valid = mean(loss_epoch[n_exp_train+1:end])
        loss_p = mae(p_pred, p_true)
        push!(history_p_pred, deepcopy(p_pred))
        push!(losses_y_train, loss_y_train)
        push!(losses_y_valid, loss_y_valid)
        push!(losses_p, loss_p)

        # Inside your training loop, replace the plotting section every 5 epochs with:
        if epoch % 5 == 0
            # Compute the reconstruction for all training examples
            y_recon = rober_reconstruct_solution(u0_list, tsteps_numbers, recon_f, p_pred)
            
            # Save the results for the current epoch to CSV
            csv_filename = "figures/train_interval_epoch_$(epoch).csv"
            save_to_csv(csv_filename, tsteps, y_true, y_recon)
        end
        
        # Set the description for the progress bar (optional for logging):
        set_description(epochs, string(@sprintf("Loss ytrain %.3e yvalid %.3e p %.3e gnorm %.3e",
                    loss_y_train, loss_y_valid, loss_p, mean(grad_norm))))
    end
    
    # Save final losses to a CSV file
    final_losses_df = DataFrame(
        epoch = n_epoch,
        losses_y_train = losses_y_train,
        losses_y_valid = losses_y_valid,
        losses_p = losses_p
    )
    
    CSV.write("figures/final_losses.csv", final_losses_df)
    
    return p_pred   # Return final parameters
end

###############################
# 5. Main Script: Data Generation, Training, and Final Plotting
###############################

# Global settings
abstol    = 1e-9
reltol    = 1e-12
p_true    = [1, 1, 1]
y0        = [1.0, 0.0, 0.0]
tsteps_numbers = 50
tsteps    = 10 .^ range(log10(1e-1), log10(1e5), length=tsteps_numbers)
tspan     = (0.0, tsteps[end] + 1e-3)
solver    = KenCarp4()  # Alternatively, use Rosenbrock23()

# Define the ODE problem and compute the true solution.
prob      = ODEProblem(f, y0, tspan, p_true)
y_true    = predict_ode(y0, p_true; sample=tsteps_numbers)

# Optionally add noise (set noise_level > 0 for noisy data)
noise_level = 0
rng         = MersenneTwister(Int32(floor(1e7*noise_level)))
scale       = vec(maximum(y_true, dims=2))
y_noise   = y_true + noise_level .* (rand(rng, length(y0), length(tsteps)) .- 0.5) .* scale

# Plot initial reconstruction using p_init
weights   = [1, 2e4, 1]
xscale    = :log10
n_epoch   =  50
n_exp     = tsteps_numbers
n_exp_train = 40
n_exp_valid = tsteps_numbers - n_exp_train
qssa_time_end = 1e4

p_init    = [0.5, 2, 0.9]
y1_true, y2_true, y3_true = y_true[1, :], y_true[2, :], y_true[3, :]

# Compute reconstruction on noisy data using the initial parameters.
y_recon  = rober_reconstruct_solution(y_noise, tsteps_numbers, recon_f, p_init)
h_init   = valid(tsteps, y_true, y_recon; xscale=xscale)
Plots.savefig(h_init, "figures/recon_init.png")

# Set training dataset; each column in u0_list is one training example.
u0_list = y_true

# Prepare history containers for training metrics.
losses_y_train  = Vector{Float64}()
losses_y_valid  = Vector{Float64}()
losses_p        = Vector{Float64}()
history_p_pred  = Vector{Array{Float64}}()

# Set optimizer
opt = ADAMW(0.001, (0.9, 0.999), 1.f-6)
# opt = ADAMW(0.05, (0.9, 0.999), 1.f-6)

# Run the training process.
p_final = train(opt; n_epoch=n_epoch, batchsize=10, reltol=reltol, abstol=abstol)

# Final reconstruction and comparison plot.
y_recon_final = rober_reconstruct_solution(u0_list, tsteps_numbers, recon_f, p_final)
# Here we compare: true data, initial reconstruction, and final reconstruction.
h_final = compare(tsteps, y_true, y_recon, y_recon_final; xscale=xscale)
Plots.savefig(h_final, "figures/final_comparison.png")
include("header.jl")
include("Robertson.jl")

gr()  # 切换绘图后端到 GR
mkpath("figures")  # 确保保存目录存在


solver = Rosenbrock23();

# sensalg = InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true)) # 0.034 s
# sensalg = ForwardDiffSensitivity()                              # 0.00733 s
# sensalg = DiffEqSensitivity.ForwardDiffSensitivity()
# sensalg = DiffEqFlux.ForwardDiffSensitivity()
# sensalg = InterpolatingAdjoint()                                # 0.51 s
# sensalg = QuadratureAdjoint()                                   # 0.53 s


function loss(p; y_train=y_true)
    return sum(abs2, (predict(p) - y_train) ./ scale) / length(y_train);
end

function train(p_init, y_noise; n_epoch=100, opt = ADAMW(0.1,(0.9,0.999),1e-6))
    p_pred = deepcopy(p_init);
    y_pred = predict(p_pred)
    losses_y = Vector{Float64}([loss(p_init; y_train=y_noise)]);
    history_p = Vector{Array{Float64}}([p_init]);
    epochs = ProgressBar(1:n_epoch);
    for epoch in epochs

        grad = Flux.gradient(x -> loss(x; y_train=y_noise), p_pred)[1]
        update!(opt, p_pred, grad)

        loss_y = loss(p_pred; y_train=y_noise)
        push!(losses_y, deepcopy(loss_y))
        push!(history_p, deepcopy(p_pred))
        set_description(epochs, string(@sprintf("loss_y = %.3e, loss_p %.3e grad %.3e",
                    loss_y, norm(p_pred.-1)^2, norm(grad))))
    end 
    return losses_y, history_p;
end

noise_level = 5e-2;
rng = MersenneTwister(Int32(floor(1e7*noise_level)));
y_noise = y_true + noise_level .* (rand(rng, length(y0), length(tsteps)).-0.5) .* scale;

# initialize


p_init = exp.((rand(rng, length(p_true)) .* 2 .- 1) .* 2);

println("cal_p_init=", p_init)

# p_init = [1.0, 1.0, 1.0];

y_init = predict(p_init);

# ===== 修改下方绘图代码 =====
# 初始化结果图
h = valid(tsteps, y_noise, y_init; xscale=xscale);
Plots.savefig(h, string(@sprintf("figures/%s_noise=%.0e_init.png", casename, noise_level)));

# 训练过程
losses_y, history_p = train(p_init, y_noise; n_epoch=200);
losses_p = vcat(sum(abs2.(hcat(history_p...) .- p_true), dims=1)...);

# 损失曲线图
h1 = plot(xlabel="Epochs", ylabel="y loss", size=(400,200), legend=false);
plot!(losses_y, yscale=:log10, line=(1, :solid), color=:black);
h2 = plot(xlabel="Epochs", ylabel="p loss", size=(400,200), legend=false);
plot!(losses_p, yscale=:log10, line=(1, :solid), color=:black);
h = plot(h1, h2, layout=(2,1), size=(400,400), framestyle=:box);
Plots.savefig(h, string(@sprintf("figures/%s_noise=%.0e_loss.png", casename, noise_level)));

# 预测结果图
p_pred = history_p[end];
y_pred = predict(p_pred);
h = valid(tsteps, y_noise, y_pred; xscale=xscale);
Plots.savefig(h, string(@sprintf("figures/%s_noise=%.0e_pred.png", casename, noise_level)));

# 对比图（已修改保存方式）


# 保存图像
h = compare_Robertson(tsteps, y_noise, y_init, y_pred; xscale=xscale)
Plots.savefig(h, string(@sprintf("figures/%s_noise=%.0e_comp.png", casename, noise_level)))
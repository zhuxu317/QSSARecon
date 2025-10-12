casename = "Robertson"

# the ODE function
function f(y, p, t)
    k = [0.04, 3e7, 1e4] .* p
    dydt1 = -k[1]*y[1]+k[3]*y[2]*y[3]
    dydt2 =  k[1]*y[1]-k[2]*y[2]^2-k[3]*y[2]*y[3]
    dydt3 =  k[2]*y[2]^2
    return [dydt1, dydt2, dydt3]
end

# the recon ODE function
function recon_f(y, p, t)
    k = [0.04, 3e7, 1e4] .* p
    dydt1 =  0
    dydt2 =  k[1]*y[1]-k[2]*y[2]^2-k[3]*y[2]*y[3]
    dydt3 =  0
    return [dydt1, dydt2, dydt3]
end

p_true = [1, 1, 1];
y0 = [1.0, 0.0, 0.0];
datasize = 50;#50;
# tsteps = 10 .^ range(log10(1e-4), log10(1e8), length=datasize);
tsteps = 10 .^ range(log10(1e-4), log10(1e3), length=datasize);

tspan = (0.0, tsteps[end]+1e-3);

# solver = DifferentialEquations.CVODE_BDF()
# solver = KenCarp4();
solver = Rosenbrock23();
function predict(p)
    _sol = solve(prob, solver, p=p, saveat=tsteps)
    # println("ODE solver status=", _sol.retcode)  # 打印实际返回值
    _pred = Array(_sol)
    if _sol.retcode == ReturnCode.Success
        return _pred
    else
        println("ODE solver failed with code: ", _sol.retcode)
        return nothing  
    end
end


function loss(p; y_train=y_true)
    return sum(abs2, (predict(p) - y_train) ./ scale) / length(y_train);
end

# get data with p_true
println("tspan = ", tspan)  # 确认时间跨度合理（如非负）

prob = ODEProblem(f, y0, tspan, p_true);
y_true = predict(p_true);

scale = vec(maximum(y_true, dims=2));

# for valid plot
weights = [1, 2e4, 1];
xscale = :log10;

# for training
n_epoch = 200;

# fnt = Plots.font("DejaVu Sans", 10.0)
# default(titlefont=fnt, guidefont=fnt, tickfont=fnt, legendfont=fnt)


function compare_Robertson(tsteps, y_noise, y_init, y_pred; xscale=:identity)
    h = plot(xlabel="t", ylabel="Concentration", legend=:left)
    
    # 绘制真实数据
    plot!(tsteps, y_true[1,:], label="True A", line=(2, :dash, :black))
    plot!(tsteps, y_true[2,:], label="True B", line=(2, :dash, :black))
    plot!(tsteps, y_true[3,:], label="True C", line=(2, :dash, :black))
    
    # 绘制噪声数据
    scatter!(tsteps, y_noise[1,:], label="Noisy A", markersize=3, color=1)
    scatter!(tsteps, y_noise[2,:], label="Noisy B", markersize=3, color=2)
    scatter!(tsteps, y_noise[3,:], label="Noisy C", markersize=3, color=3)
    
    # 绘制预测数据
    plot!(tsteps, y_pred[1,:], label="Pred A", line=3, color=1)
    plot!(tsteps, y_pred[2,:], label="Pred B", line=3, color=2)
    plot!(tsteps, y_pred[3,:], label="Pred C", line=3, color=3)
    
    # 设置坐标轴（关键修改）
    plot!(
        xscale=xscale,
        # yscale=:log10,
        ylims=(-0.1, 1),        # 根据数据范围调整
        # yticks=10.0.^(-6:2:2),   # 自定义对数刻度
        minorticks=true           # 启用次刻度
    )
    return h
end
println("[info] Robertson ODE loaded")

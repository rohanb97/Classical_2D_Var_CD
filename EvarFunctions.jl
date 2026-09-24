include("PolynomialClass.jl")
using DifferentialEquations, OrdinaryDiffEq, LinearAlgebra, Optim
using Plots
using Base.Threads
using Statistics
using Profile
using Interpolations
using StaticArrays
using DelimitedFiles
using Random
using Dates
# --------------------------
# Ensemble Sampling function
# --------------------------

struct HamDerivs
    dHdx::Function
    dHdpx::Function
    dHdy::Function
    dHdpy::Function
end

function build_ham_derivs(H_Symb)
    HamDerivs(
        polynomial_to_function(Polynomial_Derivative(H_Symb, :x)),
        polynomial_to_function(Polynomial_Derivative(H_Symb, :px)),
        polynomial_to_function(Polynomial_Derivative(H_Symb, :y)),
        polynomial_to_function(Polynomial_Derivative(H_Symb, :py))
    )
end

# Sample N points from the microcanonical ensemble
function MicroSample(N, V, E0;βinitial = 0.0)
    sampleList =Vector{Vector{Float64}}()
    for i in 1:N
        while true
            x = (2*rand()-1) * sqrt(2 * E0)
            y = (2*rand()-1)* sqrt(2 * E0)
            if V(x,y,βinitial)<= E0
                phi = 2 * pi * rand()
                p = sqrt(2 * (E0 - V(x,y,βinitial)))
                px = p * cos(phi)
                py = p * sin(phi)
                push!(sampleList,[x, px, y, py])
                break
            end
        end
    end
    println("Sampled $N points.")
    return sampleList
end

# ----------------------------
# Equation of motion functions
# ----------------------------

function beta(t,τ;βrange=(0,1))
    return βrange[1] + (βrange[2] - βrange[1]) * sin(pi/2*sin(pi/2*t/τ)^2)^2
end

function betadot(t,τ;βrange=(0,1))
    return pi^2/(4*τ)*sin(pi*t/τ)*sin(pi*sin(pi/2*t/τ)^2)*(βrange[2] - βrange[1])
end

function EoM!(du, u, p, t)
    τ, _, _,βrange, Hderivs = p
    x, px, y, py = u
    β = beta(t, τ, βrange=βrange)

    du[1] = Hderivs.dHdpx(x, px, y, py, β)
    du[2] = -Hderivs.dHdx(x, px, y, py, β)
    du[3] = Hderivs.dHdpy(x, px, y, py, β)
    du[4] = -Hderivs.dHdy(x, px, y, py, β)
end

function prob_func(prob, i, repeat)
    remake(prob, u0=ensemble[i])
end

function AGP_EoM!(du, u, p, t)
    τ, γlist, dAGP_dx_func, βrange, Hderivs = p
    x, px, y, py = u

    β = beta(t, τ, βrange=βrange)
    βdot = betadot(t, τ, βrange=βrange)

    du[1] = Hderivs.dHdpx(x, px, y, py, β)
    du[2] = -Hderivs.dHdx(x, px, y, py, β)
    du[3] = Hderivs.dHdpy(x, px, y, py, β)
    du[4] = -Hderivs.dHdy(x, px, y, py, β)

    if !isempty(γlist)
        for i in eachindex(γlist)
            γ = γlist[i](β)
            du[1] += βdot * γ * dAGP_dx_func[i,1](x, px, y, py, β)
            du[2] -= βdot * γ * dAGP_dx_func[i,2](x, px, y, py, β)
            du[3] += βdot * γ * dAGP_dx_func[i,3](x, px, y, py, β)
            du[4] -= βdot * γ * dAGP_dx_func[i,4](x, px, y, py, β)
        end
    end
end

function Ramp_EoM!(du, u, p, t)
    τ, γlist, dAGP_dx_func, βrange, Hderivs = p
    x, px, y, py = u

    β = beta(t, τ, βrange=βrange)
    βdot = betadot(t, τ, βrange=βrange)

    du[1] = Hderivs.dHdpx(x, px, y, py, β)
    du[2] = -Hderivs.dHdx(x, px, y, py, β)
    du[3] = Hderivs.dHdpy(x, px, y, py, β)
    du[4] = -Hderivs.dHdy(x, px, y, py, β)

    if !isempty(γlist)
        for i in eachindex(γlist)
            γ = γlist[i](β)
            du[1] += βdot * γ * dAGP_dx_func[i,1](x, px, y, py, β)
            du[2] -= βdot * γ * dAGP_dx_func[i,2](x, px, y, py, β)
            du[3] += βdot * γ * dAGP_dx_func[i,3](x, px, y, py, β)
            du[4] -= βdot * γ * dAGP_dx_func[i,4](x, px, y, py, β)
        end
    end
end

function Explore_EoM!(du, u, p, t)
    τ, γlist, dAGP_dx_func, βrange, Hderivs = p
    x, px, y, py = u

    β = last(βrange)
    βdot = 0

    du[1] = Hderivs.dHdpx(x, px, y, py, β)
    du[2] = -Hderivs.dHdx(x, px, y, py, β)
    du[3] = Hderivs.dHdpy(x, px, y, py, β)
    du[4] = -Hderivs.dHdy(x, px, y, py, β)

    if !isempty(γlist)
        for i in eachindex(γlist)
            γ = γlist[i](β)
            du[1] += βdot * γ * dAGP_dx_func[i,1](x, px, y, py, β)
            du[2] -= βdot * γ * dAGP_dx_func[i,2](x, px, y, py, β)
            du[3] += βdot * γ * dAGP_dx_func[i,3](x, px, y, py, β)
            du[4] -= βdot * γ * dAGP_dx_func[i,4](x, px, y, py, β)
        end
    end
end

function Deramp_EoM!(du, u, p, t)
    τ, γlist, dAGP_dx_func, βrange, Hderivs = p
    x, px, y, py = u

    β = beta(t, τ, βrange=βrange)
    βdot = betadot(t, τ, βrange=βrange)

    du[1] = Hderivs.dHdpx(x, px, y, py, β)
    du[2] = -Hderivs.dHdx(x, px, y, py, β)
    du[3] = Hderivs.dHdpy(x, px, y, py, β)
    du[4] = -Hderivs.dHdy(x, px, y, py, β)

    if !isempty(γlist)
        for i in eachindex(γlist)
            γ = γlist[i](β)
            du[1] += βdot * γ * dAGP_dx_func[i,1](x, px, y, py, β)
            du[2] -= βdot * γ * dAGP_dx_func[i,2](x, px, y, py, β)
            du[3] += βdot * γ * dAGP_dx_func[i,3](x, px, y, py, β)
            du[4] -= βdot * γ * dAGP_dx_func[i,4](x, px, y, py, β)
        end
    end
end

function EvolvedTrajectory(tau,EoM,ensemble,H_Symb; τ_drive = nothing, output_trajectory::Bool=false, saveat=nothing, γlist=[], dAGP_dx_func=[],βrange = (0, 1), batch_size::Int=8)
    tspan = (0.0, tau)
    Hderiv = build_ham_derivs(H_Symb)
    param = (τ_drive, γlist, dAGP_dx_func, βrange, Hderiv)
    n = length(ensemble)
    u0 = ensemble[1] 
    if length(ensemble) < batch_size
        batch_size = 1
    end
    prob = ODEProblem(EoM, u0, tspan, param)
    ensemble_prob = EnsembleProblem(prob, prob_func=(prob, i, repeat) -> remake(prob, u0=ensemble[i]))

    kwargs = (reltol=1e-10, abstol=1e-10, maxiters=10000000, save_everystep=output_trajectory)
    if saveat !== nothing
        sol = solve(ensemble_prob, Vern9(), EnsembleThreads(); trajectories=n, saveat=saveat, batch_size=batch_size, kwargs...)
    else
        sol = solve(ensemble_prob, Vern9(), EnsembleThreads(); trajectories=n, batch_size=batch_size, kwargs...)
    end

    return [sol[i] for i in 1:n]
end

function EnergyVariance(tau, EoM ,ensemble, H_Symb; Var_γ=[], Var_AGP=[],βrange=(0,1),τ_drive = nothing)
    τ_drive = isnothing(τ_drive) ? tau : τ_drive
    trajectories = EvolvedTrajectory(tau,EoM, ensemble,H_Symb,output_trajectory =false, γlist=Var_γ, dAGP_dx_func=Var_AGP, βrange=βrange,τ_drive=τ_drive)
    n = length(trajectories)
    energies = Vector{Float64}(undef, n)
    H = polynomial_to_function(H_Symb)
    @inbounds for i in 1:n
        sol = trajectories[i]
        state = sol.u[end]
        x, px, y, py = state
        energies[i] = H(x, px, y, py, beta(tau, τ_drive, βrange=βrange))
    end
    println("For τ", τ_drive ," Energy variance: ", var(energies))
    return var(energies)
end

function Energy(tau, EoM ,point, H_Symb; Var_γ=[], Var_AGP=[],βrange=(0,1),τ_drive = nothing)
    τ_drive = isnothing(τ_drive) ? tau : τ_drive
    trajectories = EvolvedTrajectory(tau,EoM, point,H_Symb,output_trajectory =false, γlist=Var_γ, dAGP_dx_func=Var_AGP, βrange=βrange, τ_drive=τ_drive)
    n = length(trajectories)
    energies = Vector{Float64}(undef, n)
    H = polynomial_to_function(H_Symb)
    @inbounds for i in 1:n
        sol = trajectories[i]
        state = sol.u[end]
        x, px, y, py = state
        energies[i] = H(x, px, y, py, beta(tau, τ_drive, βrange=βrange))
    end
    println("For τ", τ_drive ," Energy: ", energies[1])
    return energies[1]
end

function EnergyVarianceSequence(τ_drive, ensemble, H_symb; Var_γ=[],Var_AGP=[],βrange=(0,1), wait_time_list = [])
    println("For τ", τ_drive)
    
    ramped_ensemble = RampEnsemble(τ_drive, Ramp_EoM!, ensemble, H_symb, Var_γ = Var_γ, Var_AGP = Var_AGP, βrange = βrange)
    ramped_energy = compute_energy(ramped_ensemble, H_symb, βrange[2])
    println("  Energy variance at τ: ", var(ramped_energy))

    explored_ensemble = ExploreEnsemble(wait_time_list, Explore_EoM!, ramped_ensemble, H_symb, Var_γ = Var_γ, Var_AGP = Var_AGP, βrange = βrange)
    explored_energy = compute_energy(explored_ensemble, H_symb, βrange[2])
    println("  Energy variance at τ + wait_time: ", var(explored_energy))

    deramped_ensemble = DerampEnsemble(τ_drive, Deramp_EoM!, explored_ensemble, H_symb, Var_γ = Var_γ, Var_AGP = Var_AGP, βrange = (βrange[2], βrange[1]))
    deramped_energy = compute_energy(deramped_ensemble, H_symb, βrange[1])
    println("  Energy variance at 2τ + wait_time: ", var(deramped_energy))

    return [var(ramped_energy), var(explored_energy), var(deramped_energy)]
end

function RampEnsemble(τ_drive, EoM, ensemble, H_Symb;Var_γ=[],Var_AGP=[],βrange=(0,1))
    ramped_ensemble = []
    for (i, point) in enumerate(ensemble)
        trajectories = EvolvedTrajectory(τ_drive, EoM, [point], H_Symb, output_trajectory =false, γlist=Var_γ, dAGP_dx_func=Var_AGP, βrange=βrange, τ_drive=τ_drive)
        state = trajectories[1].u[end]
        push!(ramped_ensemble, state)
    end
    return ramped_ensemble
end

function ExploreEnsemble(wait_time_list, EoM, ensemble, H_Symb;Var_γ=[],Var_AGP=[],βrange=(0,1))
    explored_ensemble = []
    βconstant = βrange[2]
    for (i, point) in enumerate(ensemble)
        wait_time = wait_time_list[i]
        trajectories = EvolvedTrajectory(wait_time, EoM, [point], H_Symb, output_trajectory =false, γlist=Var_γ, dAGP_dx_func=Var_AGP, βrange=(βconstant,βconstant), τ_drive=wait_time)
        state = trajectories[1].u[end]
        push!(explored_ensemble, state)
    end
    return explored_ensemble
end

function DerampEnsemble(τ_drive, EoM, ensemble, H_Symb;Var_γ=[],Var_AGP=[],βrange=(1,0))
    deramped_ensemble = []
    for (i, point) in enumerate(ensemble)
        trajectories = EvolvedTrajectory(τ_drive, EoM, [point], H_Symb, output_trajectory =false, γlist=Var_γ, dAGP_dx_func=Var_AGP, βrange=βrange, τ_drive=τ_drive)
        state = trajectories[1].u[end]
        push!(deramped_ensemble, state)
    end
    return deramped_ensemble
end

function compute_energy(ensemble, H_symb, β)
    energies = Vector{Float64}(undef,length(ensemble))
    H = polynomial_to_function(H_symb)
    for (i, point) in enumerate(ensemble)
        x, px, y, py = point
        energies[i] = H(x, px, y, py, β)
    end
    return energies
end

function AGP_terms(order,H_Symb)
    terms =[]
    poly = Polynomial_Derivative(H_Symb, :β)
    for i in 1:order
        new_poly = nestedPB(poly,2*i-1)
        push!(terms, new_poly)
    end
    return terms
end

function VariationalParameter(order,ensemble,Hamil)
    n_samples = length(ensemble)
    βList = range(0,1,21)
    times = [200*asin(sqrt(2/pi*asin(sqrt(i))))/pi for i in βList]
    trajectories = EvolvedTrajectory(100.0,EoM!, ensemble, Hamil, saveat=times)

    Gfunctions = []
    base_poly = Polynomial([Term(0.5, (2,0,2,0,0))])
    Gf0 = polynomial_to_function(base_poly)
    push!(Gfunctions, (x, px, y, py, β) -> Gf0(x, px, y, py, β)) 
    for i in 1:order
        approx = nestedPB(base_poly, 2*i-1)
        pb_result = PB(approx, Hamil)
        if !(pb_result isa Polynomial)
            pb_result = Polynomial([])
        end
        Gf = polynomial_to_function(pb_result)
        G = (x, px, y, py, β) -> -1 * Gf(x, px, y, py, β)
        push!(Gfunctions, G)
    end

    n_g = length(Gfunctions)
    gamma_results = []
    for (i, β) in enumerate(βList)
        time =  200*asin(sqrt(2/pi*asin(sqrt(β))))/pi
        time_index = findfirst(t -> isapprox(t, time; atol=1e-8), times)
        Gmatrix = zeros(n_samples, n_g)
        for (sample_idx, sol) in enumerate(trajectories)
            state = sol.u[time_index]
            x, px, y, py = state
            for (g_idx, G) in enumerate(Gfunctions)
                Gmatrix[sample_idx, g_idx] = G(x, px, y, py, β)
            end
        end
        function quadform(γ)
            v = vcat(1.0, γ)
            gv = Gmatrix * v
            μ = mean(gv)
            return sum((gv .- μ).^2) / (n_samples - 1)
        end
        res = optimize(quadform, zeros(n_g-1))
        γ_opt = Optim.minimizer(res)
        γ_opt_vec = γ_opt isa Number ? [γ_opt] : γ_opt
        push!(gamma_results, γ_opt_vec)
    end
    γmat = hcat(gamma_results...)
    γ_interp = []
    for i in 1:order
        interp_obj = LinearInterpolation(βList, γmat[i, :])
        interp = β -> interp_obj(β) 
        push!(γ_interp, interp)
    end
    return γ_interp
end
#when plugging into Energy_Variance make sure AGP is alternatng
function PI_VariationalParameter(order,long_time_ensemble, H_Symb; βlist = range(0, 1, 21))
    n_samples = length(long_time_ensemble[1,:,1])  
    
    base_poly = Polynomial_Derivative(H_Symb, :β)
    Tlist = Vector{Any}(undef, 2*order+1)
    new_term = base_poly
    Tlist[1] = polynomial_to_function(new_term)

    for i in 2:(2*order+1)
        old_term  = new_term
        new_term = PB(H_Symb, old_term)
        Tlist[i] = polynomial_to_function(new_term)
    end
    ufuncs = Vector{Any}(undef, order+1) 
    vfuncs = Vector{Any}(undef,order)
    for i in 1:(order+1)
        ufuncs[i] = Tlist[2*i-1]
    end
    for i in 1:order
        vfuncs[i] = Tlist[2*i]
    end
    Tfuncs = Tlist[2:order+1]
    gamma_results = []
    for (i, β) in enumerate(βlist)
        Uvals = zeros(n_samples, order+1)
        Vvals = zeros(n_samples, order)
        Tvals = zeros(n_samples, order)
        for (j, state) in enumerate(eachrow(long_time_ensemble[i, :, :]))
            x, px, y, py = state
            for (u_idx, u) in enumerate(ufuncs)
                Uvals[j, u_idx] = u(x, px, y, py, β)
            end
            for (v_idx, v) in enumerate(vfuncs)
                Vvals[j, v_idx] = v(x, px, y, py, β)
            end
            for (t_idx, t) in enumerate(Tfuncs)
                Tvals[j, t_idx] = t(x, px, y, py, β)
            end
        end
        Σ_u = cov(Uvals)
        Σ_v = cov(Vvals)
        Σ_t = cov(Tvals)
        b_vec = -1 * Diagonal(Σ_t)
        μ = 0.00001
        aug_Σ_v = pad_matrix_with_zeros(Σ_v)
        Σ = Σ_u + μ*aug_Σ_v
        # PI_γ = inv(Σ)* b_vec
        # push!(gamma_results, PI_γ)
        # Σ = cov(Gmatrix)
        # λ = 1e-8
        # Σ_pinv = inv(Σ' * Σ + λ * I)*Σ
        function quadform(γ)
            v = vcat(1.0, γ)
            return v' * Σ * v
        end
        res = optimize(quadform, zeros(order))
        γ_opt = Optim.minimizer(res)
        γ_opt_vec = γ_opt isa Number ? [γ_opt] : γ_opt
        push!(gamma_results, γ_opt_vec)
        
    end
    γmat = hcat(gamma_results...)
    return γmat
end

function interpolater(order, γmat ;βlist = range(0,1,21))
    γ_interp = Vector{Function}(undef, order)
    for i in 1:order
        interp_obj = CubicSplineInterpolation(βlist, γmat[i, :])
        interp = β -> interp_obj(β)
        γ_interp[i] = interp
    end
    return γ_interp
end

function EnergyVariance_AGP(τlist,ensemble, order,H,H_Symb;long_time_τ=100.0,βlist = range(0,1,21))
    time = string(Dates.now())
    mkdir("results/$(time)")
    if order ==0
        return [EnergyVariance(τ, AGP_EoM!, ensemble, H) for τ in τlist]
    end
    βvals = 0:0.01:1
    agp_terms = AGP_terms(order,H_Symb)
    long_time_ensemble = Long_Time_Ensemble(ensemble, H_Symb, Long_time_τ = long_time_τ, βlist=βlist)
    γlist = PI_VariationalParameter(order, long_time_ensemble, H_Symb; βlist=βlist)
    plot()
    for (i, γ) in enumerate(γlist)
        plot!(βvals, [γ(β) for β in βvals], label="γ_$i(β)")
    end
    xlabel!("β")
    ylabel!("γ(β)")
    title!("Variational Parameters vs β")
    display(current())
    vars = [:px, :x, :py, :y]
    dAGP_dx_func = [polynomial_to_function(Polynomial_Derivative(poly, i)) for poly in agp_terms, i in vars]
    return [EnergyVariance(τ, AGP_EoM!, ensemble, H, Var_γ=γlist, Var_AGP=dAGP_dx_func) for τ in τlist]
end

function EnergyVariance_PI(τlist,ensemble, order,H_Symb;long_time_τ=100.0,βlist = range(0,1,21),time = nothing,τ_drive_list = nothing)
    variancelist= []
    βrange = [βlist[1],βlist[end]]
    if order ==0
        for (i,τ) in enumerate(τlist)
            τ_ramp = τ_drive_list[i]
            variance = EnergyVariance(τ, AGP_EoM!, ensemble, H_Symb, τ_drive=τ_ramp,βrange=βrange)
            push!(variancelist, variance)
        end
        SaveData("EnergyVariance", τ_drive_list, variancelist, order; Time = time)
        return variancelist
    end
    factor = τlist[1]/τ_drive_list[1]
    @show factor
    agp_terms = AGP_terms(order,H_Symb)
    long_time_ensemble = Long_Time_Ensemble(ensemble,H_Symb, Long_time_τ = factor*long_time_τ, βlist=βlist,factor = 1)
    γmat = PI_VariationalParameter(order, long_time_ensemble, H_Symb; βlist=βlist)
    SaveData("VariationalParameters", βlist, γmat',order,Time = time)
    γlist = interpolater(order, γmat; βlist=βlist)
    vars = [:px, :x, :py, :y]
    dAGP_dx_func = Array{Function}(undef, length(agp_terms), length(vars))
    for i in eachindex(agp_terms), j in eachindex(vars)
        dAGP_dx_func[i, j] = polynomial_to_function(Polynomial_Derivative(agp_terms[i], vars[j]))
    end
    for (i,τ) in enumerate(τlist)
        τ_ramp = τ_drive_list[i]
        variance = EnergyVariance(τ, AGP_EoM!, ensemble, H_Symb, Var_γ=γlist, Var_AGP=dAGP_dx_func, τ_drive=τ_ramp,βrange=βrange)
        push!(variancelist, variance)
    end
    SaveData("EnergyVariance", τ_drive_list, variancelist, order; Time = time)
    return variancelist
end


# -----------------------------------
# Orthonormalization of AGP functions
# -----------------------------------

# -------------------------------
# Chebychev Polynomial Basis
# -------------------------------

function AGP_Chevy_terms(order, H_Symb)
    total_terms = 2 * order
    Chebyshev_terms = Vector{Any}(undef, total_terms)
    Chebyshev_terms[1] = Polynomial_Derivative(H_Symb, :β)
    Chebyshev_terms[2] = 2*PB(H_Symb, Chebyshev_terms[1])
    for i in 3:total_terms
        Chebyshev_terms[i] = 2*PB(H_Symb, Chebyshev_terms[i-1]) - Chebyshev_terms[i-2]
    end
    return Chebyshev_terms[2:2:end]
end

function Cross_Correlation_Inner_Product(f, g, β, LongTimeEnsemble)
    x = LongTimeEnsemble[:, 1]
    px = LongTimeEnsemble[:, 2]
    y = LongTimeEnsemble[:, 3]
    py = LongTimeEnsemble[:, 4]
    f_func = isa(f, Function) ? f : polynomial_to_function(f)
    g_func = isa(g, Function) ? g : polynomial_to_function(g)
    vals_f = [f_func(x[i], px[i], y[i], py[i], β) for i in eachindex(x)]
    vals_g = [g_func(x[i], px[i], y[i], py[i], β) for i in eachindex(x)]
    innerprod = mean(vals_f .* vals_g) - mean(vals_g)*mean(vals_f)
    return innerprod
end


function normalizer(terms, Long_Time_Ensemble; βlist = range(0, 1, 21))
    term_norms = []
    for poly in terms
        normlist = Vector{Float64}(undef, length(βlist))
        for (i, β) in enumerate(βlist)
            norm = sqrt(Cross_Correlation_Inner_Product(poly, poly, β, Long_Time_Ensemble[i, :, :]))
            normlist[i] = norm
        end
        norm_interp = CubicSplineInterpolation(βlist, normlist)
        push!(term_norms, β -> norm_interp(β))
    end
    return term_norms
end

function normal_normalizer(terms, Long_Time_Ensemble; βlist = range(0, 1, 21))
    term_norms = []
    for poly in terms
        normlist = Vector{Float64}(undef, length(βlist))
        for (i, β) in enumerate(βlist)
            norm = sqrt(Inner_Product(poly, poly, β, Long_Time_Ensemble[i, :, :]))
            normlist[i] = norm
        end
        norm_interp = LinearInterpolation(βlist, normlist)
        push!(term_norms, β -> norm_interp(β))
    end
    return term_norms
end

function Chebychev_VariationalParameter(order,long_time_ensemble, H_Symb; βlist = range(0, 1, 21))
    n_samples = length(long_time_ensemble[1,:,1])  
    
    vterms = AGP_Chevy_terms(order, H_Symb)
    norms = normalizer(vterms, long_time_ensemble, βlist=βlist)
    vfuncs = [(x,px,y,py,β) -> polynomial_to_function(poly)(x,px,y,py,β)/norms[i](β) for (i,poly) in enumerate(vterms)]
    ufuncs = Vector{Any}(undef, order+1) 
    ufuncs[1] = (x,px,y,py,β) -> polynomial_to_function(Polynomial_Derivative(H_Symb, :β))(x,px,y,py,β)
    for (i,term) in enumerate(vterms)
        ufuncs[i+1] = (x,px,y,py,β) -> polynomial_to_function(-1*PB(term,H_Symb))(x,px,y,py,β)/norms[i](β)
    end
    #can check by taking inner product
    gamma_results = []
    for (i, β) in enumerate(βlist)
        Uvals = zeros(n_samples, order+1)
        Vvals = zeros(n_samples, order)
        for (j, state) in enumerate(eachrow(long_time_ensemble[i, :, :]))
            x, px, y, py = state
            for (u_idx, u) in enumerate(ufuncs)
                Uvals[j, u_idx] = u(x, px, y, py, β)
            end
            for (v_idx, v) in enumerate(vfuncs)
                Vvals[j, v_idx] = v(x, px, y, py, β)
            end
        end
        Σ_u = cov(Uvals)
        Σ_v = cov(Vvals)
        μ = 0.0001
        aug_Σ_v = pad_matrix_with_zeros(Σ_v)
        Σ = Σ_u + μ*aug_Σ_v
        # PI_γ = inv(Σ)* b_vec
        # push!(gamma_results, PI_γ)
        # Σ = cov(Gmatrix)
        # λ = 1e-8
        # Σ_pinv = inv(Σ' * Σ + λ * I)*Σ
        #function quadform(γ)
        #    v = vcat(1.0, γ)
        #    return v' * Σ * v
        #end
        #res = optimize(quadform, zeros(order))
        #γ_opt = Optim.minimizer(res)
        #γ_opt_vec = γ_opt isa Number ? [γ_opt] : γ_opt
	Σ_11 = Σ[2:end, 2:end]   # γ-γ block
        Σ_1γ = Σ[2:end, 1]       # γ-constant coupling

        γ_opt_vec = -pinv(Σ_11; atol=1e-8,rtol=1e-8) * Σ_1γ
        push!(gamma_results, γ_opt_vec)
        
    end
    γmat = hcat(gamma_results...)
    return γmat
end

function EOM_chebyshev_Constructor(order, H_Symb, long_time_ensemble; βlist = range(0, 1, 21))
    agp_terms = AGP_Chevy_terms(order, H_Symb)
    vars = [:px, :x, :py, :y]
    dAGP_dx_func = Array{Function}(undef, length(agp_terms), length(vars))
    norms  = normalizer(agp_terms, long_time_ensemble, βlist=βlist)
    for i in eachindex(agp_terms), j in eachindex(vars)
        dAGP_dx_func[i, j] = (x,px,y,py,β) -> polynomial_to_function(Polynomial_Derivative(agp_terms[i], vars[j]))(x,px,y,py,β) / norms[i](β)
    end
    return dAGP_dx_func
end

function EnergyVariance_Cheby(τlist,ensemble,long_time_ensemble, order,H_Symb;long_time_τ=100.0,βlist = range(0,1,21),wait_time_list= [],time = nothing)
    variancelist= []
    if order == 0
        for (i,τ) in enumerate(τlist)
            variance = EnergyVarianceSequence(τ, ensemble, H_Symb, βrange=(first(βlist), last(βlist)), wait_time_list = wait_time_list)
            push!(variancelist, variance)
        end
        SaveData("EnergyVariance", τlist, variancelist, order; Time = time)
        return variancelist
    end
    println("Computing Variational Parameters...")
    γmat = Chebychev_VariationalParameter(order, long_time_ensemble, H_Symb; βlist=βlist)
    SaveData("VariationalParameters", βlist, γmat',order,Time = time)
    γlist = interpolater(order, γmat; βlist=βlist)
    
    println("Computing AGP terms...")
    dAGP_dx_func = EOM_chebyshev_Constructor(order, H_Symb, long_time_ensemble; βlist=βlist)
    println("Computing Energy Variance...")
    for (i,τ) in enumerate(τlist)
        variance = EnergyVarianceSequence(τ, ensemble, H_Symb, Var_γ=γlist, Var_AGP=dAGP_dx_func, βrange=(first(βlist), last(βlist)), wait_time_list = wait_time_list)
        push!(variancelist, variance)
    end
    SaveData("EnergyVariance", τlist, variancelist, order; Time = time)
    return variancelist
end
# -------------------------------
# More Functions
# -------------------------------

function SaveData(filename, xdata, ydata, order;Time= nothing)
    if Time === nothing
        new_filename = "results/$(filename)_order$(order).txt"
    else
        new_filename = "results/$(Time)/$(filename)_order$(order).txt"
    end
    if ndims(ydata) == 1
        data = hcat(xdata, ydata)
    else
        data = hcat(xdata, eachcol(ydata)...)
    end
    writedlm(new_filename, data)
end

# ---------------------------------
# save results dyanamically
# ---------------------------------

function save_results(filename, xdata, ydata, order, τ, seed)
    if !isdir("seed_$(seed)")
        mkdir("seed_$(seed)")
    end
    new_filename = "results/seed_$(seed)/$(filename)_order$(order)_τ$(τ).txt"
    if ndims(ydata) == 1
        data = hcat(xdata, ydata)
    else
        data = hcat(xdata, eachcol(ydata)...)
    end
    writedlm(new_filename, data)
end

function combine_results(filename, order, τ, seed)
    dir = "results/seed_$(seed)"
    pattern = r"^$(filename)_order$(order)_τ([0-9eE\.\-]+)\.txt$"
    files = filter(f -> occursin(pattern, f), readdir(dir))
    combined_data = []
    for f in files
        fullpath = joinpath(dir, f)
        data = readdlm(fullpath)
        m = match(pattern, f)
        τ = m !== nothing ? parse(Float64, m.captures[1]) : missing
        if ndims(data) == 1
            push!(combined_data, hcat(fill(τ, length(data)), data))
        else
            push!(combined_data, hcat(fill(τ, size(data,1)), data))
        end
    end
    # Concatenate all data vertically
    all_data = vcat(combined_data...)
    writedlm(output_file, all_data)
    return all_data
end

function Long_Time_Ensemble(ensemble, H_Symb; Long_time_τ = 100.0, βlist = range(0, 1, 21),τ_wait = 0.0)
    τ_drive = Long_time_τ
    total_time = Long_time_τ
    βinitial = minimum(βlist)
    βfinal = maximum(βlist)
    βrange = (βinitial, βfinal)
    Saveβ = collect(βlist)
    times = [2 * τ_drive * asin(sqrt(2 / pi * asin(sqrt((i-βinitial)/(βfinal-βinitial))))) / pi for i in Saveβ]
    
    # Ensure all times are within the integration span
    times = filter(t -> t <= total_time, times)
    # Add the final time if not already included
    if !isapprox(times[end], total_time, atol=1e-10)
        push!(times, total_time)
    end

    trajectories = EvolvedTrajectory(total_time, EoM!, ensemble, H_Symb, saveat=times, βrange=βrange, τ_drive=τ_drive)
    LongTimeEnsemble = Array{Float64}(undef, length(times), length(trajectories), 4)
    
    for (time_index, time) in enumerate(times)
        for (sample_index, sol) in enumerate(trajectories)
            # The solution should now have values at all requested times
            state = sol.u[time_index]
            x, px, y, py = state
            LongTimeEnsemble[time_index, sample_index, :] = [x, px, y, py]
        end
    end
    
    return LongTimeEnsemble
end

function pad_matrix_with_zeros(A)
    n = size(A, 1)
    B = zeros(n+1, n+1)
    B[2:end, 2:end] .= A
    return B
end

function Krylov_Term_Selector(n; len::Union{Nothing,Int}=nothing)
    result = []
    function helper(current, last)
        
        if length(current) > 1 && (isnothing(len) || length(current) == len)
            push!(result, copy(current))
        end
        for nextval in (last-2):-1:1
            push!(current, nextval)
            helper(current, nextval)
            pop!(current)
        end
    end
    helper([n], n)
    return result
end



include("startup.jl")
  
function Energy_Variance_ACD()
    Job_ID = parse(Int,get(ENV, "MAIN_JOB_ID", "1"))
    order = parse(Int,get(ENV, "ORDER", "1"))
    Task_ID = parse(Int,get(ENV, "SGE_TASK_ID", "1"))
    points = parse(Int,get(ENV, "POINTS", "10"))

    dir = joinpath(homedir(),"classical-cd","reverse-cd","Julia_Files","results","Run_$(Job_ID)")
    #dir = joinpath(homedir(),"Documents","College","Research","Counter Diabatic Driving","2D simple case","Julia Files","Energy Variance","results","Run_1")


    @load joinpath(dir,"taumin.jld2") taumin
    @load joinpath(dir,"taumax.jld2") taumax
    @load joinpath(dir,"Ensemble.jld2") ensemble
    @load joinpath(dir,"Long_Time_Ensemble.jld2") long_time_ensemble
    @load joinpath(dir,"βlist.jld2") βlist
    @load joinpath(dir,"Hamiltonian.jld2") Hamiltonian
    @load joinpath(dir,"waittimes.jld2") wait_time_list

    if points == 1
        τ = taumin
    else 
        τ = 10.0^((taumax-taumin)/(points - 1)*(Task_ID - 1) + taumin)
    end
    if log10(τ) > taumax || log10(τ) < taumin
        @warn "Task_ID $(Task_ID) out of range for τ: $(τ)"
        return
    end

    if order == 0
       Evar = EnergyVarianceSequence(τ, ensemble, Hamiltonian, βrange=(first(βlist), last(βlist)), wait_time_list = wait_time_list)
    else
        @load joinpath(dir,"Variational_Parameters_order_$(order).jld2") γmat
        
        γlist = interpolater(order, γmat; βlist=βlist)

        dAGP_dx_func = EOM_chebyshev_Constructor(order, Hamiltonian, long_time_ensemble, βlist=βlist)
        Evar = EnergyVarianceSequence(τ, ensemble, Hamiltonian, Var_γ=γlist, Var_AGP=dAGP_dx_func, βrange=(first(βlist), last(βlist)), wait_time_list = wait_time_list)
    end
    
    
    @save joinpath(dir,"Energy_Variance_order_$(order)_$(Task_ID).jld2") τ Evar
end

Energy_Variance_ACD()


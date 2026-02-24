include("startup.jl")

function compute_Var_Param_Points()
    Job_ID = parse(Int,get(ENV, "MAIN_JOB_ID", "1"))
    order = parse(Int,get(ENV, "ORDER", "1"))
    if order == 0 
        return println("nah")  
    end
    dir = joinpath(homedir(),"classical-cd","reverse-cd","Julia_Files","results","Run_$(Job_ID)")
    @load joinpath(dir,"Long_Time_Ensemble.jld2") long_time_ensemble
    @load joinpath(dir,"βlist.jld2") βlist
    @load joinpath(dir,"Hamiltonian.jld2") Hamiltonian
    #dir = joinpath(homedir(),"Documents","College","Research","Counter Diabatic Driving","2D simple case","Julia Files","Energy Variance","results","Run_1")
    γmat = Chebychev_VariationalParameter(order, long_time_ensemble, Hamiltonian; βlist=βlist)
    @save joinpath(dir,"Variational_Parameters_order_$(order).jld2") γmat
end

compute_Var_Param_Points()


include("startup.jl")
  
function Energy_Variance_ACD()
    Job_ID = parse(Int, get(ENV, "MAIN_JOB_ID", "1"))
    Task_ID = parse(Int, get(ENV, "TASK_ID", "1"))
    
    # Get individual tau and beta values from environment
    τ = parse(Float64, get(ENV, "TAU_VALUE", "1.0"))
    β_final = parse(Float64, get(ENV, "BETA_VALUE", "1.0"))

    dir = joinpath(homedir(), "classical-cd", "bear-drive", "Julia_Files", "results", "Run_$(Job_ID)")

    @load joinpath(dir, "Ensemble.jld2") ensemble
    @load joinpath(dir, "Hamiltonian.jld2") Hamiltonian
    @load joinpath(dir, "waittimes.jld2") wait_time_list

    # Use the specific beta value for this job
    βlist = range(0, β_final, 101)

    # Since you're only doing order 0
    Evar = EnergyVarianceSequence(τ, ensemble, Hamiltonian, 
                                βrange=(0.0, β_final), 
                                wait_time_list = wait_time_list)
    
    # Save with tau and beta info in filename
    @save joinpath(dir, "Energy_Variance_tau_$(τ)_beta_$(β_final).jld2") τ β_final Evar
end

Energy_Variance_ACD()


 
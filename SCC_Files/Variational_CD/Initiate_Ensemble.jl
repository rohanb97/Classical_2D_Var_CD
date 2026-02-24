include("startup.jl")

function Initialize_System()
    Job_ID = parse(Int,get(ENV, "MAIN_JOB_ID", "1"))
    seed = 1234;
    Population = parse(Int,get(ENV, "SAMPLE_POINTS", "1000"));
    AGP_Population = 10000
    βlist = range(0,1,21);
    E0 = 1.0;
    long_time_τ = 100.0;
    taumin = -3
    taumax = 0
    Hamiltonian = Polynomial([
        Term(0.5, (2, 0, 0, 0, 0)),
        Term(0.5, (0, 2, 0, 0, 0)),
        Term(0.5, (0, 0, 2, 0, 0)),
        Term(0.5, (0, 0, 0, 2, 0)),
        Term(0.5, (2, 0, 2, 0, 1))
    ])

    Random.seed!(seed) 
    H = polynomial_to_function(Hamiltonian)
    V(x,y,β) = H(x,0,y,0,β)
    ensemble = MicroSample(Population,V,E0; βinitial = βlist[1])
    AGP_ensemble = MicroSample(AGP_Population,V,E0; βinitial = βlist[1])
    wait_time_list = [rand()*(70-30)+30 for i in eachindex(ensemble)]
    long_time_ensemble = Long_Time_Ensemble(AGP_ensemble,Hamiltonian, Long_time_τ = long_time_τ, βlist=βlist);
    new_dir = joinpath(homedir(),"classical-cd","reverse-cd","Julia_Files","results","Run_$(Job_ID)")
    #new_dir = joinpath(homedir(),"Documents","College","Research","Counter Diabatic Driving","2D simple case","Julia Files","Energy Variance","results","Run_1")
    if !isdir(new_dir)
        mkdir(new_dir)
    end
    @save joinpath(new_dir,"taumin.jld2") taumin
    @save joinpath(new_dir,"taumax.jld2") taumax
    @save joinpath(new_dir,"Long_Time_Ensemble.jld2") long_time_ensemble
    @save joinpath(new_dir,"Ensemble.jld2") ensemble
    @save joinpath(new_dir,"Hamiltonian.jld2") Hamiltonian
    @save joinpath(new_dir,"βlist.jld2") βlist
    @save joinpath(new_dir,"waittimes.jld2") wait_time_list
    writedlm(joinpath(new_dir,"seed.txt"), [seed])

end

Initialize_System()
include("startup.jl")

function Initialize_System()
    Job_ID = parse(Int, get(ENV, "MAIN_JOB_ID", "1"))
    seed = 1234;
    Population = parse(Int, get(ENV, "SAMPLE_POINTS", "1000"));
    
    # Read tau and beta lists from text file in home directory
    new_dir = joinpath(homedir(), "classical-cd", "bear-drive", "Julia_Files", "results", "Run_$(Job_ID)")
    βlist = range(0, 1, 101);
    E0 = 1.0;
    
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
    ensemble = MicroSample(Population, V, E0; βinitial = 0)
    wait_time_list = [0 for i in eachindex(ensemble)]
    
    if !isdir(new_dir)
        mkdir(new_dir)
    end
    
    @save joinpath(new_dir, "Ensemble.jld2") ensemble
    @save joinpath(new_dir, "Hamiltonian.jld2") Hamiltonian
    @save joinpath(new_dir, "βlist.jld2") βlist
    @save joinpath(new_dir, "waittimes.jld2") wait_time_list
    writedlm(joinpath(new_dir, "seed.txt"), [seed])
end

Initialize_System()
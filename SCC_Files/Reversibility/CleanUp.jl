include("startup.jl")

function presentable()
    Job_ID = parse(Int, get(ENV, "MAIN_JOB_ID", "1"))
    
    # Read tau and beta lists from text file in home directory
    read_dir = joinpath(homedir(), "classical-cd", "bear-drive", "Julia_Files")
    config_file = joinpath(read_dir, "tau_beta_config.txt")
    config_lines = readlines(config_file)
    
    # Find lines that aren't comments (don't start with #)
    data_lines = filter(line -> !startswith(strip(line), "#") && !isempty(strip(line)), config_lines)
    
    # Parse tau and beta lists
    taulist_str = strip(data_lines[1])
    betalist_str = strip(data_lines[2])
    
    taulist = [parse(Float64, strip(x)) for x in split(taulist_str, ",")]
    betalist = [parse(Float64, strip(x)) for x in split(betalist_str, ",")]
    
    dir = joinpath(homedir(), "classical-cd", "bear-drive", "Julia_Files", "results", "Run_$(Job_ID)")
    
    Variance_List = []
    
    # Loop through all tau-beta combinations
    for tau_val in taulist
        for beta_val in betalist
            filename = "Energy_Variance_tau_$(tau_val)_beta_$(beta_val).jld2"
            filepath = joinpath(dir, filename)
            
            if isfile(filepath)
                @load filepath τ β_final Evar
                push!(Variance_List, (τ, β_final, Evar))
                rm(filepath)  # Remove individual files after combining
            else
                println("File $filename does not exist.")
            end
        end
    end
    
    # Write combined results
    writedlm(joinpath(dir, "Energy_Variance_combined.txt"), Variance_List, ',')
    println("Collected $(length(Variance_List)) files and wrote to Energy_Variance_combined.txt")
end

presentable()
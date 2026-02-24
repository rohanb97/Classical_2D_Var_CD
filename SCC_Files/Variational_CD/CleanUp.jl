include("startup.jl")


function presentable()
    Job_ID = parse(Int,get(ENV, "MAIN_JOB_ID", "1"))
    order = parse(Int,get(ENV, "ORDER", "0"))
    points = parse(Int,get(ENV, "POINTS", "10"))
    Variance_List = []
    for i in 1:points
        if isfile("results/Run_$(Job_ID)/Energy_Variance_order_$(order)_$(i).jld2")
            @load "results/Run_$(Job_ID)/Energy_Variance_order_$(order)_$(i).jld2" τ Evar
            push!(Variance_List, (τ, Evar))
            rm("results/Run_$(Job_ID)/Energy_Variance_order_$(order)_$(i).jld2")
        else
            println("File results/Run_$(Job_ID)/Energy_Variance_order_$(order)_$(i).jld2 does not exist.")
        end
    end
    writedlm("results/Run_$(Job_ID)/Energy_Variance_order_$(order).txt", Variance_List, ',')
end

presentable()

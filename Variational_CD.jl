using DifferentialEquations, OrdinaryDiffEq, LinearAlgebra, Optim
using Statistics
using Interpolations
using StaticArrays
using DelimitedFiles
using Dates
using Random
using Profile
using ProfileView
using JLD2
using Plots

include("PolynomialClass.jl")
include("EvarFunctions.jl")

function Variational_CD(taulist,order,Time)
    Population = 1000;
    E0 = 1.0;

    Hamiltonian = Polynomial([
        Term(0.5, (2, 0, 0, 0, 0)),
        Term(0.5, (0, 2, 0, 0, 0)),
        Term(0.5, (0, 0, 2, 0, 0)),
        Term(0.5, (0, 0, 0, 2, 0)),
        Term(0.5, (2, 0, 2, 0, 1)),
        Term(0.25, (4, 0, 0, 0, 1)),
        Term(0.25, (0, 0, 4, 0, 1))
    ])

    βrange = range(0,0.229,21)
    wait_times = [rand()*(70-30) + 30 for i in 1:Population]
    H = polynomial_to_function(Hamiltonian)
    V(x,y,β) = H(x,0,y,0,β);
    ensemble = MicroSample(Population,V,E0,βinitial = first(βrange));
    long_time_ensemble = Long_Time_Ensemble(ensemble,Hamiltonian, Long_time_τ = 100, βlist=βrange)
    
    [EnergyVariance_Cheby(taulist,ensemble,long_time_ensemble,order,Hamiltonian; βlist = βrange, time = Time,wait_time_list = wait_times) for order in 0:(order)]
end

function plotdata(taulist,data, Time,order)
    p1 = plot()
    p2 = plot()
    p3 = plot()
    p = [p1,p2,p3]
    rampdata = [[ log10(set[i][1]) for i in eachindex(taulist)] for set in data]
    exploredata = [[ log10(set[i][2]) for i in eachindex(taulist)] for set in data]
    derampdata = [[ log10(set[i][3]) for i in eachindex(taulist)] for set in data]
    logdata =[rampdata,exploredata,derampdata]
    plotnames = ["rampplot","exploreplot","derampplot"]
    for i in 1:3
        for j in 1:(order+1)
            plot!(p[i], log10.(taulist), logdata[i][j], label="Order $(j-1)",marker=:o)
        end
        xlabel!(p[i],"log(τ)")
        ylabel!(p[i],"log(Variance)")
        title!(p[i],"Energy Variance vs τ ($(plotnames[i]))")
        savefig(p[i],"results/$(Time)/$(plotnames[i]).png")
        display(p[i])
    end
end


Time = string(Dates.now())  
mkdir("results/$(Time)")
taumin = -3;
taumax = 2;
points = 3;
taulist = 10 .^ range(taumin, taumax, length=points); 
order = 1;
@time data = Variational_CD(taulist,order,Time)
plotdata(taulist,data, Time, order)




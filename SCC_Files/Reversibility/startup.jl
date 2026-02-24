using Pkg

Pkg.activate(joinpath(homedir(),"classical-cd","reverse-cd","Julia_Files","2D_Energy_Variance"))


using DifferentialEquations, OrdinaryDiffEq, LinearAlgebra, Optim
using Statistics
using Interpolations
using StaticArrays
using DelimitedFiles
using Dates
using Random
using JLD2

include(joinpath(homedir(),"classical-cd","reverse-cd","Julia_Files","EvarFunctions.jl"))
include(joinpath(homedir(),"classical-cd","reverse-cd","Julia_Files","PolynomialClass.jl"))
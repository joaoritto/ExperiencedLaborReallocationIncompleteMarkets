# Skilled Labor Reallocation in an Incomplete Markets Economy

path="C:\\Users\\joaor\\Dropbox\\Economics\\ThirdYearPaper\\Code\\ExperiencedLaborReallocationIncompleteMarkets\\"

cd(path)

include(path*"ValueFunctionsIteration.jl")
include(path*"StationaryEquilibrium.jl")
include(path*"AnalyzingResults.jl")

using Statistics,LinearAlgebra,Plots,SparseArrays,Interpolations,Optim

# Parameters
β=1/1.03 # Discount factor
σ=2.0  # Inverse IES
ρ=0.04 # Exogenous separation
δ=0.025 # Separation with loss of skill
α=0.08 # Probability of becoming skilled
b=0.025 # Unemployment benefits; b > -̲a*r or c<0 at lowest wealth
σ_ϵ=0.4 # s.d. of taste shocks
ξ=0.5 # Unemployed share in matching technology
m=0.45 # Productivity of matching technology
κ=0.01 # Vacancy cost
γ=0.3 # Share of inexperienced workers
ν=1.5 # Elasticity of substitution between intermediate goods
a_min=-0.5
# Prices
w=[0.1 1.2; 0.1 1.2]
r=0.015

z=[1.0; 1.0]

# Grids
n_i=2
n_s=2
n_a=100
n_μ=30

grid_i=1:n_i
grid_s=1:n_s
grid_a=LinRange(a_min,50,n_a)
grid_μ=LinRange(0.7,1-1e-2,n_μ)


grids=(grid_i,grid_s,grid_a,grid_μ)


V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ=ValueFunctions(grids,w)

pol_val_functions=(V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)

Φ=ComputeDistribution(grids,pol_val_functions)

Y,E_I,E_E,U_I,U_E=ComputeAggregates(grids,pol_val_functions,Φ,z)

#wages(Y,z,i,e)=((1/n_i)*γ*Y^(1/ν)*z^(1-(1/ν))*i^(γ*(1-(1/ν))-1)*e^((1-γ)*(1-(1/ν))),
#                    (1/n_i)*(1-γ)*Y^(1/ν)*z^(1-(1/ν))*i^(γ*(1-(1/ν)))*e^((1-γ)*(1-(1/ν))-1))


#pol_val_functions,Φ,Y,I,E,U_I,U_E=GeneralEquilibrium(grids,z)

#=
V_E,V_U,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ=ValueFunctions(para,grids;Guess=pol_val_functions)
=#

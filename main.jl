# Skilled Labor Reallocation in an Incomplete Markets Economy

path="C:\\Users\\joaor\\Dropbox\\Economics\\ThirdYearPaper\\Code\\ExperiencedLaborReallocationIncompleteMarkets\\"

cd(path)

include(path*"ValueFunctionsIteration.jl")
include(path*"StationaryEquilibrium.jl")
include(path*"AnalyzingResults.jl")

using Statistics,LinearAlgebra,Plots,SparseArrays

# Parameters
β=0.97 # Discount factor
σ=2.0  # Inverse IES
ρ=0.05 # Exogenous separation
δ=0.01 # Separation with loss of skill
α=0.02 # Probability of becoming skilled
b=0.02 # Unemployment benefits; b > -̲a*r or c<0 at lowest wealth
σ_ϵ=0.2 # s.d. of taste shocks
ξ=0.5 # Unemployed share in matching technology
m=0.45 # Productivity of matching technology
κ=0.07 # Vacancy cost
γ=0.33 # Share of inexperienced workers
ν=1.5 # Elasticity of substitution between intermediate goods

# Prices
w=[0.1 0.2; 0.12 0.24]
r=0.02

z=[3.0; 4.0]

# Grids
n_i=2
n_s=2
n_a=120
n_μ=120

grid_i=1:n_i
grid_s=1:n_s
grid_a=LinRange(-0.2,4,n_a)
grid_μ=LinRange(0.5,1-1e-2,n_μ)


grids=(grid_i,grid_s,grid_a,grid_μ)


V_E,V_U,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ=ValueFunctions(grids)

pol_val_functions=(V_E,V_U,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)

Φ=ComputeDistribution(grids,pol_val_functions)

Y,E_I,E_E,U_I,U_E=ComputeAggregates(grids,pol_val_functions,Φ,z)

wages(Y,z,i,e)=((1/n_i)*γ*Y^(1/ν)*z^(1-(1/ν))*i^(γ*(1-(1/ν))-1)*e^((1-γ)*(1-(1/ν))),
                    (1/n_i)*(1-γ)*Y^(1/ν)*z^(1-(1/ν))*i^(γ*(1-(1/ν)))*e^((1-γ)*(1-(1/ν))-1))

#=
V_E,V_U,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ=ValueFunctions(para,grids;Guess=pol_val_functions)
=#

# Skilled Labor Reallocation in an Incomplete Markets Economy

path="C:\\Users\\joaor\\Dropbox\\Economics\\ThirdYearPaper\\Code\\ExperiencedLaborReallocationIncompleteMarkets\\"

cd(path)

include(path*"ValueFunctionsIteration2.jl")
include(path*"StationaryEquilibrium2.jl")
include(path*"AnalyzingResults.jl")

using Statistics,LinearAlgebra,Plots,SparseArrays,Interpolations,Optim

# Parameters
β=1/1.02 # Discount factor
σ=2.0  # Inverse IES
ρ=0.04 # Exogenous separation
δ=0.01 # Separation with loss of skill
α=0.05 # Probability of becoming skilled
b=0.05 # Unemployment benefits; b > -̲a*r or c<0 at lowest wealth
σ_ϵ=0.5 # s.d. of taste shocks
ξ=0.5 # Unemployed share in matching technology
m=0.45 # Productivity of matching technology
κ=0.06 # Vacancy cost
γ=0.42 # Share of inexperienced workers
ν=1.5 # Elasticity of substitution between intermediate goods
a_min=-0.5
a_max=80
# Prices
w=[0.1 1.3; 0.1 1.3]
r=0.02

z=[2.2; 2.2]

# Grids
n_i=2
n_s=2
n_μ=30

grid_i=1:n_i
grid_s=1:n_s
grid_μ=LinRange(0.7,1-1e-2,n_μ)

n_a=20
using_multigrid=1 # If set to 0, code runs just once with n_a grid points. If set to 1 it starts with n_a and increases the grid
nGrids_a=[n_a,50,100]

Results=multigrid(nGrids_a)

Check=(maximum(pol_a_E)<a_max)

(V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=Results

grid_a=LinRange(a_min,a_max,nGrids_a[end])
grids=(grid_i,grid_s,grid_a,grid_μ)
pol_a_Ei,pol_a_Ui=transformPola(pol_a_E,pol_a_U,grids)

pol_functions=(pol_a_Ei,pol_a_Ui,pol_μ_U,pol_σ_E,pol_σ_U)
Φ=ComputeDistribution(grids,pol_functions)

Y,E_I,E_E,U_I,U_E=ComputeAggregates(grids,pol_functions,Φ,z)

#wages(Y,z,i,e)=((1/n_i)*γ*Y^(1/ν)*z^(1-(1/ν))*i^(γ*(1-(1/ν))-1)*e^((1-γ)*(1-(1/ν))),
#                    (1/n_i)*(1-γ)*Y^(1/ν)*z^(1-(1/ν))*i^(γ*(1-(1/ν)))*e^((1-γ)*(1-(1/ν))-1))


#pol_val_functions,Φ,Y,I,E,U_I,U_E=GeneralEquilibrium(grids,z)

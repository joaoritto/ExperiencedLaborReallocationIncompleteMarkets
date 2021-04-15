# Skilled Labor Reallocation in an Incomplete Markets Economy

path="C:\\Users\\joaor\\Dropbox\\Economics\\ThirdYearPaper\\Code\\ExperiencedLaborReallocationIncompleteMarkets\\"

cd(path)

include(path*"ValueFunctionsIteration2.jl")
include(path*"StationaryEquilibrium2.jl")
include(path*"AnalyzingResults.jl")

using Statistics,LinearAlgebra,Plots,SparseArrays,Interpolations,Optim

# Choices: i) Partial equilibrium or General Equilibrium, ii) Use multigrid?

PE=1 # If set to 0, code runs the GE, if set to 1 it runs the PE
using_multigrid=1 # If set to 0, code runs just once with n_a grid points. If set to 1 it starts with n_a and increases the grid

# Calibration of 2 months

# Parameters
β=1/1.03 #0.9935 # Discount factor
σ=2.0  # Inverse IES
ρ=0.04#/6 # Exogenous separation
δ=0.01#/6 # Separation with loss of skill
α=0.08#/6 # Probability of becoming skilled
b=0.025 # Unemployment benefits; b > -̲a*r or c<0 at lowest wealth - Calibrate for ratio to wage
σ_ϵ=0.4 # s.d. of taste shocks
ξ=0.5 # Unemployed share in matching technology
m=0.48 # Productivity of matching technology
κ=0.01 # Vacancy cost - Calibrate to get unemployment rate
γ=0.3 # Productivity share of inexperienced workers
ν=1.5 # Elasticity of substitution between intermediate goods

a_min=-0.5
a_max=55

# Prices
if PE==1
    w=[0.06 0.7; 0.06 0.7]
end
r=0.015#/6

z=[1.0; 1.0]

# Grids
n_i=2
n_s=2
n_μ=30

grid_i=1:n_i
grid_s=1:n_s
grid_μ=LinRange(0.7,1-1e-2,n_μ)

n_a=100
nGrids_a=[n_a]
n_anew=100 #nGrids_a[end]


if PE==1
    if using_multigrid==1
        pol_val_functions=multigrid(nGrids_a,w)
        (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions
    elseif using_multigrid==0
        grid_a=LinRange(a_min,a_max,n_a)
        grids=(grid_i,grid_s,grid_a,grid_μ)
        V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ=ValueFunctions(grids,w)
        pol_val_functions=(V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)
    end

    Check=(maximum(pol_a_E)<a_max)

    grid_a=LinRange(a_min,a_max,nGrids_a[end])
    grids=(grid_i,grid_s,grid_a,grid_μ)
    if n_anew!=nGrids_a[end]
        grid_a=LinRange(a_min,a_max,n_anew)
        pol_val_functions=transformVPolFunctions(pol_val_functions,grids,grid_a)
        (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions
        grids=(grid_i,grid_s,grid_a,grid_μ)
    end

    pol_a_Ei,pol_a_Ui=transformPola(pol_a_E,pol_a_U,grids)

    pol_functions=(pol_a_Ei,pol_a_Ui,pol_μ_U,pol_σ_E,pol_σ_U,θ)
    Φ=ComputeDistribution(grids,pol_functions)

    Y,E_I,E_E,U_I,U_E=ComputeAggregates(grids,pol_functions,Φ,z)

    wages(Y,z,i,e)=((1/n_i)*γ*Y^(1/ν)*z^(1-(1/ν))*i^(γ*(1-(1/ν))-1)*e^((1-γ)*(1-(1/ν))),
                        (1/n_i)*(1-γ)*Y^(1/ν)*z^(1-(1/ν))*i^(γ*(1-(1/ν)))*e^((1-γ)*(1-(1/ν))-1))

    wages_result=wages(Y,z[1],E_I[1],E_E[1])
    display(wages_result)
    display(E_I[1])
    display(E_E[1])

elseif PE==0
    pol_val_functions,Φ,Y,E_I,E_E,U_I,U_E=GeneralEquilibrium(z)
end

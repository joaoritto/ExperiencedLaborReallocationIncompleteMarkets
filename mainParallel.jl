# Skilled Labor Reallocation in an Incomplete Markets Economy

using Distributed, SharedArrays

addprocs(6-length(procs()))

path="C:\\Users\\joaor\\Dropbox\\Economics\\ThirdYearPaper\\Code\\ExperiencedLaborReallocationIncompleteMarkets\\"
cd(path)

include(path*"ValueFunctionsIterationParallel.jl")
include(path*"StationaryEquilibriumParallel.jl")
include(path*"TransitionParallel.jl")
include(path*"AnalyzingResults.jl")

@everywhere using Statistics,LinearAlgebra,Plots,SparseArrays,Interpolations,Optim

# Choices: i) Partial equilibrium or General Equilibrium, ii) Use multigrid?

PE=0 # If set to 0, code runs the GE, if set to 1 it runs the PE
using_multigrid=1 # If set to 0, code runs just once with n_a grid points. If set to 1 it starts with n_a and increases the grid

# Calibration of 2 months

# Parameters
@everywhere β=0.9935 # Discount factor
@everywhere σ=1.2  # Inverse IES
@everywhere ρ=0.04/6 # Exogenous separation
@everywhere δ=0.025/6 # Separation with loss of skill
@everywhere α=0.08/6 # Probability of becoming skilled
@everywhere b=0.025 # Unemployment benefits; b > -̲a*r or c<0 at lowest wealth - Calibrate for ratio to wage
@everywhere σ_ϵ=0.4 # s.d. of taste shocks
@everywhere ξ=0.5 # Unemployed share in matching technology
@everywhere m=0.48 # Productivity of matching technology
@everywhere κ=0.01 # Vacancy cost - Calibrate to get unemployment rate
@everywhere γ=0.3 # Productivity share of inexperienced workers
@everywhere ν=1.5 # Elasticity of substitution between intermediate goods

@everywhere a_min=-0.5
@everywhere a_max=30

# Prices
if PE==1
    @everywhere w=[0.0466926 0.577147; 0.0466926 0.577147]
end
@everywhere r=0.015/6

@everywhere z=[1.0; 1.0]

# Grids
@everywhere n_i=2
@everywhere n_s=2
@everywhere n_μ=100

@everywhere grid_i=1:n_i
@everywhere grid_s=1:n_s
@everywhere grid_μ=LinRange(0.7,1-1e-2,n_μ)

n_a=50
nGrids_a=[n_a,100,200]
n_anew=800 #nGrids_a[end]


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
    (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions
end

zt=zeros(2,16)
zt[:,1]=[0.8;1.0]
for t in 2:15
    zt[:,t]=0.2*[1;1]+0.8*zt[:,t-1]
end
zt[:,16]=[1.0;1.0]

# WANT THE TRANSITION TO EXPLORE INTERPOLATION SO THAT ONE CAN USE LOWER NUMBER OF GRID POINTS

#StatEq=(V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ,Φ,Y,E_I,E_E,U_I,U_E)
#NewI,NewE,NewU_I,NewU_E=Transition(grids,StatEq,zt;Guess=false)
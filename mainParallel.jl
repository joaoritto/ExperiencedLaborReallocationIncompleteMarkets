# Skilled Labor Reallocation in an Incomplete Markets Economy

using Distributed, SharedArrays

addprocs(6-length(procs()))

path="C:\\Users\\joaor\\Dropbox\\Economics\\ThirdYearPaper\\Code\\ExperiencedLaborReallocationIncompleteMarkets\\"
cd(path)

include(path*"ValueFunctionsIterationParallel.jl")
include(path*"StationaryEquilibriumParallel.jl")
include(path*"TransitionParallel2.jl")
include(path*"AnalyzingResults.jl")

@everywhere using Statistics,LinearAlgebra,Plots,SparseArrays,Interpolations,Optim,StatsBase

# Choices: i) Partial equilibrium or General Equilibrium, ii) Use multigrid?

PE=0 # If set to 0, code runs the GE, if set to 1 it runs the PE
small_grid=0
comp_transition=1
using_multigrid=1 # If set to 0, code runs just once with n_a grid points. If set to 1 it starts with n_a and increases the grid

# Calibration (1 period=2 months)

# Parameters
@everywhere β=0.9935 # Discount factor
@everywhere σ=2.0  # Inverse IES
@everywhere ρ=0.032 # Exogenous separation
@everywhere δ=0.005 # Separation with loss of skill
@everywhere α=0.125/6 # Probability of becoming skilled
@everywhere b=0.45*0.94 # Replacement rate: Unemployment benefits; b*w > -̲a*r or c<0 at lowest wealth
@everywhere σ_ϵ=0.1 # s.d. of taste shocks
@everywhere ξ=0.5 # Unemployed share in matching technology
@everywhere m=0.48 # Productivity of matching technology
@everywhere κ_e=0.02 # Vacancy cost experienced
@everywhere κ_i=0.02 # Vacancy cost inexperienced
@everywhere γ=0.18 # Productivity share of inexperienced workers
@everywhere ν=10.0 # Elasticity of substitution between intermediate goods
@everywhere F_e=0.95 # Fixed cost of hiring experienced
@everywhere F_i=0.95 # Fixed cost of hiring inexperienced

@everywhere a_min=0.0
@everywhere a_max=5.0

# Prices
if PE==1
    @everywhere w=[0.54018 0.6442396; 0.54018 0.6442396]
end
@everywhere r=0.015/6

@everywhere z=[2.0; 2.0]

# Grids
@everywhere n_i=2
@everywhere n_s=2
@everywhere n_μ=100

@everywhere grid_i=1:n_i
@everywhere grid_s=1:n_s
@everywhere grid_μ=LinRange(0.6,1.0-1e-2,n_μ)

if small_grid==1
    n_a=20
    nGrids_a=[n_a,30,50,100]
else
    n_a=30
    nGrids_a=[n_a,50,100,200]
end

n_anew=1500 #nGrids_a[end]


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

# Computing unemployment duration
# Simulation??

if comp_transition==1

    ###################################################
    #           Transition Exercises
    ###################################################

    # 1) Transitory shock, 3 periods (6 months)
    # 2) Persistent shock, 9 periods (1.5 years)
    # 3) Permanent shock, new stationary eq

    shockdur=9
    zt=z*ones(1,shockdur+1)
    zt[:,1]=[1.6;2.0]
    for t in 2:shockdur
        zt[:,t]=zt[:,1]
    end


    grid_a=LinRange(a_min,a_max,n_anew)
    grids=(grid_i,grid_s,grid_a,grid_μ)

    StatEq=(V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ,Φ,Y,E_I,E_E,U_I,U_E)
    NewI,NewE,NewU_I,NewU_E,pol_val_results=Transition(grids,StatEq,zt)


    aggregates_transition=(NewI,NewE,NewU_I,NewU_E)


    PlotResultsStatEq(grids,StatEq)
    PlotResultsTransition(grids,zt,pol_val_results,aggregates_transition)

    (V_E_Tr,V_U_Tr,W_E_Tr,W_U_Tr,pol_a_Ei_Tr,pol_a_Ui_Tr,pol_μ_U_Tr,pol_σ_E_Tr,pol_σ_U_Tr,J_Tr,θ_Tr,Φ_Tr)=pol_val_results

    plot(pol_σ_U_Tr[1501:3000,1,1:5])
    plot!(pol_σ_U_Tr[1501:3000,1,end])

    plot(p.(θ_Tr[150000+501:150000+600,1:3]),legend=false)
    plot!(p.(θ_Tr[150000+501:150000+600,end]),legend=false)
end

# Skilled Labor Reallocation in an Incomplete Markets Economy

using Distributed, SharedArrays

addprocs(6-length(procs()))

path="C:\\Users\\joaor\\Dropbox\\Economics\\ThirdYearPaper\\Code\\ExperiencedLaborReallocationIncompleteMarkets\\"
cd(path)

include(path*"NewValueFunctionsIterationParallel.jl")
include(path*"NewStationaryEquilibriumParallel.jl")
include(path*"NewTransitionParallel.jl")
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
@everywhere δ=1/(45*6) # Death
@everywhere α=[1/12;1/18;1/18;0] # Probability of being promoted
@everywhere χ=[0;0.2;0.2;0.2] # Probability of losing skill
@everywhere b=0.45*0.94 # Replacement rate: Unemployment benefits; b*w > -̲a*r or c<0 at lowest wealth
@everywhere σ_ϵ=0.1 # s.d. of taste shocks
@everywhere ξ=0.5 # Unemployed share in matching technology
@everywhere m=0.48 # Productivity of matching technology
@everywhere κ=[0.28;0.28;0.28;0.28] # Vacancy cost
@everywhere γ=[0.055563;0.098816;0.159657;0.685964] # Productivity share workers
@everywhere ν=10.0 # Elasticity of substitution between intermediate goods


@everywhere a_min=0.0
@everywhere a_max=5.0

# Prices
if PE==1
    @everywhere w=[0.342383 0.360871 0.38348 0.400614; 0.342383 0.360871 0.38348 0.400614]
end
@everywhere r=0.015/6

@everywhere z=2.0*ones(2) #n_i=2

# Grids
@everywhere n_i=2
@everywhere n_s=4
@everywhere n_μ=1 #100

@everywhere grid_i=1:n_i
@everywhere grid_s=1:n_s
@everywhere grid_μ=0.94 # LinRange(0.6,1.0-1e-2,n_μ)

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

    Y,E,U=ComputeAggregates(grids,pol_functions,Φ,z)

    wages(Y,z,E,s_i)=(1/n_i)*γ[s_i]*Y^(1/ν)*z^(1-(1/ν))*E[s_i]^(γ[s_i]*(1-(1/ν))-1)*prod(E[1:end .!=s_i].^((γ[1:end .!=s_i])*(1-(1/ν))))

    wages_result=zeros(n_s)
    for s_i in 1:n_s
        wages_result[s_i]=wages(Y,z[1],E[1,:],s_i)
    end
    display(wages_result)
    display(E[1,:])

elseif PE==0
    pol_val_functions,Φ,Y,E,U=GeneralEquilibrium(z)
    (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions
end


if comp_transition==1

    ###################################################
    #           Transition Exercises
    ###################################################

    # 1) Transitory shock, 3 periods (6 months)
    # 2) Persistent shock, 9 periods (1.5 years)
    # 3) Permanent shock, new stationary eq

    shockdur=40
    zt=z*ones(1,shockdur+1)
    zt[:,1]=[1.7;2.0]
    for t in 2:shockdur
        zt[:,t]=0.9*zt[:,t-1]+0.1*z
    end


    grid_a=LinRange(a_min,a_max,n_anew)
    grids=(grid_i,grid_s,grid_a,grid_μ)

    StatEq=(V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ,Φ,Y,E,U)
    NewE,NewU,pol_val_results=Transition(grids,StatEq,zt;Guess=NewE)


    aggregates_transition=(NewE,NewU)

    PlotResultsStatEq(grids,StatEq)
    #PlotResultsTransition(grids,zt,pol_val_results,aggregates_transition)

    (V_E_Tr,V_U_Tr,W_E_Tr,W_U_Tr,pol_a_Ei_Tr,pol_a_Ui_Tr,pol_μ_U_Tr,pol_σ_E_Tr,pol_σ_U_Tr,J_Tr,θ_Tr,Φ_Tr)=pol_val_results

    plot(pol_σ_U_Tr[1501:6000,1,1:5])
    plot!(pol_σ_U_Tr[1501:3000,1,end])

end


#####################
p1=plot(grid_a/(0.94*w[1,1]),pol_σ_U_Tr[1501:3000,1,1:3],title="2 years of tenure",label=["period 1" "period 2" "period 3"],linestyle=[:dash :dashdot :solid],lc=[:black :red :blue])
xlabel!("asset holdings")
ylabel!("probability of searching in occupation")

p2=plot(grid_a/(0.94*w[1,1]),pol_σ_U_Tr[3001:4500,1,1:3],title="5 years of tenure",legend=false,linestyle=[:dash :dashdot :solid],lc=[:black :red :blue])
xlabel!("asset holdings")
ylabel!("probability of searching in occupation")

p3=plot(grid_a/(0.94*w[1,1]),pol_σ_U_Tr[4501:6000,1,1:3],title="8 years of tenure",legend=false,linestyle=[:dash :dashdot :solid],lc=[:black :red :blue])
xlabel!("asset holdings")
ylabel!("probability of searching in occupation")

plot(p1,p2,p3,layout=(1,3))


## Computing wealth inequality

# Gini coefficient

cumulx=zeros(1500)
cumuly=zeros(1500)
for a_i in 1:1500
    cumulx[a_i]=Φ'*kron(ones(n_i*n_s+1+n_i*(n_s-1)),[ones(a_i);zeros(1500-a_i)])
    cumuly[a_i]=(Φ'*kron(ones(n_i*n_s+1+n_i*(n_s-1)),[grid_a[1:a_i];zeros(1500-a_i)]))/(Φ'*kron(ones(n_i*n_s+1+n_i*(n_s-1)),grid_a))
end

plot(cumulx,cumuly,legend=false)
plot!(cumulx,cumulx,legend=false)

Gini=(cumulx[1]-cumuly[1])*cumulx[1]
for a_i in 2:1500
    Gini+=(cumulx[a_i]-cumuly[a_i])*(cumulx[a_i]-cumulx[a_i-1])
end
Gini=Gini/0.5


# 90-10 ratio

q10=findmin(abs.(cumulx.-0.1))[2]
q90=findmin(abs.(cumulx.-0.9))[2]

ratio9010=grid_a[q90]/grid_a[q10]

# Computing income inequality


cumulx=zeros(n_s)
cumuly=zeros(n_s)
for s_i in 1:n_s
    cumulx[s_i]=sum(E[:,1:s_i])/sum(E)
    cumuly[s_i]=sum(E[:,1:s_i].*w[:,1:s_i])/sum(E.*w)
end

plot(cumulx,cumuly,legend=false)
plot!(cumulx,cumulx,legend=false)

Gini=(cumulx[1]-cumuly[1])*cumulx[1]
for s_i in 2:n_s
    Gini+=(cumulx[s_i]-cumuly[s_i])*(cumulx[s_i]-cumulx[s_i-1])
end
Gini=Gini/0.5
=#

# Skilled Labor Reallocation in an Incomplete Markets Economy

using Distributed, SharedArrays

addprocs(4-length(procs()))

path="C:\\Users\\joaor\\Dropbox\\Economics\\ThirdYearPaper\\Code\\ExperiencedLaborReallocationIncompleteMarkets\\New\\"
cd(path)

include(path*"ValueFunctionsIteration.jl")
include(path*"StationaryEquilibrium.jl")
include(path*"Transition.jl")
include(path*"EarningsLoss.jl")
include(path*"SeqSpaceJacobian.jl")
include(path*"bisection_derivative.jl")

@everywhere using Statistics,LinearAlgebra,Plots,SparseArrays,Interpolations,Optim,StatsBase,Distributions#,JLD

# Choices: i) Partial equilibrium or General Equilibrium, ii) Use multigrid?

PE=0 # If set to 0, code runs the GE, if set to 1 it runs the PE
small_grid=0
comp_transition=1
using_multigrid=1 # If set to 0, code runs just once with n_a grid points. If set to 1 it starts with n_a and increases the grid

# Calibration (1 period=2 months)

# Parameters
@everywhere β=0.9935 # Discount factor
@everywhere σ=2.0  # Inverse IES
@everywhere ρ=0.042 # Exogenous separation
@everywhere δ=1/(45*6) # Death
@everywhere α=[1/12;1/18;1/18;0] # Probability of being promoted
@everywhere χ1=0.6
@everywhere χ=[0;χ1;2/3*χ1;2/3*χ1] # Probability of losing skill
@everywhere φ=0.955 # Portion of good paid to worker
@everywhere b=0.45*φ # Replacement rate: Unemployment benefits; b*p > -̲a*r or c<0 at lowest wealth
@everywhere σ_ϵ=0.06# s.d. of taste shocks
@everywhere ξ=0.5 # Unemployed share in matching technology
@everywhere m=0.5 # Productivity of matching technology
@everywhere κ=0.1482 # Vacancy cost
@everywhere γ=[0.09592;0.14432;0.22198;0.53778] # Productivity share workers
@everywhere ω=-1e-7 # Elasticity of substitution between worker types (CES prod function)
@everywhere ν=3.0 # Elasticity of substitution between intermediate goods
@everywhere λ_2=4.23

@everywhere O=200
@everywhere ϕ=[1/(1+(O-1)^(1/ν));((O-1)^(1/ν))/(1+(O-1)^(1/ν))] #  weight of each occupation

@everywhere p_low=0.40505

@everywhere a_min=-9*p_low*φ
@everywhere a_max=20.0

# Prices
if PE==1
    #@everywhere p=[1 1.089 1.2 1.279; 1 1.089 1.2 1.279]*p_low
    @everywhere p=[0.4075 0.4440 0.4886 0.5205; 0.4075 0.4440 0.4886 0.5205]
end
@everywhere r=0.04/6

@everywhere z=2.0*ones(2) #n_o=2

# Grids
@everywhere n_beq=3
@everywhere n_o=2 # Should be equal to 2, because code is written for the case of 2 occupations, occupation 2 being large
@everywhere n_e=4

@everywhere grid_beq=[0.0;-190.0;-310.0]
@everywhere weight_beq=[0.28;0.56;0.16]

@everywhere grid_o=1:n_o
@everywhere grid_e=1:n_e


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
        pol_val_functions0=multigrid(nGrids_a,p)
        (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions0
    elseif using_multigrid==0
        grid_a=LinRange(a_min,a_max,n_a)
        grids=(grid_beq,grid_o,grid_e,grid_a)
        V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ=ValueFunctions(grids,p)
        pol_val_functions0=(V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ)
    end

    Check=(maximum(pol_a_E)<a_max)

    grid_a=LinRange(a_min,a_max,nGrids_a[end])
    grids=(grid_beq,grid_o,grid_e,grid_a)
    if n_anew!=nGrids_a[end]
        grid_a=LinRange(a_min,a_max,n_anew)
        pol_val_functions=transformVPolFunctions(pol_val_functions0,grids,grid_a)
        (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions
        grids=(grid_beq,grid_o,grid_e,grid_a)
    end

    pol_a_Ei,pol_a_Ui=transformPola(pol_a_E,pol_a_U,grids)

    pol_functions=(pol_a_Ei,pol_a_Ui,pol_σ_E,pol_σ_U,θ)
    Φ,Tr=ComputeDistribution(grids,pol_functions)

    Y,E,U=ComputeAggregates(grids,Φ,z)

    prices(ϕ,Y,z,E,e_i)=ϕ*γ[e_i]*Y^(1/ν)*z^(1-(1/ν))*E[e_i]^(ω-1)*sum(γ.*(E.^ω))^((1/ω)*(1-(1/ν))-1)

    prices_result=zeros(n_o,n_e)
    for o_i in 1:n_o
        for e_i in 1:n_e
            prices_result[o_i,e_i]=prices(ϕ[o_i],Y,z[o_i],E[o_i,:],e_i)
        end
    end
    display(prices_result)
    display(E[:,:])

    grid_a=LinRange(a_min,a_max,nGrids_a[end])
    grids=(grid_beq,grid_o,grid_e,grid_a)

    V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ=ValueFunctions(grids,prices_result;Guess=pol_val_functions0)
    pol_val_functions0=(V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ)

    grid_a=LinRange(a_min,a_max,n_anew)
    pol_val_functions=transformVPolFunctions(pol_val_functions0,grids,grid_a)
    (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions
    grids=(grid_beq,grid_o,grid_e,grid_a)

    pol_a_Ei,pol_a_Ui=transformPola(pol_a_E,pol_a_U,grids)

    StatEq=(V_E,V_U,W_E,W_U,pol_a_Ei,pol_a_Ui,pol_σ_E,pol_σ_U,J,θ,Φ,Y,E,U)

elseif PE==0
    StatEq,Tr=GeneralEquilibrium(z)

    grid_a=LinRange(a_min,a_max,n_anew)
    grids=(grid_beq,grid_o,grid_e,grid_a)
end

(V_E,V_U,W_E,W_U,pol_a_Ei,pol_a_Ui,pol_σ_E,pol_σ_U,J,θ,Φ,Y,E,U)=StatEq
pol_functions=(pol_a_Ei,pol_a_Ui,pol_σ_E,pol_σ_U,θ)

prices(ϕ,Y,z,E,e_i)=ϕ*γ[e_i]*Y^(1/ν)*z^(1-(1/ν))*E[e_i]^(ω-1)*sum(γ.*(E.^ω))^((1/ω)*(1-(1/ν))-1)

p=zeros(n_o,n_e)
for o_i in 1:n_o
    for e_i in 1:n_e
        p[o_i,e_i]=prices(ϕ[o_i],Y,z[o_i],E[o_i,:],e_i)
    end
end

Φnew=NewDistribution(grids,Φ)

StatEq=(V_E,V_U,W_E,W_U,pol_a_Ei,pol_a_Ui,pol_σ_E,pol_σ_U,J,θ,Φnew,Y,E,U)


earningsloss,wageloss=EarningsLoss(grids,pol_functions,Tr,p)

Jacobian=SeqSpaceJacobian(grids,StatEq,p,Tr)

shockdur=40
zt=z*ones(1,shockdur+1)
zt[:,1]=[1.7;2.0]
for t in 2:shockdur
    zt[:,t]=0.9*zt[:,t-1]+0.1*z
end

dZ=zeros(n_o,200)
dZ[:,1:shockdur+1]=zt.-z

dU,dUdZ=ComputingdU(grids,StatEq,Jacobian,dZ)

#=
# simulate
shock=rand(Normal(0,1),10000,2)
logzt=zeros(n_o,200)
dtotal=[zeros(n_o,n_e) for t in 1:10000]
for s in 1:10000-200+1
    for t in s:s-1+200
        if t==s
            logzt[:,t-s+1]=log.(z)+0.05*shock[s,:]
        else
            logzt[:,t-s+1]=0.9*logzt[:,t-s]+0.1*log.(z)
        end
    end
    zt=exp.(logzt)
    dZ=zt.-z
    dZ=dZ[:]
    dE=dUdZ*dZ
    for t in s:s-1+200
        for o_i in 1:n_o
            for e_i in 1:n_e
                dtotal[t][o_i,e_i]+=dE[(t-s)*n_o*n_e+(o_i-1)*n_e+e_i]
            end
        end
    end
end
Ett=zeros(n_o*n_e,10000)
unemp=zeros(10000)
for t in 1:10000
    for o_i in 1:n_o
        for e_i in 1:n_e
            Ett[(o_i-1)*n_e+e_i,t]=E[o_i,e_i]+dtotal[t][o_i,e_i]
        end
    end
    unemp[t]=max(1-sum(Ett[:,t]),0.0)
end
meanu=mean(unemp)
display(meanu)
=#

if comp_transition==1

    ###################################################
    #           Transition Exercises
    ###################################################


    T=200
    Eold=[zeros(n_o,n_e) for t in 1:T]
    for t in 1:T
        for o_i in 1:n_o
            for e_i in 1:n_e
                Eold[t][o_i,e_i]=E[o_i,e_i]+dU[(t-1)*n_o*n_e+(o_i-1)*n_e+e_i]
            end
        end
    end

    Uold=zeros(T)
    for t in 1:T
        Uold[t]=1-sum(Eold[t])
    end


    Eplot=zeros(n_e*n_o,T)
    for t in 1:T
        for o_i in 1:n_o
            for e_i in 1:n_e
                Eplot[(o_i-1)*n_e+e_i,t]=Eold[t][o_i,e_i]
            end
        end
    end

    plot0=plot(1:T,Eplot',legend=false)
    display(plot0)

    shockdur=40
    zt=z*ones(1,shockdur+1)
    zt[:,1]=[1.7;2.0]
    for t in 2:shockdur
        zt[:,t]=0.9*zt[:,t-1]+0.1*z
    end

    NewE,NewU,pol_val_results=Transition(grids,StatEq,zt;Guess=Eold)

    aggregates_transition=(NewE,NewU)

    #PlotResultsStatEq(grids,StatEq)
    #PlotResultsTransition(grids,zt,pol_val_results,aggregates_transition)

    (V_E_Tr,V_U_Tr,W_E_Tr,W_U_Tr,pol_a_Ei_Tr,pol_a_Ui_Tr,pol_σ_E_Tr,pol_σ_U_Tr,J_Tr,θ_Tr,Φ_Tr)=pol_val_results

    #plot(pol_σ_U_Tr[1501:6000,1,1:5])
    #plot!(pol_σ_U_Tr[1501:3000,1,end])

end


pt=[zeros(n_o,n_e) for t=1:T]
for t in 1:T
    Y=0.0
    y=zeros(n_o)
    for o_i in 1:n_o
        if t<shockdur
            z_o=zt[o_i,t]
        else
            z_o=zt[o_i,end]
        end
        y[o_i]=z_o*sum(γ.*(Eold[t][o_i,:].^ω))^(1/ω)
        Y+=ϕ[o_i]*y[o_i]^((ν-1)/ν)
    end
    Y=Y^(ν/(ν-1))

    for o_i in 1:n_o
        for e_i in 1:n_e
            if t<shockdur
                z_o=zt[o_i,t]
            else
                z_o=zt[o_i,end]
            end
            pt[t][o_i,e_i]=prices(ϕ[o_i],Y,z_o,Eold[t][o_i,:],e_i)
        end
    end
end

save(path*"results.jld", "StatEq", StatEq, "pol_val_results", pol_val_results,"NewE",NewE,"NewU",NewU)
#StatEq=load(path*"results.jld", "StatEq")
#pol_val_results=load(path*"results.jld", "pol_val_results")
#NewE=load(path*"results.jld", "NewE")
#NewU=load(path*"results.jld", "NewU")

#=
plot(grid_a[1:100]/(6*φ*p[1,1]),pol_σ_U_Tr[4501:4600,1,1],legend=false,lc=[:black],ylabel="probability of searching in occupation 1",xlabel="wealth-to-earnings ratio")
plot(P.(θ_Tr[4501,:]),title="2 years of tenure",legend=false,lc=[:black])


#####################
p1=plot(grid_a/(0.94*p[1,1]),pol_σ_U_Tr[1501:3000,1,1:3],title="2 years of tenure",label=["period 1" "period 2" "period 3"],linestyle=[:dash :dashdot :solid],lc=[:black :red :blue])
xlabel!("asset holdings")
ylabel!("probability of searching in occupation")

p2=plot(grid_a/(0.94*p[1,1]),pol_σ_U_Tr[3001:4500,1,1:3],title="5 years of tenure",legend=false,linestyle=[:dash :dashdot :solid],lc=[:black :red :blue])
xlabel!("asset holdings")
ylabel!("probability of searching in occupation")

p3=plot(grid_a/(0.94*p[1,1]),pol_σ_U_Tr[4501:6000,1,1:3],title="8 years of tenure",legend=false,linestyle=[:dash :dashdot :solid],lc=[:black :red :blue])
xlabel!("asset holdings")
ylabel!("probability of searching in occupation")

plot(p1,p2,p3,layout=(1,3))


=#

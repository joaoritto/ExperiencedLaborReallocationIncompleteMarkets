# Skilled Labor Reallocation in an Incomplete Markets Economy

using Statistics,LinearAlgebra,Plots,SparseArrays

# Parameters

β=0.96 # Discount factor
σ=2.0  # Inverse IES
h=0.3  # Value of leisure
ρ=0.04 # Exogenous separation
δ=0.01 # Separation with loss of skill
α=0.1 # Probability of becoming skilled
b=0.2
σ_ϵ=0.5
ξ=0.5
m=0.35
κ=0.9

# Prices
w=[0.5 1; 0.6 1.2]
r=0.03

# Job finding
η=0.8


# Grids
n_i=2
n_s=2
n_a=150
n_μ=150

grid_i=1:n_i
grid_s=1:n_s
grid_a=LinRange(-2,5,n_a)
grid_μ=LinRange(1e-2,1-1e-2,n_μ)

grids=(grid_i,grid_s,grid_a,grid_μ)

para=(β,σ,h,ρ,δ,α,w,η)


function VFunctionIter(para,grids,θ;Vguess=false)

    (β,σ,h,ρ,δ,α,w,η)=para
    (grid_i,grid_s,grid_a,grid_μ)=grids
    ϵ=1e-5

    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)

    nstates_E=n_μ*n_a*n_s*n_i
    ngrids_vars_E=[n_i,n_s,n_a,n_μ]
    nstates_U=n_a*n_s*n_i
    ngrids_vars_U=[n_i,n_s,n_a]
    nstates_S=n_a

    nsvars_E=4
    nsvars_U=3
    nsvars_S=1

    statestogrid_E=zeros(Int64,nstates_E,nsvars_E)
    for v in 1:nsvars_E
        statestogrid_E[:,v]=kron(ones(prod(ngrids_vars_E[1:v-1]),1),kron(1:ngrids_vars_E[v],ones(prod(ngrids_vars_E[v+1:nsvars_E]),1)))
    end

    statestogrid_U=zeros(Int64,nstates_U,nsvars_U)
    for v in 1:nsvars_U
        statestogrid_U[:,v]=kron(ones(prod(ngrids_vars_U[1:v-1]),1),kron(1:ngrids_vars_U[v],ones(prod(ngrids_vars_U[v+1:nsvars_U]),1)))
    end

    statestogrid_S=1:n_a

    if Vguess==false
        V_E_old=zeros(nstates_E)
        V_U_old=zeros(nstates_U)
        V_S_old=zeros(nstates_S)
        W_E_old=zeros(nstates_E)
        W_U_old=zeros(nstates_U)
        W_S_old=zeros(nstates_S)
    else
        V_E_old,V_U_old,V_S_old,W_E_old,W_U_old,W_S_old=Vguess
    end

    V_E_aux=zeros(n_a,nstates_E)
    V_U_aux=zeros(n_a,nstates_U)
    V_S_aux=zeros(n_a,nstates_S)
    V_E=zeros(nstates_E)
    V_U=zeros(nstates_U)
    V_S=zeros(nstates_S)
    W_E=zeros(nstates_E)
    W_U=zeros(nstates_U)
    W_S=zeros(nstates_S)

    pol_a_E=zeros(Int64,nstates_E)
    pol_a_U=zeros(Int64,nstates_U)
    pol_a_S=zeros(Int64,nstates_S)

    pol_μ_U=zeros(Int64,nstates_U)

    pol_σ_E=zeros(nstates_E)
    pol_σ_U=zeros(nstates_U)


    u(c)=if c>0 c^(1-σ)/(1-σ) else -Inf end
    p(θ)=min(m*θ^(1-ξ),1)
    error=1000

    iter=0
    while error>ϵ
        iter+=1

        # 1) V_E
        for ind in eachindex(V_E)
            μ_i=statestogrid_E[ind,4]
            a_i=statestogrid_E[ind,3]
            s_i=statestogrid_E[ind,2]
            i_i=statestogrid_E[ind,1]

            wage=w[i_i,s_i]

            for a1_i in 1:n_a
                if (1+r)*grid_a[a_i]+grid_μ[μ_i]*wage-grid_a[a1_i]<0
                    V_E_aux[a1_i,ind]=-Inf
                else
                    if s_i==1
                        ind1_u=a1_i
                        ind1_e1=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        ind1_e2=[i_i-1,2-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        V_E_aux[a1_i,ind]=u((1+r)*grid_a[a_i]+grid_μ[μ_i]*wage-grid_a[a1_i])+β*((ρ-δ)*W_U_old[ind1_u]+δ*W_S_old[a1_i]+(1-ρ)*((1-α)*W_E_old[ind1_e1]+α*W_E_old[ind1_e2]))
                    elseif s_i==2
                        ind1_u=[i_i-1,s_i-1,a_i]'*[n_s*n_a,n_a,1]
                        ind1_e=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        V_E_aux[a1_i,ind]=u((1+r)*grid_a[a_i]+grid_μ[μ_i]*wage-grid_a[a1_i])+β*((ρ-δ)*W_U_old[ind1_u]+δ*W_S_old[a1_i]+(1-ρ)*W_E_old[ind1_e])
                    end
                end
            end
            V_E[ind],pol_a_E[ind]=findmax(V_E_aux[:,ind])
        end

        # 2) V_U
        for ind in eachindex(V_U)
            a_i=statestogrid_U[ind,3]
            s_i=statestogrid_U[ind,2]
            i_i=statestogrid_U[ind,1]

            for a1_i in 1:n_a
                if (1+r)*grid_a[a_i]+b-grid_a[a1_i]<0
                    V_U_aux[a1_i,ind]=-Inf
                else
                    ind1=[i_i-1,s_i-1,a1_i]'*[n_s*n_a,n_a,1]
                    V_U_aux[a1_i,ind]=u((1+r)*grid_a[a_i]+b-grid_a[a1_i])+β*(W_U_old[ind1])
                end
            end
            V_U[ind],pol_a_U[ind]=findmax(V_U_aux[:,ind])
        end

        # 3) V_S
        for ind in eachindex(V_S)
            a_i=statestogrid_S[ind,1]

            for a1_i in 1:n_a
                if (1+r)*grid_a[a_i]+b-grid_a[a1_i]<0
                    V_S_aux[a1_i,ind]=-Inf
                else
                    ind1=a1_i
                    V_S_aux[a1_i,ind]=u((1+r)*grid_a[a_i]+b-grid_a[a1_i])+β*(W_S_old[ind1])
                end
            end
            V_S[ind],pol_a_S[ind]=findmax(V_S_aux[:,ind])
        end

        # Now compute Ws
        # 1) W_E
        for ind in eachindex(W_E)
            μ_i=statestogrid_E[ind,4]
            a_i=statestogrid_E[ind,3]
            s_i=statestogrid_E[ind,2]
            i_i=statestogrid_E[ind,1]

            W_E[ind]=σ_ϵ*log(exp(V_E[ind]/σ_ϵ)+exp(W_S_old[a_i]/σ_ϵ))
            pol_σ_E[ind]=exp(V_E[ind]/σ_ϵ)/(exp(V_E[ind]/σ_ϵ)+exp(V_S[a_i]/σ_ϵ))
        end

        subiter=0
        error1=1000
        while error1>ϵ
            subiter+=1
            # 2) W_U
            for ind in eachindex(W_U)
                a_i=statestogrid_U[ind,3]
                s_i=statestogrid_U[ind,2]
                i_i=statestogrid_U[ind,1]

                Val_temp=zeros(n_μ)
                for μ_i in 1:n_μ
                    ste=[i_i-1,s_i-1,a_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    prob=p(θ[ste])
                    Val_temp[μ_i]=prob*V_E[ste]+(1-prob)*V_U[ind]
                end
                V_high,pol_μ_U[ind]=findmax(Val_temp)

                W_U[ind]=σ_ϵ*log(exp(V_high/σ_ϵ)+exp(W_S_old[a_i]/σ_ϵ))
                pol_σ_U[ind]=exp(V_high/σ_ϵ)/(exp(V_high/σ_ϵ)+exp(V_S[a_i]/σ_ϵ))
            end

            # 3) W_S
            for a_i in 1:n_a
                ind=a_i

                ind1=zeros(Int64,n_i)
                for i_i in 1:n_i
                    ind1[i_i]=[i_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]
                end
                W_S[ind]=mean(W_U[ind1])
            end
            error1=maximum((W_S-W_S_old).^2)
            W_S,W_S_old=W_S_old,W_S
            #println("subiter ",subiter," error ",error1)
        end

        error=maximum((W_E-W_E_old).^2)+maximum((W_U-W_U_old).^2)
        println("iter ",iter,": error=",error)
        W_E,W_E_old=W_E_old,W_E
        W_U,W_U_old=W_U_old,W_U
        #W_S,W_S_old=W_S_old,W_S
    end

    W_E,W_E_old=W_E_old,W_E
    W_U,W_U_old=W_U_old,W_U
    W_S,W_S_old=W_S_old,W_S

    return V_E,V_U,V_S,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U
end



function JFunctionIter(para,grids,policyfunctions_W; Jguess=false)

    (β,σ,h,ρ,δ,α,w,η)=para
    (grid_i,grid_s,grid_a,grid_μ)=grids
    (pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U)=policyfunctions_W
    ϵ=1e-6

    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)


    nstates=n_μ*n_a*n_s*n_i
    nsvars=4
    ngrids_vars=[n_i,n_s,n_a,n_μ]


    statestogrid=zeros(Int64,nstates,nsvars)
    for v in 1:nsvars
        statestogrid[:,v]=kron(ones(prod(ngrids_vars[1:v-1]),1),kron(1:ngrids_vars[v],ones(prod(ngrids_vars[v+1:nsvars]),1)))
    end

    if Jguess==false
        J_old=zeros(nstates)
    else
        J_old=zeros(nstates)
        J_old[:]=Jguess
    end

    J=zeros(nstates)

    error=1000
    iter=0
    while error>ϵ
        iter+=1
        for ind in eachindex(J)
            μ_i=statestogrid[ind,4]
            a_i=statestogrid[ind,3]
            s_i=statestogrid[ind,2]
            i_i=statestogrid[ind,1]

            wage=w[i_i,s_i]
            a1_i=pol_a_E[ind]

            if s_i==1
                ind1_0=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_1=[i_i-1,2-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                σ_prob0=pol_σ_E[ind1_0]
                σ_prob1=pol_σ_E[ind1_1]
                J[ind]=(1-grid_μ[μ_i])*wage+β*(1-ρ)*((1-α)*σ_prob0*J_old[ind1_0]+α*σ_prob1*J_old[ind1_1])
            elseif s_i==2
                ind1=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                σ_prob=pol_σ_E[ind1]
                J[ind]=(1-grid_μ[μ_i])*wage+β*(1-ρ)*σ_prob*J_old[ind1]
            end
        end

        error=maximum((J-J_old).^2)
        println("iter ",iter,": error=",error)
        J,J_old=J_old,J
    end
    J,J_old=J_old,J

    return J
end


function ValueFunctions(para,grids)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)

    q_inv(y)=(y/m)^(-1/ξ)

    θ=zeros(n_μ*n_a*n_s*n_i)
    J_old=5.0*ones(n_μ*n_a*n_s*n_i)
    Vfunctions=false
    policyfunctions=false

    ϵ=1e-6
    error=1000

    iter=0
    while error>ϵ
        iter+=1
        for ind in eachindex(θ)
            θ[ind]=q_inv(κ/J_old[ind])
        end
        V_E,V_U,V_S,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U=VFunctionIter(para,grids,θ,Vguess=Vfunctions)

        policyfunctions=(pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U)
        Vfunctions=(V_E,V_U,V_S,W_E,W_U,W_S)

        J=JFunctionIter(para,grids,policyfunctions,Jguess=J_old)

        error=maximum((J-J_old).^2)
        println("iter ",iter," in outward loop, error of ",error)
        J,J_old=J_old,J
    end
    J,J_old=J_old,J

    for ind in eachindex(θ)
        θ[ind]=q_inv(κ/J[ind])
    end
    (V_E,V_U,V_S,W_E,W_U,W_S)=Vfunctions
    (pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U)=policyfunctions

    return V_E,V_U,V_S,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U,J,θ
end

function ComputeDistribution(para,grids,pol_val_functions)

    p(θ)=min(m*θ^(1-ξ),1)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    (V_E,V_U,V_S,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions

    nstates=n_i*n_s*n_a*n_μ+n_i*n_s*n_a+n_a
    nsvars=5
    ngrids_vars=[3,n_i,n_s,n_a,n_μ]

    nstates_E=n_i*n_s*n_a*n_μ
    statestogrid_E=ones(Int64,nstates_E,nsvars)
    nstates_U=n_i*n_s*n_a
    statestogrid_U=ones(Int64,nstates_U,nsvars)
    nstates_S=n_a
    statestogrid_S=ones(Int64,nstates_S,nsvars)
    for v in 1:nsvars
        if v==1
            statestogrid_E[:,v]=1*ones(nstates_E,1)
            statestogrid_U[:,v]=2*ones(nstates_U,1)
            statestogrid_S[:,v]=3*ones(nstates_S,1)
        else
            statestogrid_E[:,v]=kron(ones(prod(ngrids_vars[2:v-1]),1),kron(1:ngrids_vars[v],ones(prod(ngrids_vars[v+1:nsvars]),1)))
            if v==nsvars
            else
                statestogrid_U[:,v]=kron(ones(prod(ngrids_vars[2:v-1]),1),kron(1:ngrids_vars[v],ones(prod(ngrids_vars[v+1:nsvars-1]),1)))
            end
            if v==4
                statestogrid_S[:,v]=1:n_a
            end
        end
    end
    statestogrid=[statestogrid_E;statestogrid_U;statestogrid_S]

    # MUST USE SPARSE MATRICES
    # Construct Transition matrix
    T1=spzeros(nstates,nstates)
    T2=spzeros(nstates,nstates)

    for ind in 1:nstates
        μ_i=statestogrid[ind,5]
        a_i=statestogrid[ind,4]
        s_i=statestogrid[ind,3]
        i_i=statestogrid[ind,2]
        e_i=statestogrid[ind,1]

        if e_i==1

            indE=[i_i-1,s_i-1,a_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
            a1_i=pol_a_E[indE]
            if s_i==1
                ind1_e0=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_e1=[i_i-1,2-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_u=n_i*n_s*n_a*n_μ+[i_i-1,s_i-1,a1_i]'*[n_s*n_a,n_a,1]


                T1[ind,ind1_e0]=(1-ρ)*(1-α)
                T1[ind,ind1_e1]=(1-ρ)*α
                T1[ind,ind1_u]=ρ-δ

                ind1_s=zeros(Int64,n_i)
                for i1_i in 1:n_i
                    ind1_s[i1_i]=n_i*n_s*n_a*n_μ+[i1_i-1,1-1,a1_i]'*[n_s*n_a,n_a,1]
                end

                for i1_i in 1:n_i
                    T1[ind,ind1_s[i1_i]]+=δ*(1/n_i)
                end

            elseif s_i==2
                ind1_e=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_u=n_i*n_s*n_a*n_μ+[i_i-1,s_i-1,a1_i]'*[n_s*n_a,n_a,1]
                ind1_s=n_i*n_s*n_a*n_μ+n_i*n_s*n_a+a1_i

                T1[ind,ind1_e]=(1-ρ)
                T1[ind,ind1_u]=ρ-δ

                ind1_s=zeros(Int64,n_i)
                for i1_i in 1:n_i
                    ind1_s[i1_i]=n_i*n_s*n_a*n_μ+[i1_i-1,1-1,a1_i]'*[n_s*n_a,n_a,1]
                end

                for i1_i in 1:n_i
                    T1[ind,ind1_s[i1_i]]+=δ*(1/n_i)
                end
            end


            ind1_s_u=zeros(Int64,n_i)
            ind1_s_e=zeros(Int64,n_i)
            for i1_i in 1:n_i
                ind1_s_u[i1_i]=n_i*n_s*n_a*n_μ+[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]
                ind1_s_e[i1_i]=[i1_i-1,1-1,a_i-1,pol_μ_U[[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
            end

            T2[ind,ind]=pol_σ_E[indE]
            for i1_i in 1:n_i
                T2[ind,ind1_s_u[i1_i]]+=(1-pol_σ_E[indE])*(1/n_i)*(1-p(θ[[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]]))
                T2[ind,ind1_s_e[i1_i]]+=(1-pol_σ_E[indE])*(1/n_i)*p(θ[[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]])
            end

        elseif e_i==2

            indU=[i_i-1,s_i-1,a_i]'*[n_s*n_a,n_a,1]
            a1_i=pol_a_U[indU]

            ind1_u=n_i*n_s*n_a*n_μ+[i_i-1,s_i-1,a1_i]'*[n_s*n_a,n_a,1]

            T1[ind,ind1_u]=1

            ind1_s_u=zeros(Int64,n_i)
            ind1_s_e=zeros(Int64,n_i)
            for i1_i in 1:n_i
                ind1_s_u[i1_i]=n_i*n_s*n_a*n_μ+[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]
                ind1_s_e[i1_i]=[i1_i-1,1-1,a_i-1,pol_μ_U[[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
            end
            ind1_e=[i_i-1,s_i-1,a_i-1,pol_μ_U[indU]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]


            T2[ind,ind1_u]=pol_σ_U[indU]*(1-p(θ[indU]))
            T2[ind,ind1_e]=pol_σ_U[indU]*p(θ[indU])
            for i1_i in 1:n_i
                T2[ind,ind1_s_u[i1_i]]+=(1-pol_σ_U[indU])*(1/n_i)*(1-p(θ[[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]]))
                T2[ind,ind1_s_e[i1_i]]+=(1-pol_σ_U[indU])*(1/n_i)*p(θ[[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]])
            end

        end
    end
    T=spzeros(nstates,nstates)
    T=T1'*T2'

    Φ=zeros(nstates)
    ind0=200
    Φ[ind0]=1.0

    for i in 1:300
        Φ=T*Φ
    end

    return Φ
end

function ComputeAggregates(para,grids,pol_val_functions,Φ)
    p(θ)=min(m*θ^(1-ξ),1)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    (V_E,V_U,V_S,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions

    nstates=n_i*n_s*n_a*n_μ+n_i*n_s*n_a+n_a
    nsvars=5
    ngrids_vars=[3,n_i,n_s,n_a,n_μ]

    # Compute unemployed and employed of each type
    E_U=zeros(n_i)
    E_E=zeros(n_i)
    U_U=zeros(n_i)
    U_E=zeros(n_i)
    for i_i in 1:n_i
        E_U[i_i]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+1:(i_i-1)*n_s*n_a*n_μ+n_μ*n_a])
        E_E[i_i]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+n_μ*n_a+1:(i_i-1)*n_s*n_a*n_μ+2*n_μ*n_a])
        U_U[i_i]=sum(Φ[n_i*n_s*n_a*n_μ+(i_i-1)*n_s*n_a+1:n_i*n_s*n_a*n_μ+(i_i-1)*n_s*n_a+n_a])
        U_E[i_i]=sum(Φ[n_i*n_s*n_a*n_μ+(i_i-1)*n_s*n_a+n_a+1:n_i*n_s*n_a*n_μ+(i_i-1)*n_s*n_a+2*n_a])
    end

end

function PlotResults(para,grids,pol_val_functions,Φ)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    (V_E,V_U,V_S,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions

    nstates=n_i*n_s*n_a*n_μ+n_i*n_s*n_a+n_a
    nsvars=5
    ngrids_vars=[3,n_i,n_s,n_a,n_μ]

    statestogrid=zeros(nstates)

    nstates_E=n_i*n_s*n_a*n_μ
    statestogrid_E=ones(Int64,nstates_E,nsvars)
    nstates_U=n_i*n_s*n_a
    statestogrid_U=ones(Int64,nstates_U,nsvars)
    nstates_S=n_a
    statestogrid_S=ones(Int64,nstates_S,nsvars)
    for v in 1:nsvars
        if v==1
            statestogrid_E[:,v]=1*ones(nstates_E,1)
            statestogrid_U[:,v]=2*ones(nstates_U,1)
            statestogrid_S[:,v]=3*ones(nstates_S,1)
        else
            statestogrid_E[:,v]=kron(ones(prod(ngrids_vars[2:v-1]),1),kron(1:ngrids_vars[v],ones(prod(ngrids_vars[v+1:nsvars]),1)))
            if v==nsvars
            else
                statestogrid_U[:,v]=kron(ones(prod(ngrids_vars[2:v-1]),1),kron(1:ngrids_vars[v],ones(prod(ngrids_vars[v+1:nsvars-1]),1)))
            end
            if v==4
                statestogrid_S[:,v]=1:n_a
            end
        end
    end
    statestogrid=[statestogrid_E;statestogrid_U]

    plot(pol_μ_U[1:n_a])
    plot!(pol_μ_U[n_a+1:2*n_a])
    plot!(pol_μ_U[2*n_a+1:3*n_a])
    plot!(pol_μ_U[3*n_a+1:4*n_a])

    μ_a=zeros(n_a)
    for a_i in 1:n_a
        μ_a[a_i]=Φ[(a_i-1)*n_μ+1:a_i*n_μ]'*grid_μ[statestogrid[(a_i-1)*n_μ+1:a_i*n_μ,5]]/sum(Φ[(a_i-1)*n_μ+1:a_i*n_μ])
    end

    plot(μ_a)

    μ_a2=zeros(n_a)
    for a_i in 1:n_a
        μ_a2[a_i]=Φ[n_a*n_μ+(a_i-1)*n_μ+1:n_a*n_μ+a_i*n_μ]'*grid_μ[statestogrid[n_a*n_μ+(a_i-1)*n_μ+1:n_a*n_μ+a_i*n_μ,5]]/sum(Φ[n_a*n_μ+(a_i-1)*n_μ+1:n_a*n_μ+a_i*n_μ])
    end

    plot!(μ_a2)
    yaxis!([0,1])

    d_a=zeros(n_a)
    for a_i in 1:n_a
        d_a[a_i]=sum(Φ[(a_i-1)*n_μ+1:a_i*n_μ])
    end
    plot(d_a)


end

V_E,V_U,V_S,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U,J,θ=ValueFunctions(para,grids)

pol_val_functions=(V_E,V_U,V_S,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)

Φ=ComputeDistribution(para,grids,pol_val_functions)

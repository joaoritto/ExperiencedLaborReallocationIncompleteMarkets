
# Workers' value function iteration

function VFunctionIter(grids,θ;Vguess=false,tol=false)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    if tol==false
        ϵ=1e-5
    else
        ϵ=tol
    end

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
        W_E_old=zeros(nstates_E)
        W_U_old=zeros(nstates_U)
        W_S_old=zeros(nstates_S)
    else
        V_E_old,V_U_old,W_E_old,W_U_old,W_S_old=Vguess
    end

    V_E_aux=zeros(n_a,nstates_E)
    V_U_aux=zeros(n_a,nstates_U)
    V_S_aux=zeros(n_a,nstates_S)
    V_E=zeros(nstates_E)
    V_U=zeros(nstates_U)
    W_E=zeros(nstates_E)
    W_U=zeros(nstates_U)
    W_S=zeros(nstates_S)

    pol_a_E=zeros(Int64,nstates_E)
    pol_a_U=zeros(Int64,nstates_U)

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
                if a1_i>1 && V_E_aux[a1_i,ind]<V_E_aux[a1_i-1,ind]
                    V_E[ind]=V_E_aux[a1_i-1,ind]
                    pol_a_E[ind]=a1_i-1
                    break
                elseif a1_i==n_a
                    V_E[ind]=V_E_aux[a1_i,ind]
                    pol_a_E[ind]=a1_i
                end
            end
            #V_E[ind],pol_a_E[ind]=V_E_aux[a1_i-1,ind]
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
                if a1_i>1 && V_U_aux[a1_i,ind]<V_U_aux[a1_i-1,ind]
                    V_U[ind]=V_U_aux[a1_i-1,ind]
                    pol_a_U[ind]=a1_i-1
                    break
                elseif a1_i==n_a
                    V_U[ind]=V_U_aux[a1_i,ind]
                    pol_a_U[ind]=a1_i
                end
            end
            #V_U[ind],pol_a_U[ind]=findmax(V_U_aux[:,ind])
        end


        # Now compute Ws
        # 1) W_E
        for ind in eachindex(W_E)
            μ_i=statestogrid_E[ind,4]
            a_i=statestogrid_E[ind,3]
            s_i=statestogrid_E[ind,2]
            i_i=statestogrid_E[ind,1]

            W_E[ind]=σ_ϵ*log(exp(V_E[ind]/σ_ϵ)+exp(W_S_old[a_i]/σ_ϵ))
            pol_σ_E[ind]=exp(V_E[ind]/σ_ϵ)/(exp(V_E[ind]/σ_ϵ)+exp(W_S_old[a_i]/σ_ϵ))
        end

        subiter=0
        error1=1000
        while error1>1e-4
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
                pol_σ_U[ind]=exp(V_high/σ_ϵ)/(exp(V_high/σ_ϵ)+exp(W_S_old[a_i]/σ_ϵ))
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
        end

        error=maximum((W_E-W_E_old).^2)+maximum((W_U-W_U_old).^2)
        println("iter ",iter,": error=",error)
        W_E,W_E_old=W_E_old,W_E
        W_U,W_U_old=W_U_old,W_U
    end

    W_E,W_E_old=W_E_old,W_E
    W_U,W_U_old=W_U_old,W_U
    W_S,W_S_old=W_S_old,W_S

    return V_E,V_U,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U
end

# Firms' value function iteration

function JFunctionIter(grids,policyfunctions_W; Jguess=false,tol=false)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    (pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U)=policyfunctions_W
    if tol==false
        ϵ=1e-6
    else
        ϵ=tol
    end

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

# Iteration of firms and workers' value functions

function ValueFunctions(grids;Guess=false)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)

    q_inv(y)=(y/m)^(-1/ξ)

    θ=zeros(n_μ*n_a*n_s*n_i)
    if Guess==false
        J_old=5.0*ones(n_μ*n_a*n_s*n_i)
        Vfunctions=false
        policyfunctions=false
    else
        (V_E_old,V_U_old,V_S_old,W_E_old,W_U_old,W_S_old,pol_a_E_old,pol_a_U_old,pol_a_S_old,pol_μ_U_old,pol_σ_E_old,pol_σ_U_old,J_old,θ_old)=Guess
        Vfunctions=(V_E_old,V_U_old,V_S_old,W_E_old,W_U_old,W_S_old)
        policyfunctions=(pol_a_E_old,pol_a_U_old,pol_a_S_old,pol_μ_U_old,pol_σ_E_old,pol_σ_U_old)
    end


    ϵ=1e-6
    error=1000

    iter=0
    while error>ϵ && iter<100
        iter+=1
        for ind in eachindex(θ)
            θ[ind]=q_inv(κ/J_old[ind])
        end

        V_E,V_U,V_S,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U=VFunctionIter(grids,θ,Vguess=Vfunctions)

        policyfunctions=(pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U)
        Vfunctions=(V_E,V_U,V_S,W_E,W_U,W_S)

        J=JFunctionIter(grids,policyfunctions,Jguess=J_old)

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

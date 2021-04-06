
# Workers' value function iteration

function VFunctionIterEq(grids,w,θ;Vguess=false,tol=false)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    if tol==false
        ϵ=1e-5
    else
        ϵ=tol
    end

    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)

    nstates_E=n_μ*n_a*n_s*n_i
    ngrids_vars_E=[n_i,n_s,n_a,n_μ]
    nstates_U=n_a+n_i*n_a
    ngrids_vars_U=[n_i,n_s,n_a]

    nsvars_E=4
    nsvars_U=3

    statestogrid_E=zeros(Int64,nstates_E,nsvars_E)
    for v in 1:nsvars_E
        statestogrid_E[:,v]=kron(ones(prod(ngrids_vars_E[1:v-1]),1),kron(1:ngrids_vars_E[v],ones(prod(ngrids_vars_E[v+1:nsvars_E]),1)))
    end

    statestogrid_U=zeros(Int64,nstates_U,nsvars_U)
    statestogrid_U[1:n_a,:]=hcat(ones(n_a,2),1:n_a)
    for i_i in 1:n_i
        statestogrid_U[n_a+(i_i-1)*n_a+1:n_a+i_i*n_a,:]=hcat(i_i*ones(n_a,1),2*ones(n_a,1),1:n_a)
    end


    if Vguess==false
        V_E_old=zeros(nstates_E)
        V_U_old=zeros(nstates_U)
        W_E_old=zeros(nstates_E)
        W_U_old=zeros(nstates_U)
    else
        V_E_old,V_U_old,W_E_old,W_U_old=Vguess
    end

    V_E_aux=zeros(n_a,nstates_E)
    V_U_aux=zeros(n_a,nstates_U)
    V_E=zeros(nstates_E)
    V_U=zeros(nstates_U)
    W_E=zeros(nstates_E)
    W_U=zeros(nstates_U)

    pol_a_E=zeros(Int64,nstates_E)
    pol_a_U=zeros(Int64,nstates_U)

    pol_μ_U=zeros(Int64,nstates_U,n_i)

    pol_σ_E=zeros(nstates_E)
    pol_σ_U=zeros(nstates_U,n_i)


    u(c)=if c>0 c^(1-σ)/(1-σ) else -Inf end
    p(θ)=min(m*θ^(1-ξ),1)
    error=1000

    iter=0
    while error>ϵ && iter<1000
        iter+=1

        # 1) V_E
        for ind in eachindex(V_E)
            μ_i=statestogrid_E[ind,4]
            a_i=statestogrid_E[ind,3]
            s_i=statestogrid_E[ind,2]
            i_i=statestogrid_E[ind,1]

            wage=w[i_i,s_i]

            #if a_i==1
                a1_start=1
            #else
            #    a1_start=max(pol_a_E[[i_i-1,s_i-1,(a_i-1)-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]],1)
            #end

            for a1_i in a1_start:n_a
                c=(1+r)*grid_a[a_i]+grid_μ[μ_i]*wage-grid_a[a1_i]
                if c<0
                    V_E_aux[a1_i,ind]=-Inf
                else
                    if s_i==1
                        ind1_u=a1_i
                        ind1_e1=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        ind1_e2=[i_i-1,2-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        V_E_aux[a1_i,ind]=u(c)+β*(ρ*V_U_old[ind1_u]+(1-ρ)*((1-α)*W_E_old[ind1_e1]+α*W_E_old[ind1_e2]))
                    elseif s_i==2
                        ind1_u=n_a+[i_i-1,a1_i]'*[n_a,1]
                        ind1_e=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        V_E_aux[a1_i,ind]=u(c)+β*((ρ-δ)*V_U_old[ind1_u]+δ*V_U_old[a1_i]+(1-ρ)*W_E_old[ind1_e])
                    end
                end
            #    if a1_i>a1_start && V_E_aux[a1_i,ind]<V_E_aux[a1_i-1,ind]
            #        V_E[ind]=V_E_aux[a1_i-1,ind]
            #        pol_a_E[ind]=a1_i-1
            #        break
            #    elseif a1_i==n_a
            #        V_E[ind]=V_E_aux[a1_i,ind]
            #        pol_a_E[ind]=a1_i
            #    end
            end
            V_E[ind],pol_a_E[ind]=findmax(V_E_aux[a1_start:end,ind])
        end

        # 2) V_U
        for ind in eachindex(V_U)
            a_i=statestogrid_U[ind,3]
            s_i=statestogrid_U[ind,2]
            i_i=statestogrid_U[ind,1]

            #if a_i==1
                a1_start=1
            #else
            #    a1_start=max(pol_a_U[ind-1],1)
            #end

            for a1_i in a1_start:n_a
                c=(1+r)*grid_a[a_i]+b-grid_a[a1_i]
                if c<0
                    V_U_aux[a1_i,ind]=-Inf
                else
                    if s_i==1
                        ind1=a1_i
                    elseif s_i==2
                        ind1=n_a+[i_i-1,a1_i]'*[n_a,1]
                    end
                    V_U_aux[a1_i,ind]=u(c)+β*W_U_old[ind1]
                end
            #    if a1_i>1 && V_U_aux[a1_i,ind]<V_U_aux[a1_i-1,ind]
            #        V_U[ind]=V_U_aux[a1_i-1,ind]
            #        pol_a_U[ind]=a1_i-1
            #        break
            #    elseif a1_i==n_a
            #        V_U[ind]=V_U_aux[a1_i,ind]
            #        pol_a_U[ind]=a1_i
            #    end
            end
            V_U[ind],pol_a_U[ind]=findmax(V_U_aux[a1_start:end,ind])
        end


        # Now compute Ws

        # 1) W_E
        for ind in eachindex(W_E)
            μ_i=statestogrid_E[ind,4]
            a_i=statestogrid_E[ind,3]
            s_i=statestogrid_E[ind,2]
            i_i=statestogrid_E[ind,1]

            if s_i==1
                W_E[ind]=σ_ϵ*((V_E[ind]/σ_ϵ)+log(1+exp((V_U[a_i]-V_E[ind])/σ_ϵ)))
                pol_σ_E[ind]=(1+exp((V_U[a_i]-V_E[ind])/σ_ϵ))^(-1)
            elseif s_i==2
                W_E[ind]=σ_ϵ*((V_E[ind]/σ_ϵ)+log(1+exp((V_U[i_i*n_a+a_i]-V_E[ind])/σ_ϵ)))
                pol_σ_E[ind]=(1+exp((V_U[i_i*n_a+a_i]-V_E[ind])/σ_ϵ))^(-1)
            end
        end

        # 2) W_U
        V_job_s=zeros(length(W_U),n_i)
        for ind in eachindex(W_U)
            a_i=statestogrid_U[ind,3]
            s_i=statestogrid_U[ind,2]
            i_i=statestogrid_U[ind,1]

            Val_temp=zeros(n_μ)


            if s_i==1
                for i1_i in 1:n_i
                    for μ1_i in 1:n_μ
                        ste=[i1_i-1,1-1,a_i-1,μ1_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        prob=p(θ[ste])
                        Val_temp[μ1_i]=prob*V_E[ste]+(1-prob)*V_U[ind]
                    end
                    V_job_s[ind,i1_i],pol_μ_U[ind,i1_i]=findmax(Val_temp)
                end

                W_U[ind]=σ_ϵ*((V_job_s[ind,1]/σ_ϵ)+log(1+sum(exp.((V_job_s[ind,2:end].-V_job_s[ind,1])/σ_ϵ))))
                for i1_i in 1:n_i
                    pol_σ_U[ind,i1_i]=(sum(exp.((V_job_s[ind,:].-V_job_s[ind,i1_i])/σ_ϵ)))^(-1)
                end
            elseif s_i==2
                for i1_i in 1:n_i
                    for μ1_i in 1:n_μ
                        if i1_i==i_i
                            s1_i=2
                        else
                            s1_i=1
                        end
                        ste=[i1_i-1,s1_i-1,a_i-1,μ1_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        prob=p(θ[ste])
                        Val_temp[μ1_i]=prob*V_E[ste]+(1-prob)*V_U[ind]
                    end
                    V_job_s[ind,i1_i],pol_μ_U[ind,i1_i]=findmax(Val_temp)
                end

                W_U[ind]=σ_ϵ*((V_job_s[ind,1]/σ_ϵ)+log(1+sum(exp.((V_job_s[ind,2:end].-V_job_s[ind,1])/σ_ϵ))))
                for i1_i in 1:n_i
                    pol_σ_U[ind,i1_i]=(sum(exp.((V_job_s[ind,:].-V_job_s[ind,i1_i])/σ_ϵ)))^(-1)
                end
            end
        end

        error=maximum((W_E-W_E_old).^2)+maximum((W_U-W_U_old).^2)
        println("iter ",iter,": error=",error)
        W_E,W_E_old=W_E_old,W_E
        W_U,W_U_old=W_U_old,W_U
    end

    W_E,W_E_old=W_E_old,W_E
    W_U,W_U_old=W_U_old,W_U

    return V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U
end

# Firms' value function iteration

function JFunctionIter(grids,w,policyfunctions_W; Jguess=false,tol=false)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    (pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U)=policyfunctions_W
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

function ValueFunctions(grids,w;Guess=false)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)

    q_inv(y)=(y/m)^(-1/ξ)

    θ=zeros(n_μ*n_a*n_s*n_i)
    if Guess==false
        J_old=0.5*ones(n_μ*n_a*n_s*n_i)
        Vfunctions=false
        policyfunctions=false
    else
        (V_E_old,V_U_old,W_E_old,W_U_old,pol_a_E_old,pol_a_U_old,pol_μ_U_old,pol_σ_E_old,pol_σ_U_old,J_old,θ_old)=Guess
        Vfunctions=(V_E_old,V_U_old,W_E_old,W_U_old)
        policyfunctions=(pol_a_E_old,pol_a_U_old,pol_μ_U_old,pol_σ_E_old,pol_σ_U_old)
    end


    ϵ=1e-6
    error=1000

    iter=0
    while error>1e-6 && iter<50
        iter+=1
        for ind in eachindex(θ)
            θ[ind]=q_inv(κ/J_old[ind])
        end

        V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U=VFunctionIterEq(grids,w,θ,Vguess=Vfunctions,tol=1e-6)

        policyfunctions=(pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U)
        Vfunctions=(V_E,V_U,W_E,W_U)

        J=JFunctionIter(grids,w,policyfunctions,Jguess=J_old,tol=1e-9)

        error=maximum((J-J_old).^2)
        println("iter ",iter," in outward loop, error of ",error)
        J,J_old=J_old,J
    end
    J,J_old=J_old,J

    for ind in eachindex(θ)
        θ[ind]=q_inv(κ/J[ind])
    end
    (V_E,V_U,W_E,W_U)=Vfunctions
    (pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U)=policyfunctions

    return V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ
end

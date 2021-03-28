# Computing a transition

function Transition(grids,StatEq,zt;Guess=false)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)

    shockdur=size(zt,2)-1
    ϵ=1e-6
    wupdate=0.5
    # Is shock permanent? (or transitory?)
    aux=zeros(n_i)
    for i_i in 1:n_i
        aux[i_i]=abs(zt[i_i,t]/zt[i_i,t]-1)
    end
    if maximum(aux)==0
        permanent=0
    else
        permanent=1
    end

    if permanent==1
        (StatEq1,StatEq2)=StatEq
        (V_E0,V_U0,W_E0,W_U0,pol_a_E0,pol_a_U0,pol_μ_U0,pol_σ_E0,pol_σ_U0,J0,θ0,Φ0,Y0,I0,E0,U_I0,U_E0)=StatEq
        (V_Ess,V_Uss,W_Ess,W_Uss,pol_a_Ess,pol_a_Uss,pol_μ_Uss,pol_σ_Ess,pol_σ_Uss,Jss,θss,Φss,Yss,Iss,Ess,U_Iss,U_Ess)=StatEq
    elseif permanent==0
        (V_Ess,V_Uss,W_Ess,W_Uss,pol_a_Ess,pol_a_Uss,pol_μ_Uss,pol_σ_Ess,pol_σ_Uss,Jss,θss,Φss,Yss,Iss,Ess,U_Iss,U_Ess)=StatEq
        V_E0,V_U0,W_E0,W_U0=V_Ess,V_Uss,W_Ess,W_Uss
        pol_a_E0,pol_a_U0,pol_μ_U0,pol_σ_E0,pol_σ_U0=pol_a_Ess,pol_a_Uss,pol_μ_Uss,pol_σ_Ess,pol_σ_Uss
        J0,θ0,Φ0,Y0,I0,E0,U_I0,U_E0=Jss,θss,Φss,Yss,Iss,Ess,U_Iss,U_Ess
    end


    if Guess==false
        T=shockdur*2
        Iold=zeros(n_i,T)
        Eold=zeros(n_i,T)

        if permanent==1
            r=range(0,1,length=T+1)
            for i_i in 1:n_i
                for t in 1:T
                    Iold[i_i,t]=I0*(1-r[t])+Iss*r[t]
                    Eold[i_i,t]=E0*(1-r[t])+Ess*r[t]
                end
            end
        elseif permanent==0
            for t in 1:T
                Iold[:,t]=I0
                Eold[:,t]=E0
            end

            for i_i in 1:n_i
                if zt[i_i,2]<zt[i_i,1]
                    shocki=i_i
                end
            end
            Iold[shocki,2]=0.7*I0[shocki,2]
            Eold[shocki,2]=0.7*E0[shocki,2]
            r=range(0,1,length=T)
            for t in 3:T
                Iold[shocki,t]=Iold[shocki,2]*(1-r[t-1])+Iss*r[t-1]
                Eold[shocki,t]=Eold[shocki,2]*(1-r[t-1])+Ess*r[t-1]
            end
        end

    else
        (Iold,Eold)=Guess
        T=size(Iold,2)
        if (Iold[:,1]!=I0 || Eold[:,1]!=E0)
            println("code error: Guess for aggregates differs from initial steady state")
            return
        end
    end

    wages(Y,z,i,e)=((1/n_i)*γ*Y^(1/ν)*z^(1-(1/ν))*i^(γ*(1-(1/ν))-1)*e^((1-γ)*(1-(1/ν))),
                        (1/n_i)*(1-γ)*Y^(1/ν)*z^(1-(1/ν))*i^(γ*(1-(1/ν)))*e^((1-γ)*(1-(1/ν))-1))

    wt=zeros(n_i,n_s,T)
    for t in 1:T
        Y=0
        y=zeros(n_i)
        for i_i in 1:n_i
            y[i_i]=z[i_i]*E_I[i_i]^γ*E_E[i_i]^(1-γ)
            Y+=(1/n_i)*y[i_i]^((ν-1)/ν)
        end
        Y=Y^(ν/(ν-1))

        for i_i in 1:n_i
            for s_i in 1:n_s
                if t<shockdur
                    z=zt[i_i,t+1]
                else
                    z=zt[i_i,end]
                end
                wt[i_i,s_i,t]=wages(Y,z,Iold[i_i,t],Eold[i_i,t])
            end
        end
    end


    nstates=length(Φss)
    nstates_E=length(V_Ess)
    nstates_U=length(V_Uss)


    V_E=zeros(nstates_E,T)
    V_U=zeros(nstates_U,T)
    W_E=zeros(nstates_E,T)
    W_U=zeros(nstates_U,T)
    pol_a_E=zeros(Int64,nstates_E,T)
    pol_a_U=zeros(Int64,nstates_U,T)
    pol_μ_U=zeros(Int64,nstates_U,n_i,T)
    pol_σ_E=zeros(nstates_E,T)
    pol_σ_U=zeros(nstates_E,n_i,T)

    J=zeros(nstates_E,T)
    θ=zeros(nstates_E,T)

    Φ=zeros(length(Φss),T)

    I=zeros(n_i,T)
    E=zeros(n_i,T)
    U_I=zeros(T)
    U_E=zeros(n_i,T)

    u(c)=if c>0 c^(1-σ)/(1-σ) else -Inf end
    p(θ)=min(m*θ^(1-ξ),1)

    error=1000

    while error>ϵ

        for t in T:-1:1

            V_E_aux=zeros(n_a,nstates_E)
            V_U_aux=zeros(n_a,nstates_U)

            if t==T
                W_E_old,W_U_old=W_Ess,W_Uss
                V_E_old,V_U_old=V_Ess,V_Uss
            else
                W_E_old,W_U_old=W_E[:,t+1],W_U[:,t+1]
                V_E_old,V_U_old=V_E[:,t+1],V_U[:,t+1]
            end

            # 1) V_E
            for ind in eachindex(V_E[:,t])
                μ_i=statestogrid_E[ind,4]
                a_i=statestogrid_E[ind,3]
                s_i=statestogrid_E[ind,2]
                i_i=statestogrid_E[ind,1]

                wage=wt[i_i,s_i,t]

                #if a_i==1
                    a1_start=1
                #else
                #    a1_start=max(pol_a_E[[i_i-1,s_i-1,(a_i-1)-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]]-1,1)
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
                #    if a1_i>1 && V_E_aux[a1_i,ind]<V_E_aux[a1_i-1,ind]
                #        V_E[ind]=V_E_aux[a1_i-1,ind]
                #        pol_a_E[ind]=a1_i-1
                #        break
                #    elseif a1_i==n_a
                #        V_E[ind]=V_E_aux[a1_i,ind]
                #        pol_a_E[ind]=a1_i
                #    end
                end
                V_E[ind,t],pol_a_E[ind,t]=findmax(V_E_aux[:,ind])
            end

            # 2) V_U
            for ind in eachindex(V_U[:,t])
                a_i=statestogrid_U[ind,3]
                s_i=statestogrid_U[ind,2]
                i_i=statestogrid_U[ind,1]

                #if a_i==1
                    a1_start=1
                #else
                #    a1_start=max(pol_a_U[[i_i-1,s_i-1,(a_i-1)]'*[n_s*n_a,n_a,1]]-1,1)
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
                V_U[ind,t],pol_a_U[ind,t]=findmax(V_U_aux[:,ind])
            end

            if t==T
                pol_σ_E_old=pol_a_Ess
                J_old=Jss
            else
                pol_σ_E_old=pol_σ_E_old[:,t+1]
                J_old=Jss[:,t+1]
            end

            for ind in eachindex(J[:,t])
                μ_i=statestogrid[ind,4]
                a_i=statestogrid[ind,3]
                s_i=statestogrid[ind,2]
                i_i=statestogrid[ind,1]

                wage=wt[i_i,s_i,t]
                a1_i=pol_a_E[ind]

                if s_i==1
                    ind1_0=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    ind1_1=[i_i-1,2-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    σ_prob0=pol_σ_E_old[ind1_0]
                    σ_prob1=pol_σ_E_old[ind1_1]
                    J[ind,t]=(1-grid_μ[μ_i])*wage+β*(1-ρ)*((1-α)*σ_prob0*J_old[ind1_0]+α*σ_prob1*J_old[ind1_1])
                elseif s_i==2
                    ind1=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    σ_prob=pol_σ_E_old[ind1]
                    J[ind,t]=(1-grid_μ[μ_i])*wage+β*(1-ρ)*σ_prob*J_old[ind1]
                end
            end

            for ind in eachindex(θ[:,t])
                θ[ind,t]=q_inv(κ/J[ind,t])
            end

            # 1) W_E
            for ind in eachindex(W_E[:,t])
                μ_i=statestogrid_E[ind,4]
                a_i=statestogrid_E[ind,3]
                s_i=statestogrid_E[ind,2]
                i_i=statestogrid_E[ind,1]

                if s_i==1
                    W_E[ind,t]=σ_ϵ*log(exp(V_E[ind,t]/σ_ϵ)+exp(V_U[a_i,t]/σ_ϵ))
                    pol_σ_E[ind,t]=exp(V_E[ind,t]/σ_ϵ)/(exp(V_E[ind,t]/σ_ϵ)+exp(V_U[a_i,t]/σ_ϵ))
                elseif s_i==2
                    W_E[ind,t]=σ_ϵ*log(exp(V_E[ind,t]/σ_ϵ)+exp(V_U[i_i*n_a+a_i,t]/σ_ϵ))
                    pol_σ_E[ind,t]=exp(V_E[ind,t]/σ_ϵ)/(exp(V_E[ind,t]/σ_ϵ)+exp(V_U[i_i*n_a+a_i,t]/σ_ϵ))
                end
            end

            # 2) W_U
            V_job_s=zeros(nstates_U,n_i)
            for ind in eachindex(W_U[:,t])
                a_i=statestogrid_U[ind,3]
                s_i=statestogrid_U[ind,2]
                i_i=statestogrid_U[ind,1]

                Val_temp=zeros(n_μ)


                if s_i==1
                    for i1_i in 1:n_i
                        for μ1_i in 1:n_μ
                            ste=[i1_i-1,1-1,a_i-1,μ1_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                            prob=p(θ[ste,t])
                            Val_temp[μ1_i]=prob*V_E[ste,t]+(1-prob)*V_U[ind,t]
                        end
                        V_job_s[ind,i1_i],pol_μ_U[ind,i1_i,t]=findmax(Val_temp)
                    end

                    W_U[ind,t]=σ_ϵ*log(sum(exp.(V_job_s[ind,:]/σ_ϵ)))
                    for i1_i in 1:n_i
                        pol_σ_U[ind,i1_i,t]=exp(V_job_s[ind,i1_i]/σ_ϵ)/sum(exp.(V_job_s[ind,:]/σ_ϵ))
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
                            prob=p(θ[ste,t])
                            Val_temp[μ1_i]=prob*V_E[ste,t]+(1-prob)*V_U[ind,t]
                        end
                        V_job_s[ind,i1_i],pol_μ_U[ind,i1_i,t]=findmax(Val_temp)
                    end

                    W_U[ind,t]=σ_ϵ*log(sum(exp.(V_job_s[ind,:]/σ_ϵ)))
                    for i1_i in 1:n_i
                        pol_σ_U[ind,i1_i,t]=exp(V_job_s[ind,i1_i]/σ_ϵ)/sum(exp.(V_job_s[ind,:]/σ_ϵ))
                    end
                end
            end
        end

        Φ[:,1]=Φ0

        for t in 1:T
            # Construct Transition matrix
            Tr=spzeros(nstates,nstates)

            for ind in 1:nstates
                μ_i=statestogrid[ind,5]
                a_i=statestogrid[ind,4]
                s_i=statestogrid[ind,3]
                i_i=statestogrid[ind,2]
                e_i=statestogrid[ind,1]

                if e_i==1

                    indE=[i_i-1,s_i-1,a_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    a1_i=pol_a_E[indE,t]
                    if s_i==1
                        ind1_e0=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        ind1_e1=[i_i-1,2-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        ind1_u0=n_i*n_s*n_a*n_μ+a1_i
                        ind1_u1=n_i*n_s*n_a*n_μ+i_i*n_a+a1_i


                        Tr[ind,ind1_e0]=(1-ρ)*(1-α)*pol_σ_E[ind1_e0,t]
                        Tr[ind,ind1_u0]=(1-ρ)*(1-α)*(1-pol_σ_E[ind1_e0,t])+ρ
                        Tr[ind,ind1_e1]=(1-ρ)*α*pol_σ_E[ind1_e1,t]
                        Tr[ind,ind1_u1]=(1-ρ)*α*(1-pol_σ_E[ind1_e1,t])


                    elseif s_i==2
                        ind1_e=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        ind1_u1=n_i*n_s*n_a*n_μ+i_i*n_a+a1_i
                        ind1_u0=n_i*n_s*n_a*n_μ+a_i

                        Tr[ind,ind1_e]=(1-ρ)*pol_σ_E[ind1_e,t]
                        Tr[ind,ind1_u1]=(1-ρ)*(1-pol_σ_E[ind1_e,t])+ρ-δ
                        Tr[ind,ind1_u0]=δ

                    end


                elseif e_i==2

                    if s_i==1
                        indU=a_i
                    elseif s_i==2
                        indU=i_i*n_a+a_i
                    end

                    ind1_u=zeros(Int64,n_i)
                    ind1_e=zeros(Int64,n_i)

                    a1_i=pol_a_U[indU,t]

                    if s_i==1

                        for i1_i in 1:n_i
                            ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                            ind1_e[i1_i]=[i1_i-1,s_i-1,a1_i-1,pol_μ_U[a1_i,i1_i,t]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]

                            Tr[ind,ind1_e[i1_i]]=pol_σ_U[a1_i,i1_i,t]*p(θ[ind1_e[i1_i],t])
                            Tr[ind,ind1_u[i1_i]]+=pol_σ_U[a1_i,i1_i,t]*(1-p(θ[ind1_e[i1_i],t]))
                        end



                    elseif s_i==2

                        for i1_i in 1:n_i
                            if i1_i==i_i
                                ind1_u[i1_i]=n_i*n_s*n_a*n_μ+i1_i*n_a+a1_i
                                ind1_e[i1_i]=[i1_i-1,s_i-1,a1_i-1,pol_μ_U[i1_i*n_a+a1_i,i1_i,t]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                            else
                                ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                                ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i,t]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                            end

                            Tr[ind,ind1_e[i1_i]]=pol_σ_U[i_i*n_a+a1_i,i1_i,t]*p(θ[ind1_e[i1_i],t])
                            Tr[ind,ind1_u[i1_i]]+=pol_σ_U[i_i*n_a+a1_i,i1_i,t]*(1-p(θ[ind1_e[i1_i],t]))
                        end
                    end
                end
            end

            # Compute the evolution of the distribution
            if t<T
                Φ[:,t+1]=Tr'*Φ[:,t]
            else
                ΦT=Tr'*Φ[:,t]
            end

            # Compute unemployed and employed of each type
            if t<T
                for i_i in 1:n_i
                    I[i_i,t]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+1:(i_i-1)*n_s*n_a*n_μ+n_μ*n_a,t])
                    E[i_i,t]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+n_μ*n_a+1:(i_i-1)*n_s*n_a*n_μ+2*n_μ*n_a,t])
                    U_E[i_i,t]=sum(Φ[n_i*n_s*n_a*n_μ+i_i*n_a+1:n_i*n_s*n_a*n_μ+i_i*n_a+n_a,t])
                end
                U_I[t]=sum(Φ[n_i*n_s*n_a*n_μ+1:n_i*n_s*n_a*n_μ+n_a,t])
            end
        end

        errors=[maximum(abs.(I-Iold));maximum(abs.(E-Eold))]
        error=maximum(errors)

        println("error: ",error)

        Iold=I*wupdate+Iold*(1-wupdate)
        Eold=E*wupdate+Eold*(1-wupdate)
    end

    errorT=maximum(abs.(ΦT-Φss))

    println("error last period: ",errorT)


    return I,E,U_I,U_E
end
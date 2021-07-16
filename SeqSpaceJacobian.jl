# Sequence space Jacobian

function SeqSpaceJacobian(grids,StatEq,w)

    grid_a0=LinRange(grid_a[1],grid_a[end],200)
    grids0=(grid_i,grid_s,grid_a0,grid_μ)

    @eval @everywhere grid_a=$grid_a
    @eval @everywhere grid_a0=$grid_a0
    @eval @everywhere grid_μ=$grid_μ

    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    n_a0=length(grid_a0)

    @eval @everywhere n_i=$n_i
    @eval @everywhere n_s=$n_s
    @eval @everywhere n_a=$n_a
    @eval @everywhere n_μ=$n_μ
    @eval @everywhere n_a0=$n_a0

    (V_Ess,V_Uss,W_Ess,W_Uss,pol_a_Ess,pol_a_Uss,pol_μ_Uss,pol_σ_Ess,pol_σ_Uss,Jss,θss,Φss,Yss,Ess,Uss)=StatEq

    T=200
    dx=1e-6

    nstates=length(Φss)
    nstates_E=length(V_Ess)
    nstates_U=length(V_Uss)
    nsvars_E=4
    ngrids_vars_E=[n_i,n_s,n_a,n_μ]
    nsvars_U=3
    ngrids_vars_U=[n_i,n_s,n_a]

    statestogrid_E=zeros(Int64,nstates_E,nsvars_E)
    for v in 1:nsvars_E
        statestogrid_E[:,v]=kron(ones(prod(ngrids_vars_E[1:v-1]),1),kron(1:ngrids_vars_E[v],ones(prod(ngrids_vars_E[v+1:nsvars_E]),1)))
    end

    statestogrid_U=zeros(Int64,nstates_U,nsvars_U)
    statestogrid_U[1:n_a,:]=hcat(ones(n_a,2),1:n_a)
    for i_i in 1:n_i
        statestogrid_U[n_a+(i_i-1)*(n_s-1)*n_a+1:n_a+i_i*(n_s-1)*n_a,:]=hcat(i_i*ones((n_s-1)*n_a,1),kron(2:n_s,ones(n_a)),kron(ones(n_s-1),1:n_a))
    end

    statestogrid=[hcat(ones(Int64,size(statestogrid_E,1),1),statestogrid_E);hcat(2*ones(Int64,size(statestogrid_U,1),1),statestogrid_U,ones(Int64,size(statestogrid_U,1),1))]

    @eval @everywhere statestogrid_E=$statestogrid_E
    @eval @everywhere statestogrid_U=$statestogrid_U
    @eval @everywhere statestogrid=$statestogrid

    nstates0=n_i*n_s*n_a0*n_μ+n_a0+n_i*n_a0
    nstates_E0=n_i*n_s*n_a0*n_μ
    nstates_U0=n_a0+n_i*(n_s-1)*n_a0
    nsvars_E0=4
    ngrids_vars_E0=[n_i,n_s,n_a0,n_μ]
    nsvars_U0=3
    ngrids_vars_U0=[n_i,n_s,n_a0]

    statestogrid_E0=zeros(Int64,nstates_E0,nsvars_E0)
    for v in 1:nsvars_E0
        statestogrid_E0[:,v]=kron(ones(prod(ngrids_vars_E0[1:v-1]),1),kron(1:ngrids_vars_E0[v],ones(prod(ngrids_vars_E0[v+1:nsvars_E0]),1)))
    end

    statestogrid_U0=zeros(Int64,nstates_U0,nsvars_U0)
    statestogrid_U0[1:n_a0,:]=hcat(ones(n_a0,2),1:n_a0)
    for i_i in 1:n_i
        statestogrid_U0[n_a0+(i_i-1)*(n_s-1)*n_a0+1:n_a0+i_i*(n_s-1)*n_a0,:]=hcat(i_i*ones((n_s-1)*n_a0,1),kron(2:n_s,ones(n_a0)),kron(ones(n_s-1),1:n_a0))
    end

    statestogrid0=[hcat(ones(Int64,size(statestogrid_E0,1),1),statestogrid_E0);hcat(2*ones(Int64,size(statestogrid_U0,1),1),statestogrid_U0,ones(Int64,size(statestogrid_U0,1),1))]

    @eval @everywhere statestogrid_E0=$statestogrid_E0
    @eval @everywhere statestogrid_U0=$statestogrid_U0
    @eval @everywhere statestogrid0=$statestogrid0

    V_E=SharedArray{Float64}(nstates_E,T)
    V_U=SharedArray{Float64}(nstates_U,T)
    W_E=SharedArray{Float64}(nstates_E,T)
    W_U=SharedArray{Float64}(nstates_U,T)
    pol_a_E=SharedArray{Float64}(nstates_E,T)
    pol_a_U=SharedArray{Float64}(nstates_U,T)
    pol_μ_U=SharedArray{Int64}(nstates_U,n_i,T)
    pol_σ_E=SharedArray{Float64}(nstates_E,T)
    pol_σ_U=SharedArray{Float64}(nstates_U,n_i,T)

    pol_a_E=SharedArray{Float64}(nstates_E,T)
    pol_a_U=SharedArray{Float64}(nstates_U,T)

    pol_a_Ei=SharedArray{Int64}(nstates_E,T)
    pol_a_Ui=SharedArray{Int64}(nstates_U,T)

    J=SharedArray{Float64}(nstates_E,T)
    θ=SharedArray{Float64}(nstates_E,T)

    Φ=zeros(length(Φss),T)

    E=[zeros(n_i,n_s) for t=1:T]
    U=[zeros(n_i,n_s) for t=1:T]

    @everywhere u(c)=if c>0 (c^(1-σ)-1)/(1-σ) else -Inf end
    @everywhere p(θ)=min(m*θ^(1-ξ),1)
    @everywhere q_inv(y)=if y>1 0.0 elseif y<0 0.0 else (y/m)^(-1/ξ) end

    Jacobiani=[zeros(T,T) for i in 1:length(E)]
    Jacobian=[Jacobiani for i in 1:length(w)]

    for i in eachindex(w)
        for s in 1:T
            wt=[w for t=1:T]
            wt[s][i]=w[i]+dx

            V_Eaux=SharedArray{Float64}(nstates_E0,T)
            V_Uaux=SharedArray{Float64}(nstates_U0,T)
            W_Eaux=SharedArray{Float64}(nstates_E0,T)
            W_Uaux=SharedArray{Float64}(nstates_U0,T)
            pol_a_Eaux=SharedArray{Float64}(nstates_E0,T)
            pol_a_Uaux=SharedArray{Float64}(nstates_U0,T)
            pol_μ_Uaux=SharedArray{Int64}(nstates_U0,n_i,T)
            pol_σ_Eaux=SharedArray{Float64}(nstates_E0,T)
            pol_σ_Uaux=SharedArray{Float64}(nstates_U0,n_i,T)

            Jaux=SharedArray{Float64}(nstates_E0,T)
            θaux=SharedArray{Float64}(nstates_E0,T)


            for t in T:-1:1

                if t==T
                    W_E_old,W_U_old=W_Ess,W_Uss
                    V_E_old,V_U_old=V_Ess,V_Uss
                    grid_a_aux=grid_a
                    n_a_aux=n_a
                else
                    W_E_old,W_U_old=W_Eaux[:,t+1],W_Uaux[:,t+1]
                    V_E_old,V_U_old=V_Eaux[:,t+1],V_Uaux[:,t+1]
                    grid_a_aux=grid_a0
                    n_a_aux=n_a0
                end

                @eval @everywhere V_E_old=$V_E_old
                @eval @everywhere V_U_old=$V_U_old
                @eval @everywhere W_E_old=$W_E_old
                @eval @everywhere W_U_old=$W_U_old

                # 1) V_E
                @sync @distributed for ind in eachindex(V_Eaux[:,t])
                    μ_i=statestogrid_E0[ind,4]
                    a_i=statestogrid_E0[ind,3]
                    s_i=statestogrid_E0[ind,2]
                    i_i=statestogrid_E0[ind,1]

                    wage=w[i_i,s_i]

                    a_guess=[grid_a0[a_i]+1e-2]

                    if s_i==1
                        interp_V_U=LinearInterpolation(grid_a_aux,V_U_old[1:n_a_aux];extrapolation_bc=Line())
                        ind1_en=[i_i-1,1-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                        interp_W_En=LinearInterpolation(grid_a_aux,W_E_old[ind1_en:n_μ:(ind1_en-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                        ind1_ep=[i_i-1,2-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                        interp_W_Ep=LinearInterpolation(grid_a_aux,W_E_old[ind1_ep:n_μ:(ind1_ep-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())

                        Veval_0(a1)= -(u((1+r)*grid_a0[a_i]+grid_μ[μ_i]*wage-a1[1])+β*(ρ*interp_V_U(a1[1])+(1-ρ)*((1-α[s_i])*interp_W_En(a1[1])+α[s_i]*interp_W_Ep(a1[1]))))
                        if Veval_0(a_min)<Veval_0(a_min+1e-12)
                            pol_a_Eaux[ind,t]=a_min
                            V_Eaux[ind,t]=-Veval_0(a_min)
                        else
                            opt=optimize(Veval_0,a_guess,BFGS())
                            pol_a_Eaux[ind,t]=opt.minimizer[1]
                            V_Eaux[ind,t]=-opt.minimum
                        end
                    elseif s_i<n_s
                        interp_V_Ud=LinearInterpolation(grid_a_aux,V_U_old[1:n_a_aux];extrapolation_bc=Line())
                        interp_V_U=LinearInterpolation(grid_a_aux,V_U_old[n_a_aux+(i_i-1)*(n_s-1)*n_a_aux+(s_i-2)*n_a_aux+1:n_a_aux+(i_i-1)*(n_s-1)*n_a_aux+(s_i-1)*n_a_aux];extrapolation_bc=Line())
                        ind1_en=[i_i-1,s_i-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                        interp_W_En=LinearInterpolation(grid_a_aux,W_E_old[ind1_en:n_μ:(ind1_en-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                        ind1_ep=[i_i-1,(s_i+1)-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                        interp_W_Ep=LinearInterpolation(grid_a_aux,W_E_old[ind1_ep:n_μ:(ind1_ep-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())


                        Veval_1(a1)= -(u((1+r)*grid_a0[a_i]+grid_μ[μ_i]*wage-a1[1])+β*((ρ-δ)*interp_V_U(a1[1])+δ*interp_V_Ud(a1[1])+(1-ρ)*((1-α[s_i])*interp_W_En(a1[1])+α[s_i]*interp_W_Ep(a1[1]))))
                        if Veval_1(a_min)<Veval_1(a_min+1e-12)
                            pol_a_Eaux[ind,t]=a_min
                            V_Eaux[ind,t]=-Veval_1(a_min)
                        else
                            opt=optimize(Veval_1,a_guess,BFGS())
                            pol_a_Eaux[ind,t]=opt.minimizer[1]
                            V_Eaux[ind,t]=-opt.minimum
                        end
                    elseif s_i==n_s
                        interp_V_Ud=LinearInterpolation(grid_a_aux,V_U_old[1:n_a_aux];extrapolation_bc=Line())
                        interp_V_U=LinearInterpolation(grid_a_aux,V_U_old[n_a_aux+(i_i-1)*(n_s-1)*n_a_aux+(s_i-2)*n_a_aux+1:n_a_aux+(i_i-1)*(n_s-1)*n_a_aux+(s_i-1)*n_a_aux];extrapolation_bc=Line())
                        ind1_e=[i_i-1,s_i-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                        interp_W_E=LinearInterpolation(grid_a_aux,W_E_old[ind1_e:n_μ:(ind1_e-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())

                        Veval_2(a1)= -(u((1+r)*grid_a0[a_i]+grid_μ[μ_i]*wage-a1[1])+β*((ρ-δ)*interp_V_U(a1[1])+δ*interp_V_Ud(a1[1])+(1-ρ)*interp_W_E(a1[1])))
                        if Veval_2(a_min)<Veval_2(a_min+1e-12)
                            pol_a_Eaux[ind,t]=a_min
                            V_Eaux[ind,t]=-Veval_2(a_min)
                        else
                            opt=optimize(Veval_2,a_guess,BFGS())
                            pol_a_Eaux[ind,t]=opt.minimizer[1]
                            V_Eaux[ind,t]=-opt.minimum
                        end
                    end
                end


                # 2) V_U
                @sync @distributed for ind in eachindex(V_Uaux[:,t])
                    a_i=statestogrid_U0[ind,3]
                    s_i=statestogrid_U0[ind,2]
                    i_i=statestogrid_U0[ind,1]

                    ub=b*w[i_i,s_i]

                    if s_i==1
                        interp_W_U=LinearInterpolation(grid_a_aux,W_U_old[1:n_a_aux];extrapolation_bc=Line())

                        Veval_0(a1)=-(u((1+r)*grid_a0[a_i]+ub-a1[1])+β*interp_W_U(a1[1]))

                        a_guess=[grid_a0[a_i]+1e-2]

                        if Veval_0(a_min)<Veval_0(a_min+1e-12)
                            pol_a_Uaux[ind,t]=a_min
                            V_Uaux[ind,t]=-Veval_0(a_min)
                        else
                            opt=optimize(Veval_0,a_guess,BFGS())
                            pol_a_Uaux[ind,t]=opt.minimizer[1]
                            V_Uaux[ind,t]=-opt.minimum
                        end
                    else
                        interp_W_Ud=LinearInterpolation(grid_a_aux,W_U_old[1:n_a_aux];extrapolation_bc=Line())
                        if s_i==2
                            interp_W_Ul=LinearInterpolation(grid_a_aux,W_U_old[1:n_a_aux];extrapolation_bc=Line())
                        elseif s_i>2
                            interp_W_Ul=LinearInterpolation(grid_a_aux,W_U_old[n_a_aux+(i_i-1)*(n_s-1)*n_a_aux+(s_i-3)*n_a_aux+1:n_a_aux+(i_i-1)*(n_s-1)*n_a_aux+(s_i-2)*n_a_aux];extrapolation_bc=Line())
                        end
                        interp_W_U=LinearInterpolation(grid_a_aux,W_U_old[n_a_aux+(i_i-1)*(n_s-1)*n_a_aux+(s_i-2)*n_a_aux+1:n_a_aux+(i_i-1)*(n_s-1)*n_a_aux+(s_i-1)*n_a_aux];extrapolation_bc=Line())

                        Veval_1(a1)=-(u((1+r)*grid_a0[a_i]+ub-a1[1])+β*((1-δ-χ[s_i])*interp_W_U(a1[1])+δ*interp_W_Ud(a1[1])+χ[s_i]*interp_W_Ul(a1[1])))

                        a_guess=[grid_a0[a_i]+1e-2]

                        if Veval_1(a_min)<Veval_1(a_min+1e-12)
                            pol_a_Uaux[ind,t]=a_min
                            V_Uaux[ind,t]=-Veval_1(a_min)
                        else
                            opt=optimize(Veval_1,a_guess,BFGS())
                            pol_a_Uaux[ind,t]=opt.minimizer[1]
                            V_Uaux[ind,t]=-opt.minimum
                        end
                    end
                end

                if t==T
                    pol_σ_E_old=pol_σ_Ess
                    J_old=Jss
                else
                    pol_σ_E_old=pol_σ_Eaux[:,t+1]
                    J_old=Jaux[:,t+1]
                end

                @eval @everywhere J_old=$J_old
                @eval @everywhere pol_σ_E_old=$pol_σ_E_old

                @sync @distributed for ind in eachindex(Jaux[:,t])
                    μ_i=statestogrid_E0[ind,4]
                    a_i=statestogrid_E0[ind,3]
                    s_i=statestogrid_E0[ind,2]
                    i_i=statestogrid_E0[ind,1]

                    wage=w[i_i,s_i]
                    a1=pol_a_Eaux[ind,t]

                    if s_i<n_s
                        ind1_n=[i_i-1,s_i-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                        ind1_p=[i_i-1,(s_i+1)-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                        interp_pol_σ_En=LinearInterpolation(grid_a_aux,pol_σ_E_old[ind1_n:n_μ:(ind1_n-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                        interp_pol_σ_Ep=LinearInterpolation(grid_a_aux,pol_σ_E_old[ind1_p:n_μ:(ind1_p-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                        σ_probn=interp_pol_σ_En(a1)
                        σ_probp=interp_pol_σ_Ep(a1)
                        interp_Jn=LinearInterpolation(grid_a_aux,J_old[ind1_n:n_μ:(ind1_n-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                        interp_Jp=LinearInterpolation(grid_a_aux,J_old[ind1_p:n_μ:(ind1_p-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                        Jn=interp_Jn(a1)
                        Jp=interp_Jp(a1)
                        Jaux[ind,t]=wt[t][i_i,s_i]-grid_μ[μ_i]*wage+β*(1-ρ)*((1-α[s_i])*σ_probn*Jn+α[s_i]*σ_probp*Jp)
                    elseif s_i==n_s
                        ind1=[i_i-1,s_i-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                        interp_pol_σ_E=LinearInterpolation(grid_a_aux,pol_σ_E_old[ind1:n_μ:(ind1-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                        σ_prob=interp_pol_σ_E(a1)
                        interp_J=LinearInterpolation(grid_a_aux,J_old[ind1:n_μ:(ind1-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                        Je=interp_J(a1)
                        Jaux[ind,t]=wt[t][i_i,s_i]-grid_μ[μ_i]*wage+β*(1-ρ)*σ_prob*Je
                    end
                end

                @sync @distributed for ind in eachindex(θaux[:,t])
                    s_i=statestogrid_E0[ind,2]

                    θaux[ind,t]=q_inv(κ[s_i]/(Jaux[ind,t]))
                end

                # 1) W_E
                @sync @distributed for ind in eachindex(W_Eaux[:,t])
                    μ_i=statestogrid_E0[ind,4]
                    a_i=statestogrid_E0[ind,3]
                    s_i=statestogrid_E0[ind,2]
                    i_i=statestogrid_E0[ind,1]

                    if s_i==1
                        W_Eaux[ind,t]=σ_ϵ*((V_Eaux[ind,t]/σ_ϵ)+log(1+exp((V_Uaux[a_i,t]-V_Eaux[ind,t])/σ_ϵ)))
                        pol_σ_Eaux[ind,t]=(1+exp((V_Uaux[a_i,t]-V_Eaux[ind,t])/σ_ϵ))^(-1)
                    else
                        W_Eaux[ind,t]=σ_ϵ*((V_Eaux[ind,t]/σ_ϵ)+log(1+exp((V_Uaux[n_a0+(i_i-1)*(n_s-1)*n_a0+(s_i-2)*n_a0+a_i,t]-V_Eaux[ind,t])/σ_ϵ)))
                        pol_σ_Eaux[ind,t]=(1+exp((V_Uaux[n_a0+(i_i-1)*(n_s-1)*n_a0+(s_i-2)*n_a0+a_i,t]-V_Eaux[ind,t])/σ_ϵ))^(-1)
                    end
                end

                # 2) W_U
                V_job_s=SharedArray{Float64}(nstates_U,n_i)
                @sync @distributed for ind in eachindex(W_Uaux[:,t])
                    a_i=statestogrid_U0[ind,3]
                    s_i=statestogrid_U0[ind,2]
                    i_i=statestogrid_U0[ind,1]

                    Val_temp=zeros(n_μ)

                    if s_i==1
                        for i1_i in 1:n_i
                            for μ1_i in 1:n_μ
                                ste=[i1_i-1,1-1,a_i-1,μ1_i]'*[n_s*n_a0*n_μ,n_a0*n_μ,n_μ,1]
                                prob=p(θaux[ste,t])
                                Val_temp[μ1_i]=prob*V_Eaux[ste,t]+(1-prob)*V_Uaux[ind,t]
                            end
                            V_job_s[ind,i1_i],pol_μ_Uaux[ind,i1_i,t]=findmax(Val_temp)
                        end

                        W_Uaux[ind,t]=σ_ϵ*((V_job_s[ind,1]/σ_ϵ)+log(1+sum(exp.((V_job_s[ind,2:end].-V_job_s[ind,1])/σ_ϵ))))
                        for i1_i in 1:n_i
                            pol_σ_Uaux[ind,i1_i,t]=(sum(exp.((V_job_s[ind,:].-V_job_s[ind,i1_i])/σ_ϵ)))^(-1)
                        end
                    else
                        for i1_i in 1:n_i
                            for μ1_i in 1:n_μ
                                if i1_i==i_i
                                    s1_i=s_i
                                    stu=ind
                                else
                                    s1_i=1
                                    stu=ind
                                end
                                ste=[i1_i-1,s1_i-1,a_i-1,μ1_i]'*[n_s*n_a0*n_μ,n_a0*n_μ,n_μ,1]
                                prob=p(θaux[ste,t])
                                Val_temp[μ1_i]=prob*V_Eaux[ste,t]+(1-prob)*V_Uaux[stu,t]
                            end
                            V_job_s[ind,i1_i],pol_μ_Uaux[ind,i1_i,t]=findmax(Val_temp)
                        end

                        W_Uaux[ind,t]=σ_ϵ*((V_job_s[ind,1]/σ_ϵ)+log(1+sum(exp.((V_job_s[ind,2:end].-V_job_s[ind,1])/σ_ϵ))))
                        for i1_i in 1:n_i
                            pol_σ_Uaux[ind,i1_i,t]=(sum(exp.((V_job_s[ind,:].-V_job_s[ind,i1_i])/σ_ϵ)))^(-1)
                        end
                    end
                end
                println("Value function computed for t=",t)
            end


            Φ[:,1]=Φss

            pol_a_Ei=SharedArray{Int64}(nstates_E,T)
            pol_a_Ui=SharedArray{Int64}(nstates_U,T)

            @sync @distributed for t in 1:T
                pol_val_functionsaux=(V_Eaux[:,t],V_Uaux[:,t],W_Eaux[:,t],W_Uaux[:,t],pol_a_Eaux[:,t],pol_a_Uaux[:,t],pol_μ_Uaux[:,:,t],pol_σ_Eaux[:,t],pol_σ_Uaux[:,:,t],Jaux[:,t],θaux[:,t])
                pol_val_functions=transformVPolFunctions(pol_val_functionsaux,grids0,grid_a)
                (V_E[:,t],V_U[:,t],W_E[:,t],W_U[:,t],pol_a_E[:,t],pol_a_U[:,t],pol_μ_U[:,:,t],pol_σ_E[:,t],pol_σ_U[:,:,t],J[:,t],θ[:,t])=pol_val_functions
                pol_a_Ei[:,t],pol_a_Ui[:,t]=transformPola(pol_a_E[:,t],pol_a_U[:,t],grids)
            end


            for t in 1:T-1
                # Construct Transition matrix

                i=Int64[]
                j=Int64[]
                k=Float64[]

                for ind in 1:nstates
                    μ_i=statestogrid[ind,5]
                    a_i=statestogrid[ind,4]
                    s_i=statestogrid[ind,3]
                    i_i=statestogrid[ind,2]
                    e_i=statestogrid[ind,1]

                    if e_i==1

                        indE=[i_i-1,s_i-1,a_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        a1_i=pol_a_Ei[indE,t]
                        if s_i==1
                            ind1_en=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                            ind1_ep=[i_i-1,(s_i+1)-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                            ind1_un=n_i*n_s*n_a*n_μ+a1_i
                            ind1_up=n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-1)*n_a+a1_i

                            push!(i,ind)
                            push!(j,ind1_en)
                            push!(k,(1-ρ)*(1-α[s_i])*pol_σ_E[ind1_en,t+1])

                            push!(i,ind)
                            push!(j,ind1_un)
                            push!(k,(1-ρ)*(1-α[s_i])*(1-pol_σ_E[ind1_en,t+1])+ρ)

                            push!(i,ind)
                            push!(j,ind1_ep)
                            push!(k,(1-ρ)*α[s_i]*pol_σ_E[ind1_ep,t+1])

                            push!(i,ind)
                            push!(j,ind1_up)
                            push!(k,(1-ρ)*α[s_i]*(1-pol_σ_E[ind1_ep,t+1]))

                        elseif s_i<n_s
                            ind1_en=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                            ind1_ep=[i_i-1,(s_i+1)-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                            ind1_un=n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+a1_i
                            ind1_up=n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-1)*n_a+a1_i
                            ind1_ud=n_i*n_s*n_a*n_μ+a1_i

                            push!(i,ind)
                            push!(j,ind1_en)
                            push!(k,(1-ρ)*(1-α[s_i])*pol_σ_E[ind1_en,t+1])

                            push!(i,ind)
                            push!(j,ind1_un)
                            push!(k,(1-ρ)*(1-α[s_i])*(1-pol_σ_E[ind1_en,t+1])+ρ-δ)

                            push!(i,ind)
                            push!(j,ind1_ep)
                            push!(k,(1-ρ)*α[s_i]*pol_σ_E[ind1_ep,t+1])

                            push!(i,ind)
                            push!(j,ind1_up)
                            push!(k,(1-ρ)*α[s_i]*(1-pol_σ_E[ind1_ep,t+1]))

                            push!(i,ind)
                            push!(j,ind1_ud)
                            push!(k,δ)

                        elseif s_i==n_s
                            ind1_e=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                            ind1_un=n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+a1_i
                            ind1_ud=n_i*n_s*n_a*n_μ+a1_i

                            push!(i,ind)
                            push!(j,ind1_e)
                            push!(k,(1-ρ)*pol_σ_E[ind1_e,t+1])

                            push!(i,ind)
                            push!(j,ind1_un)
                            push!(k,(1-ρ)*(1-pol_σ_E[ind1_e,t+1])+ρ-δ)

                            push!(i,ind)
                            push!(j,ind1_ud)
                            push!(k,δ)
                        end


                    elseif e_i==2

                        if s_i==1
                            indU=a_i
                        else
                            indU=n_a+(i_i-1)*(n_s-1)*n_a+a_i
                        end

                        ind1_u=zeros(Int64,n_i)
                        ind1_e=zeros(Int64,n_i)

                        a1_i=pol_a_Ui[indU,t]

                        if s_i==1

                            for i1_i in 1:n_i
                                ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                                ind1_e[i1_i]=[i1_i-1,s_i-1,a1_i-1,pol_μ_U[a1_i,i1_i,t+1]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]

                                push!(i,ind)
                                push!(j,ind1_e[i1_i])
                                push!(k,pol_σ_U[a1_i,i1_i,t]*p(θ[ind1_e[i1_i],t+1]))

                                push!(i,ind)
                                push!(j,ind1_u[i1_i])
                                push!(k,pol_σ_U[a1_i,i1_i,t]*(1-p(θ[ind1_e[i1_i],t+1])))
                            end

                        else

                            for i1_i in 1:n_i
                                if i1_i==i_i
                                    ind1_u[i1_i]=ind
                                    ind1_e[i1_i]=[i1_i-1,s_i-1,a1_i-1,pol_μ_U[n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+a1_i,i1_i,t+1]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                                else
                                    ind1_u[i1_i]=ind
                                    ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i,t+1]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                                end

                                push!(i,ind)
                                push!(j,ind1_e[i1_i])
                                push!(k,(1-δ-χ[s_i])*pol_σ_U[n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+a1_i,i1_i,t+1]*p(θ[ind1_e[i1_i],t+1]))

                                push!(i,ind)
                                push!(j,ind1_u[i1_i])
                                push!(k,(1-δ-χ[s_i])*pol_σ_U[n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+a1_i,i1_i,t+1]*(1-p(θ[ind1_e[i1_i],t+1])))

                                if s_i==2

                                    ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                                    ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i,t+1]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]

                                    push!(i,ind)
                                    push!(j,ind1_e[i1_i])
                                    push!(k,(δ+χ[s_i])*pol_σ_U[a1_i,i1_i,t+1]*p(θ[ind1_e[i1_i],t+1]))

                                    push!(i,ind)
                                    push!(j,ind1_u[i1_i])
                                    push!(k,(δ+χ[s_i])*pol_σ_U[a1_i,i1_i,t+1]*(1-p(θ[ind1_e[i1_i],t+1])))

                                else
                                    ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                                    ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i,t]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]

                                    push!(i,ind)
                                    push!(j,ind1_e[i1_i])
                                    push!(k,δ*pol_σ_U[a1_i,i1_i,t+1]*p(θ[ind1_e[i1_i],t+1]))

                                    push!(i,ind)
                                    push!(j,ind1_u[i1_i])
                                    push!(k,δ*pol_σ_U[a1_i,i1_i,t+1]*(1-p(θ[ind1_e[i1_i],t+1])))

                                    if i1_i==i_i
                                        ind1_u[i1_i]=n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-3)*n_a+a1_i
                                        ind1_e[i1_i]=[i1_i-1,(s_i-1)-1,a1_i-1,pol_μ_U[n_a+(i_i-1)*(n_s-1)*n_a+(s_i-3)*n_a+a1_i,i1_i,t+1]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                                    else
                                        ind1_u[i1_i]=n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-3)*n_a+a1_i
                                        ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i,t+1]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                                    end

                                    push!(i,ind)
                                    push!(j,ind1_e[i1_i])
                                    push!(k,χ[s_i]*pol_σ_U[n_a+(i_i-1)*(n_s-1)*n_a+(s_i-3)*n_a+a1_i,i1_i,t+1]*p(θ[ind1_e[i1_i],t+1]))

                                    push!(i,ind)
                                    push!(j,ind1_u[i1_i])
                                    push!(k,χ[s_i]*pol_σ_U[n_a+(i_i-1)*(n_s-1)*n_a+(s_i-3)*n_a+a1_i,i1_i,t+1]*(1-p(θ[ind1_e[i1_i],t+1])))
                                end
                            end
                        end
                    end
                end

                push!(i,nstates)
                push!(j,nstates)
                push!(k,0.0)

                Tr=sparse(i,j,k)

                # Compute the evolution of the distribution
                if t<T
                    Φ[:,t+1]=Tr'*Φ[:,t]
                else
                    ΦT=Tr'*Φ[:,t]
                end

                # Compute unemployed and employed of each type
                if t<=T
                    for i_i in 1:n_i
                        for s_i in 1:n_s
                            E[t][i_i,s_i]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+(s_i-1)*n_a*n_μ+1:(i_i-1)*n_s*n_a*n_μ+s_i*n_μ*n_a,t])
                            if s_i==1
                                if i_i==1
                                    U[t][i_i,s_i]=sum(Φ[n_i*n_s*n_a*n_μ+1:n_i*n_s*n_a*n_μ+n_a,t])
                                end
                            else
                                U[t][i_i,s_i]=sum(Φ[n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+1:n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-1)*n_a,t])
                            end
                        end
                    end
                end
                println("Distribution computed for t=",t)
            end

            for i_i in 1:n_i
                for s_i in 1:n_s
                    E[T][i_i,s_i]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+(s_i-1)*n_a*n_μ+1:(i_i-1)*n_s*n_a*n_μ+s_i*n_μ*n_a,T])
                    if s_i==1
                        if i_i==1
                            U[T][i_i,s_i]=sum(Φ[n_i*n_s*n_a*n_μ+1:n_i*n_s*n_a*n_μ+n_a,T])
                        end
                    else
                        U[T][i_i,s_i]=sum(Φ[n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+1:n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-1)*n_a,T])
                    end
                end
            end

            for i_i in 1:n_i
                for s_i in 1:n_s
                    o=(i_i-1)*n_s+s_i
                    for t in 1:T
                        Jacobian[i][o][s,t]=(E[t][i_i,s_i]-Ess[i_i,s_i])/dx
                        display(Jacobian[i][o][s,t])
                    end
                end
            end
        end
    end
    return Jacobian
end

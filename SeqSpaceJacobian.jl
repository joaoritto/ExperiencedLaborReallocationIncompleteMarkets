# Sequence space Jacobian

# Fake news algorithm

function SeqSpaceJacobian(grids,StatEq,p,Trss)

    (grid_o,grid_e,grid_a)=grids
    grid_a0=LinRange(grid_a[1],grid_a[end],200)
    grids0=(grid_o,grid_e,grid_a0)

    @eval @everywhere grid_a=$grid_a
    @eval @everywhere grid_a0=$grid_a0


    n_o,n_e,n_a=length(grid_o),length(grid_e),length(grid_a)
    n_a0=length(grid_a0)

    @eval @everywhere n_o=$n_o
    @eval @everywhere n_e=$n_e
    @eval @everywhere n_a=$n_a
    @eval @everywhere n_a0=$n_a0

    (V_Ess,V_Uss,W_Ess,W_Uss,pol_a_Eiss,pol_a_Uiss,pol_σ_Ess,pol_σ_Uss,Jss,θss,Φss,Yss,Ess,Uss)=StatEq

    T=200
    dx=1e-5

    nstates=length(Φss)
    nstates_E=length(V_Ess)
    nstates_U=length(V_Uss)
    nsvars_E=3
    ngrids_vars_E=[n_o,n_e,n_a]
    nsvars_U=3
    ngrids_vars_U=[n_o,n_e,n_a]

    statestogrid_E=zeros(Int64,nstates_E,nsvars_E)
    for v in 1:nsvars_E
        statestogrid_E[:,v]=kron(ones(prod(ngrids_vars_E[1:v-1]),1),kron(1:ngrids_vars_E[v],ones(prod(ngrids_vars_E[v+1:nsvars_E]),1)))
    end

    statestogrid_U=zeros(Int64,nstates_U,nsvars_U)
    statestogrid_U[1:n_a,:]=hcat(ones(n_a,2),1:n_a)
    for o_i in 1:n_o
        statestogrid_U[n_a+(o_i-1)*(n_e-1)*n_a+1:n_a+o_i*(n_e-1)*n_a,:]=hcat(o_i*ones((n_e-1)*n_a,1),kron(2:n_e,ones(n_a)),kron(ones(n_e-1),1:n_a))
    end

    statestogrid=[hcat(ones(Int64,size(statestogrid_E,1),1),statestogrid_E);hcat(2*ones(Int64,size(statestogrid_U,1),1),statestogrid_U)]

    @eval @everywhere statestogrid_E=$statestogrid_E
    @eval @everywhere statestogrid_U=$statestogrid_U
    @eval @everywhere statestogrid=$statestogrid

    nstates0=n_o*n_e*n_a0+n_a0+n_o*n_a0
    nstates_E0=n_o*n_e*n_a0
    nstates_U0=n_a0+n_o*(n_e-1)*n_a0
    nsvars_E0=3
    ngrids_vars_E0=[n_o,n_e,n_a0]
    nsvars_U0=3
    ngrids_vars_U0=[n_o,n_e,n_a0]

    statestogrid_E0=zeros(Int64,nstates_E0,nsvars_E0)
    for v in 1:nsvars_E0
        statestogrid_E0[:,v]=kron(ones(prod(ngrids_vars_E0[1:v-1]),1),kron(1:ngrids_vars_E0[v],ones(prod(ngrids_vars_E0[v+1:nsvars_E0]),1)))
    end

    statestogrid_U0=zeros(Int64,nstates_U0,nsvars_U0)
    statestogrid_U0[1:n_a0,:]=hcat(ones(n_a0,2),1:n_a0)
    for o_i in 1:n_o
        statestogrid_U0[n_a0+(o_i-1)*(n_e-1)*n_a0+1:n_a0+o_i*(n_e-1)*n_a0,:]=hcat(o_i*ones((n_e-1)*n_a0,1),kron(2:n_e,ones(n_a0)),kron(ones(n_e-1),1:n_a0))
    end

    statestogrid0=[hcat(ones(Int64,size(statestogrid_E0,1),1),statestogrid_E0);hcat(2*ones(Int64,size(statestogrid_U0,1),1),statestogrid_U0)]

    @eval @everywhere statestogrid_E0=$statestogrid_E0
    @eval @everywhere statestogrid_U0=$statestogrid_U0
    @eval @everywhere statestogrid0=$statestogrid0

    V_E=SharedArray{Float64}(nstates_E,T)
    V_U=SharedArray{Float64}(nstates_U,T)
    W_E=SharedArray{Float64}(nstates_E,T)
    W_U=SharedArray{Float64}(nstates_U,T)
    pol_a_E=SharedArray{Float64}(nstates_E,T)
    pol_a_U=SharedArray{Float64}(nstates_U,T)
    pol_σ_E=SharedArray{Float64}(nstates_E,T)
    pol_σ_U=SharedArray{Float64}(nstates_U,n_o,T)

    pol_a_E=SharedArray{Float64}(nstates_E,T)
    pol_a_U=SharedArray{Float64}(nstates_U,T)

    pol_a_Ei=SharedArray{Int64}(nstates_E,T)
    pol_a_Ui=SharedArray{Int64}(nstates_U,T)

    J=SharedArray{Float64}(nstates_E,T)
    θ=SharedArray{Float64}(nstates_E,T)

    Φ=zeros(length(Φss),T)

    E=[zeros(n_o,n_e) for t=1:T]
    U=[zeros(n_o,n_e) for t=1:T]

    @everywhere u(c)=if c>0 (c^(1-σ)-1)/(1-σ) else -Inf end
    @everywhere P(θ)=min(m*θ^(1-ξ),1)
    @everywhere q_inv(y)=if y>1 0.0 elseif y<0 0.0 else (y/m)^(-1/ξ) end


    Jacobian=[[zeros(T,T) for i in 1:length(p)] for o in 1:length(E)]

    # Step 1
    Di=[zeros(nstates,T) for i in 1:length(p)]

    Yoi=[[zeros(T) for i in 1:length(p)] for o in 1:length(E)]


    for p_i in 1:length(p)
        pt=[zeros(n_o,n_e) for t in 1:T]
        for t in 1:T
            pt[t][:]=p
        end

        o_i=ceil(Int64,p_i/n_e)
        e_i=p_i-n_e*(o_i-1)

        pt[T][o_i,e_i]=p[o_i,e_i]+dx

        V_Eaux=SharedArray{Float64}(nstates_E0,T)
        V_Uaux=SharedArray{Float64}(nstates_U0,T)
        W_Eaux=SharedArray{Float64}(nstates_E0,T)
        W_Uaux=SharedArray{Float64}(nstates_U0,T)
        pol_a_Eaux=SharedArray{Float64}(nstates_E0,T)
        pol_a_Uaux=SharedArray{Float64}(nstates_U0,T)
        pol_σ_Eaux=SharedArray{Float64}(nstates_E0,T)
        pol_σ_Uaux=SharedArray{Float64}(nstates_U0,n_o,T)

        Jaux=SharedArray{Float64}(nstates_E0,T)
        θaux=SharedArray{Float64}(nstates_E0,T)

        for t in T:-1:1 # u=T-1-t (u as defined in Auclert et al.)

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
                a_i=statestogrid_E0[ind,3]
                e_i=statestogrid_E0[ind,2]
                o_i=statestogrid_E0[ind,1]

                p_eo=p[o_i,e_i]

                a_guess=[grid_a0[a_i]+1e-2]

                if e_i==1
                    interp_V_U=LinearInterpolation(grid_a_aux,V_U_old[1:n_a_aux];extrapolation_bc=Line())
                    ind1_en=[o_i-1,1-1,1]'*[n_e*n_a_aux,n_a_aux,1]
                    interp_W_En=LinearInterpolation(grid_a_aux,W_E_old[ind1_en:(ind1_en-1)+n_a_aux];extrapolation_bc=Line())
                    ind1_ep=[o_i-1,2-1,1]'*[n_e*n_a_aux,n_a_aux,1]
                    interp_W_Ep=LinearInterpolation(grid_a_aux,W_E_old[ind1_ep:(ind1_ep-1)+n_a_aux];extrapolation_bc=Line())

                    Veval_0(a1)= -(u((1+r)*grid_a0[a_i]+φ*p_eo-a1[1])+β*(ρ*interp_V_U(a1[1])+(1-ρ)*((1-α[e_i])*interp_W_En(a1[1])+α[e_i]*interp_W_Ep(a1[1]))))
                    if Veval_0(a_min)<Veval_0(a_min+1e-12)
                        pol_a_Eaux[ind,t]=a_min
                        V_Eaux[ind,t]=-Veval_0(a_min)
                    else
                        opt=optimize(Veval_0,a_guess,BFGS())
                        pol_a_Eaux[ind,t]=opt.minimizer[1]
                        V_Eaux[ind,t]=-opt.minimum
                    end
                elseif e_i<n_e
                    interp_V_Ud=LinearInterpolation(grid_a_aux,V_U_old[1:n_a_aux];extrapolation_bc=Line())
                    interp_V_U=LinearInterpolation(grid_a_aux,V_U_old[n_a_aux+(o_i-1)*(n_e-1)*n_a_aux+(e_i-2)*n_a_aux+1:n_a_aux+(o_i-1)*(n_e-1)*n_a_aux+(e_i-1)*n_a_aux];extrapolation_bc=Line())
                    ind1_en=[o_i-1,e_i-1,1]'*[n_e*n_a_aux,n_a_aux,1]
                    interp_W_En=LinearInterpolation(grid_a_aux,W_E_old[ind1_en:(ind1_en-1)+n_a_aux];extrapolation_bc=Line())
                    ind1_ep=[o_i-1,(e_i+1)-1,1]'*[n_e*n_a_aux,n_a_aux,1]
                    interp_W_Ep=LinearInterpolation(grid_a_aux,W_E_old[ind1_ep:(ind1_ep-1)+n_a_aux];extrapolation_bc=Line())


                    Veval_1(a1)= -(u((1+r)*grid_a0[a_i]+φ*p_eo-a1[1])+β*((ρ-δ)*interp_V_U(a1[1])+δ*interp_V_Ud(a1[1])+(1-ρ)*((1-α[e_i])*interp_W_En(a1[1])+α[e_i]*interp_W_Ep(a1[1]))))
                    if Veval_1(a_min)<Veval_1(a_min+1e-12)
                        pol_a_Eaux[ind,t]=a_min
                        V_Eaux[ind,t]=-Veval_1(a_min)
                    else
                        opt=optimize(Veval_1,a_guess,BFGS())
                        pol_a_Eaux[ind,t]=opt.minimizer[1]
                        V_Eaux[ind,t]=-opt.minimum
                    end
                elseif e_i==n_e
                    interp_V_Ud=LinearInterpolation(grid_a_aux,V_U_old[1:n_a_aux];extrapolation_bc=Line())
                    interp_V_U=LinearInterpolation(grid_a_aux,V_U_old[n_a_aux+(o_i-1)*(n_e-1)*n_a_aux+(e_i-2)*n_a_aux+1:n_a_aux+(o_i-1)*(n_e-1)*n_a_aux+(e_i-1)*n_a_aux];extrapolation_bc=Line())
                    ind1_e=[o_i-1,e_i-1,1]'*[n_e*n_a_aux,n_a_aux,1]
                    interp_W_E=LinearInterpolation(grid_a_aux,W_E_old[ind1_e:(ind1_e-1)+n_a_aux];extrapolation_bc=Line())

                    Veval_2(a1)= -(u((1+r)*grid_a0[a_i]+φ*p_eo-a1[1])+β*((ρ-δ)*interp_V_U(a1[1])+δ*interp_V_Ud(a1[1])+(1-ρ)*interp_W_E(a1[1])))
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
                e_i=statestogrid_U0[ind,2]
                o_i=statestogrid_U0[ind,1]

                ub=b*p[o_i,e_i]

                if e_i==1
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
                    if e_i==2
                        interp_W_Ul=LinearInterpolation(grid_a_aux,W_U_old[1:n_a_aux];extrapolation_bc=Line())
                    elseif e_i>2
                        interp_W_Ul=LinearInterpolation(grid_a_aux,W_U_old[n_a_aux+(o_i-1)*(n_e-1)*n_a_aux+(e_i-3)*n_a_aux+1:n_a_aux+(o_i-1)*(n_e-1)*n_a_aux+(e_i-2)*n_a_aux];extrapolation_bc=Line())
                    end
                    interp_W_U=LinearInterpolation(grid_a_aux,W_U_old[n_a_aux+(o_i-1)*(n_e-1)*n_a_aux+(e_i-2)*n_a_aux+1:n_a_aux+(o_i-1)*(n_e-1)*n_a_aux+(e_i-1)*n_a_aux];extrapolation_bc=Line())

                    Veval_1(a1)=-(u((1+r)*grid_a0[a_i]+ub-a1[1])+β*((1-δ-χ[e_i])*interp_W_U(a1[1])+δ*interp_W_Ud(a1[1])+χ[e_i]*interp_W_Ul(a1[1])))

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
                a_i=statestogrid_E0[ind,3]
                e_i=statestogrid_E0[ind,2]
                o_i=statestogrid_E0[ind,1]

                p_eo=p[o_i,e_i]
                a1=pol_a_Eaux[ind,t]

                if e_i<n_e
                    ind1_n=[o_i-1,e_i-1,1]'*[n_e*n_a_aux,n_a_aux,1]
                    ind1_p=[o_i-1,(e_i+1)-1,1]'*[n_e*n_a_aux,n_a_aux,1]
                    interp_pol_σ_En=LinearInterpolation(grid_a_aux,pol_σ_E_old[ind1_n:(ind1_n-1)+n_a_aux];extrapolation_bc=Line())
                    interp_pol_σ_Ep=LinearInterpolation(grid_a_aux,pol_σ_E_old[ind1_p:(ind1_p-1)+n_a_aux];extrapolation_bc=Line())
                    σ_probn=interp_pol_σ_En(a1)
                    σ_probp=interp_pol_σ_Ep(a1)
                    interp_Jn=LinearInterpolation(grid_a_aux,J_old[ind1_n:(ind1_n-1)+n_a_aux];extrapolation_bc=Line())
                    interp_Jp=LinearInterpolation(grid_a_aux,J_old[ind1_p:(ind1_p-1)+n_a_aux];extrapolation_bc=Line())
                    Jn=interp_Jn(a1)
                    Jp=interp_Jp(a1)
                    Jaux[ind,t]=pt[t][o_i,e_i]-φ*p_eo+β*(1-ρ)*((1-α[e_i])*σ_probn*Jn+α[e_i]*σ_probp*Jp)
                elseif e_i==n_e
                    ind1=[o_i-1,e_i-1,1]'*[n_e*n_a_aux,n_a_aux,1]
                    interp_pol_σ_E=LinearInterpolation(grid_a_aux,pol_σ_E_old[ind1:(ind1-1)+n_a_aux];extrapolation_bc=Line())
                    σ_prob=interp_pol_σ_E(a1)
                    interp_J=LinearInterpolation(grid_a_aux,J_old[ind1:(ind1-1)+n_a_aux];extrapolation_bc=Line())
                    Je=interp_J(a1)
                    Jaux[ind,t]=pt[t][o_i,e_i]-φ*p_eo+β*(1-ρ)*σ_prob*Je
                end
            end

            @sync @distributed for ind in eachindex(θaux[:,t])
                e_i=statestogrid_E0[ind,2]
                θaux[ind,t]=q_inv(κ[e_i]/(Jaux[ind,t]))
            end

            # 1) W_E
            @sync @distributed for ind in eachindex(W_Eaux[:,t])
                a_i=statestogrid_E0[ind,3]
                e_i=statestogrid_E0[ind,2]
                o_i=statestogrid_E0[ind,1]

                if e_i==1
                    W_Eaux[ind,t]=σ_ϵ*((V_Eaux[ind,t]/σ_ϵ)+log(1+exp((V_Uaux[a_i,t]-V_Eaux[ind,t])/σ_ϵ)))
                    pol_σ_Eaux[ind,t]=(1+exp((V_Uaux[a_i,t]-V_Eaux[ind,t])/σ_ϵ))^(-1)
                else
                    W_Eaux[ind,t]=σ_ϵ*((V_Eaux[ind,t]/σ_ϵ)+log(1+exp((V_Uaux[n_a0+(o_i-1)*(n_e-1)*n_a0+(e_i-2)*n_a0+a_i,t]-V_Eaux[ind,t])/σ_ϵ)))
                    pol_σ_Eaux[ind,t]=(1+exp((V_Uaux[n_a0+(o_i-1)*(n_e-1)*n_a0+(e_i-2)*n_a0+a_i,t]-V_Eaux[ind,t])/σ_ϵ))^(-1)
                end
            end

            # 2) W_U
            V_job_s=SharedArray{Float64}(nstates_U,n_o)
            @sync @distributed for ind in eachindex(W_Uaux[:,t])
                a_i=statestogrid_U0[ind,3]
                e_i=statestogrid_U0[ind,2]
                o_i=statestogrid_U0[ind,1]

                if e_i==1
                    for i1_i in 1:n_o
                        ste=[i1_i-1,1-1,a_i]'*[n_e*n_a0,n_a0,1]
                        prob=P(θaux[ste,t])
                        V_job_s[ind,i1_i]=prob*V_Eaux[ste,t]+(1-prob)*V_Uaux[ind,t]
                    end

                    W_Uaux[ind,t]=σ_ϵ*((V_job_s[ind,1]/σ_ϵ)+log(1+sum(exp.((V_job_s[ind,2:end].-V_job_s[ind,1])/σ_ϵ))))
                    for i1_i in 1:n_o
                        pol_σ_Uaux[ind,i1_i,t]=(sum(exp.((V_job_s[ind,:].-V_job_s[ind,i1_i])/σ_ϵ)))^(-1)
                    end
                else
                    for i1_i in 1:n_o
                        if i1_i==o_i
                            s1_i=e_i
                            stu=ind
                        else
                            s1_i=1
                            stu=ind
                        end
                        ste=[i1_i-1,s1_i-1,a_i]'*[n_e*n_a0,n_a0,1]
                        prob=P(θaux[ste,t])
                        V_job_s[ind,i1_i]=prob*V_Eaux[ste,t]+(1-prob)*V_Uaux[stu,t]
                    end

                    W_Uaux[ind,t]=σ_ϵ*((V_job_s[ind,1]/σ_ϵ)+log(1+sum(exp.((V_job_s[ind,2:end].-V_job_s[ind,1])/σ_ϵ))))
                    for i1_i in 1:n_o
                        pol_σ_Uaux[ind,i1_i,t]=(sum(exp.((V_job_s[ind,:].-V_job_s[ind,i1_i])/σ_ϵ)))^(-1)
                    end
                end
            end
            println("Value function computed for t=",t)
        end

        pol_a_Ei=SharedArray{Int64}(nstates_E,T)
        pol_a_Ui=SharedArray{Int64}(nstates_U,T)

        @sync @distributed for t in 1:T
            pol_val_functionsaux=(V_Eaux[:,t],V_Uaux[:,t],W_Eaux[:,t],W_Uaux[:,t],pol_a_Eaux[:,t],pol_a_Uaux[:,t],pol_σ_Eaux[:,t],pol_σ_Uaux[:,:,t],Jaux[:,t],θaux[:,t])
            pol_val_functions=transformVPolFunctions(pol_val_functionsaux,grids0,grid_a)
            (V_E[:,t],V_U[:,t],W_E[:,t],W_U[:,t],pol_a_E[:,t],pol_a_U[:,t],pol_σ_E[:,t],pol_σ_U[:,:,t],J[:,t],θ[:,t])=pol_val_functions
            pol_a_Ei[:,t],pol_a_Ui[:,t]=transformPola(pol_a_E[:,t],pol_a_U[:,t],grids)
        end

        for t in 0:T-1
            # Construct Transition matrix

            i=Int64[]
            j=Int64[]
            k=Float64[]

            for ind in 1:nstates
                a_i=statestogrid[ind,4]
                e_i=statestogrid[ind,3]
                o_i=statestogrid[ind,2]
                s_i=statestogrid[ind,1]

                if s_i==1

                    indE=[o_i-1,e_i-1,a_i]'*[n_e*n_a,n_a,1]
                    if t==0
                        a1_i=pol_a_Eiss[indE]
                    else
                        a1_i=pol_a_Ei[indE,t]
                    end
                    if e_i==1
                        ind1_en=[o_i-1,e_i-1,a1_i]'*[n_e*n_a,n_a,1]
                        ind1_ep=[o_i-1,(e_i+1)-1,a1_i]'*[n_e*n_a,n_a,1]
                        ind1_un=n_o*n_e*n_a+a1_i
                        ind1_up=n_o*n_e*n_a+n_a+(o_i-1)*(n_e-1)*n_a+(e_i-1)*n_a+a1_i

                        push!(i,ind)
                        push!(j,ind1_en)
                        push!(k,(1-ρ)*(1-α[e_i])*pol_σ_E[ind1_en,t+1])

                        push!(i,ind)
                        push!(j,ind1_un)
                        push!(k,(1-ρ)*(1-α[e_i])*(1-pol_σ_E[ind1_en,t+1])+ρ)

                        push!(i,ind)
                        push!(j,ind1_ep)
                        push!(k,(1-ρ)*α[e_i]*pol_σ_E[ind1_ep,t+1])

                        push!(i,ind)
                        push!(j,ind1_up)
                        push!(k,(1-ρ)*α[e_i]*(1-pol_σ_E[ind1_ep,t+1]))

                    elseif e_i<n_e
                        ind1_en=[o_i-1,e_i-1,a1_i]'*[n_e*n_a,n_a,1]
                        ind1_ep=[o_i-1,(e_i+1)-1,a1_i]'*[n_e*n_a,n_a,1]
                        ind1_un=n_o*n_e*n_a+n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+a1_i
                        ind1_up=n_o*n_e*n_a+n_a+(o_i-1)*(n_e-1)*n_a+(e_i-1)*n_a+a1_i
                        ind1_ud=n_o*n_e*n_a+a1_i

                        push!(i,ind)
                        push!(j,ind1_en)
                        push!(k,(1-ρ)*(1-α[e_i])*pol_σ_E[ind1_en,t+1])

                        push!(i,ind)
                        push!(j,ind1_un)
                        push!(k,(1-ρ)*(1-α[e_i])*(1-pol_σ_E[ind1_en,t+1])+ρ-δ)

                        push!(i,ind)
                        push!(j,ind1_ep)
                        push!(k,(1-ρ)*α[e_i]*pol_σ_E[ind1_ep,t+1])

                        push!(i,ind)
                        push!(j,ind1_up)
                        push!(k,(1-ρ)*α[e_i]*(1-pol_σ_E[ind1_ep,t+1]))

                        push!(i,ind)
                        push!(j,ind1_ud)
                        push!(k,δ)

                    elseif e_i==n_e
                        ind1_e=[o_i-1,e_i-1,a1_i]'*[n_e*n_a,n_a,1]
                        ind1_un=n_o*n_e*n_a+n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+a1_i
                        ind1_ud=n_o*n_e*n_a+a1_i

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


                elseif s_i==2

                    if e_i==1
                        indU=a_i
                    else
                        indU=n_a+(o_i-1)*(n_e-1)*n_a+a_i
                    end

                    ind1_u=zeros(Int64,n_o)
                    ind1_e=zeros(Int64,n_o)

                    if t==0
                        a1_i=pol_a_Uiss[indU]
                    else
                        a1_i=pol_a_Ui[indU,t]
                    end

                    if e_i==1

                        for i1_i in 1:n_o
                            ind1_u[i1_i]=n_o*n_e*n_a+a1_i
                            ind1_e[i1_i]=[i1_i-1,e_i-1,a1_i]'*[n_e*n_a,n_a,1]

                            push!(i,ind)
                            push!(j,ind1_e[i1_i])
                            push!(k,pol_σ_U[a1_i,i1_i,t+1]*P(θ[ind1_e[i1_i],t+1]))

                            push!(i,ind)
                            push!(j,ind1_u[i1_i])
                            push!(k,pol_σ_U[a1_i,i1_i,t+1]*(1-P(θ[ind1_e[i1_i],t+1])))
                        end

                    else

                        for i1_i in 1:n_o
                            if i1_i==o_i
                                ind1_u[i1_i]=ind
                                ind1_e[i1_i]=[i1_i-1,e_i-1,a1_i]'*[n_e*n_a,n_a,1]
                            else
                                ind1_u[i1_i]=ind
                                ind1_e[i1_i]=[i1_i-1,1-1,a1_i]'*[n_e*n_a,n_a,1]
                            end

                            push!(i,ind)
                            push!(j,ind1_e[i1_i])
                            push!(k,(1-δ-χ[e_i])*pol_σ_U[n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+a1_i,i1_i,t+1]*P(θ[ind1_e[i1_i],t+1]))

                            push!(i,ind)
                            push!(j,ind1_u[i1_i])
                            push!(k,(1-δ-χ[e_i])*pol_σ_U[n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+a1_i,i1_i,t+1]*(1-P(θ[ind1_e[i1_i],t+1])))

                            if e_i==2

                                ind1_u[i1_i]=n_o*n_e*n_a+a1_i
                                ind1_e[i1_i]=[i1_i-1,1-1,a1_i]'*[n_e*n_a,n_a,1]

                                push!(i,ind)
                                push!(j,ind1_e[i1_i])
                                push!(k,(δ+χ[e_i])*pol_σ_U[a1_i,i1_i,t+1]*P(θ[ind1_e[i1_i],t+1]))

                                push!(i,ind)
                                push!(j,ind1_u[i1_i])
                                push!(k,(δ+χ[e_i])*pol_σ_U[a1_i,i1_i,t+1]*(1-P(θ[ind1_e[i1_i],t+1])))

                            else
                                ind1_u[i1_i]=n_o*n_e*n_a+a1_i
                                ind1_e[i1_i]=[i1_i-1,1-1,a1_i]'*[n_e*n_a,n_a,1]

                                push!(i,ind)
                                push!(j,ind1_e[i1_i])
                                push!(k,δ*pol_σ_U[a1_i,i1_i,t+1]*P(θ[ind1_e[i1_i],t+1]))

                                push!(i,ind)
                                push!(j,ind1_u[i1_i])
                                push!(k,δ*pol_σ_U[a1_i,i1_i,t+1]*(1-P(θ[ind1_e[i1_i],t+1])))

                                if i1_i==o_i
                                    ind1_u[i1_i]=n_o*n_e*n_a+n_a+(o_i-1)*(n_e-1)*n_a+(e_i-3)*n_a+a1_i
                                    ind1_e[i1_i]=[i1_i-1,(e_i-1)-1,a1_i]'*[n_e*n_a,n_a,1]
                                else
                                    ind1_u[i1_i]=n_o*n_e*n_a+n_a+(o_i-1)*(n_e-1)*n_a+(e_i-3)*n_a+a1_i
                                    ind1_e[i1_i]=[i1_i-1,1-1,a1_i]'*[n_e*n_a,n_a,1]
                                end

                                push!(i,ind)
                                push!(j,ind1_e[i1_i])
                                push!(k,χ[e_i]*pol_σ_U[n_a+(o_i-1)*(n_e-1)*n_a+(e_i-3)*n_a+a1_i,i1_i,t+1]*P(θ[ind1_e[i1_i],t+1]))

                                push!(i,ind)
                                push!(j,ind1_u[i1_i])
                                push!(k,χ[e_i]*pol_σ_U[n_a+(o_i-1)*(n_e-1)*n_a+(e_i-3)*n_a+a1_i,i1_i,t+1]*(1-P(θ[ind1_e[i1_i],t+1])))
                            end
                        end
                    end
                end
            end

            push!(i,nstates)
            push!(j,nstates)
            push!(k,0.0)

            Tr=sparse(i,j,k)

            Di[p_i][:,T-t]=(Tr'*Φss-Φss)/dx

            for o_i in 1:n_o
                for e_i in 1:n_e
                    o=(o_i-1)*n_e+e_i
                    A=zeros(nstates)
                    A[(o_i-1)*n_e*n_a+(e_i-1)*n_a+1:(o_i-1)*n_e*n_a+e_i*n_a].=1
                    Yoi[o][p_i][T-t]=((Tr*A)'*Φss-Ess[o_i,e_i])/dx
                end
            end
        end
    end

    # Step 2
    Pmatrix=[zeros(nstates,T-1) for o in 1:length(E)]
    for o_i in 1:n_o
        for e_i in 1:n_e
            o=(o_i-1)*n_e+e_i
            A=zeros(nstates)
            A[(o_i-1)*n_e*n_a+(e_i-1)*n_a+1:(o_i-1)*n_e*n_a+e_i*n_a].=1
            Pmatrix[o][:,1]=Trss*A
            for u=1:T-2
                Pmatrix[o][:,u+1]=Trss*Pmatrix[o][:,u]
            end
        end
    end

    # Step 3
    Foi=[[zeros(T,T) for i in 1:length(p)] for o in 1:length(E)]

    for o in 1:length(E)
        for i in 1:length(p)
            Foi[o][i][1,:]=Yoi[o][i]
            Foi[o][i][2:end,:]=Pmatrix[o]'*Di[i]
        end
    end

    # Step 4
    for o in 1:length(E)
        for i in 1:length(p)
            for t in 0:T-1
                for s in 0:T-1
                    if t==0 || s==0
                        Jacobian[o][i][t+1,s+1]=Foi[o][i][t+1,s+1]
                    else
                        Jacobian[o][i][t+1,s+1]=Jacobian[o][i][t,s]+Foi[o][i][t+1,s+1]
                    end
                end
            end
        end
    end

    return Jacobian
end


function ComputingdU(grids,StatEq,Jacobian,dZ)

    Y,E,U=StatEq[12],StatEq[13],StatEq[14]
    (grid_o,grid_e,grid_a)=grids
    n_o,n_e,n_a=length(grid_o),length(grid_e),length(grid_a)

    # derivative1U is the derivative of dp_{e1,o1}/dL_{e2,o2} with e1!=e2 and o1!=o2
    # derivative2U is the derivative of dp_{e1,o1}/dL_{e2,o2} with e1!=e2 and o1=o2
    # derivative3U is the derivative of dp_{e1,o1}/dL_{e2,o2} with e1=e2 and o1=o2
    # derivative1Z is the derivative of dp_{e1,o1}/dz_{o2} with o1!=o2
    # derivative2Z is the derivative of dp_{e1,o1}/dz_{o2} with o1=o2

    derivative1U(ϕ,Y,z,E,o_1,o_2,e_1,e_2)=ϕ[o_1]*γ[e_1]*z[o_1]^(1-(1/ν))*E[o_1,e_1]^(ω-1)*sum(γ.*(E[o_1,:].^ω))^((1/ω)*(1-(1/ν))-1)*
    (1/ν)*Y^((1/ν)-1)*(Y^(1/ν)*(z[o_2]*sum(γ.*(E[o_2,:].^ω))^(1/ω))^(-1/ν)*sum(γ.*(E[o_2,:].^ω))^((1/ω)-1)*ϕ[o_2]*
    z[o_2]*γ[e_2]*ω*E[o_2,e_2]^(ω-1))
    derivative2U(ϕ,Y,z,E,o_1,o_2,e_1,e_2)=derivative1U(ϕ,Y,z,E,o_1,o_2,e_1,e_2)+ϕ[o_1]*γ[e_1]*Y^(1/ν)*z[o_1]^(1-(1/ν))*E[o_1,e_1]^(ω-1)*
    ((1/ω)*(1-(1/ν))-1)*sum(γ.*(E[o_2,:].^ω))^((1/ω)*(1-(1/ν))-2)*γ[e_2]*ω*E[o_2,e_2]^(ω-1)
    derivative3U(ϕ,Y,z,E,o_1,o_2,e_1,e_2)=derivative2U(ϕ,Y,z,E,o_1,o_2,e_1,e_2)+ϕ[o_1]*γ[e_1]*Y^(1/ν)*z[o_1]^(1-(1/ν))*
    sum(γ.*(E[o_1,:].^ω))^((1/ω)*(1-(1/ν))-1)*(ω-1)*E[o_2,e_2]^(ω-2)

    derivative1Z(ϕ,Y,z,E,o_1,o_2,e_1)=ϕ[o_1]*γ[e_1]*z[o_1]^(1-(1/ν))*E[o_1,e_1]^(ω-1)*sum(γ.*(E[o_1,:].^ω))^((1/ω)*(1-(1/ν))-1)*(1/ν)*
    Y^((1/ν)-1)*(Y^(1/ν)*(z[o_2]*sum(γ.*(E[o_2,:].^ω))^(1/ω))^(-1/ν)*sum(γ.*(E[o_2,:].^ω))^(1/ω)*ϕ[o_2])
    derivative2Z(ϕ,Y,z,E,o_1,o_2,e_1)=derivative1Z(ϕ,Y,z,E,o_1,o_2,e_1)+ϕ[o_1]*γ[e_1]*Y^(1/ν)*E[o_1,e_1]^(ω-1)*
    sum(γ.*(E[o_1,:].^ω))^((1/ω)*(1-(1/ν))-1)*(1-(1/ν))*z[o_2]^(-1/ν)

    T=size(dZ)[2]
    n_u=n_o*n_e
    H_U=zeros(n_u*T,n_u*T)
    H_Z=zeros(n_u*T,n_o*T)
    for t in 1:T
        for s in 1:T
            for o in 1:n_u
                aux1=zeros(n_u)
                for j in 1:n_u
                    aux1[j]=Jacobian[o][j][t,s]
                end
                for i in 1:n_u
                    o_i=ceil(Int64,i/n_e)
                    e_i=i-n_e*(o_i-1)
                    aux2=zeros(n_u)
                    for j in 1:n_u
                        o_j=ceil(Int64,j/n_e)
                        e_j=j-n_e*(o_j-1)
                        if o_j==o_i
                            if e_j==e_i
                                aux2[j]=derivative3U(ϕ,Y,z,E,o_j,o_i,e_j,e_i)
                            else
                                aux2[j]=derivative2U(ϕ,Y,z,E,o_j,o_i,e_j,e_i)
                            end
                        else
                            aux2[j]=derivative1U(ϕ,Y,z,E,o_j,o_i,e_j,e_i)
                        end
                    end
                    if s==t && o==i
                        H_U[(t-1)*n_u+o,(s-1)*n_u+i]=aux1'*aux2-1
                    else
                        H_U[(t-1)*n_u+o,(s-1)*n_u+i]=aux1'*aux2
                    end
                end
                for o_i in 1:n_o
                    aux2=zeros(n_u)
                    for j in 1:n_u
                        o_j=ceil(Int64,j/n_e)
                        e_j=j-n_e*(o_j-1)
                        if o_j==o_i
                            aux2[j]=derivative2Z(ϕ,Y,z,E,o_j,o_i,e_j)
                        else
                            aux2[j]=derivative1Z(ϕ,Y,z,E,o_j,o_i,e_j)
                        end
                    end
                    H_Z[(t-1)*n_u+o,(s-1)*n_o+o_i]=aux1'*aux2
                end
            end
        end
    end

    dZ=dZ[:]

    dUdZ=-H_U\H_Z
    dU=dUdZ*dZ

    return dU,dUdZ
end

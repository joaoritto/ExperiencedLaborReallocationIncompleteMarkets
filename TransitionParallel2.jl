# Computing a transition

function Transition(grids,StatEq,zt;Guess=false,permanent=0,i_shock=1)

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

    shockdur=size(zt,2)
    ϵ=5e-4
    wupdate=0.005


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
        T=250 #shockdur*2
        Iold=Iss*ones(1,T)
        Eold=Ess*ones(1,T)

        if permanent==1
            step_aux=range(0,1,length=T+1)
            for i_i in 1:n_i
                for t in 1:T
                    Iold[i_i,t]=I0*(1-step_aux[t])+Iss*step_aux[t]
                    Eold[i_i,t]=E0*(1-step_aux[t])+Ess*step_aux[t]
                end
            end
        elseif permanent==0
            for t in 1:T
                Iold[:,t]=I0
                Eold[:,t]=E0
            end


            Iold[i_shock,2]=0.97*I0[i_shock]
            Eold[i_shock,2]=0.97*E0[i_shock]
            step_aux=range(0,1,length=min(ceil(Int64,shockdur*1.5),T))
            for t in 3:length(step_aux) #min(ceil(Int64,shockdur*1.5),T)
                Iold[i_shock,t]=Iold[i_shock,2]*(1-step_aux[t-1])+Iss[i_shock]*step_aux[t-1]
                Eold[i_shock,t]=Eold[i_shock,2]*(1-step_aux[t-1])+Ess[i_shock]*step_aux[t-1]
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
        statestogrid_U[n_a+(i_i-1)*n_a+1:n_a+i_i*n_a,:]=hcat(i_i*ones(n_a,1),2*ones(n_a,1),1:n_a)
    end

    statestogrid=[hcat(ones(Int64,size(statestogrid_E,1),1),statestogrid_E);hcat(2*ones(Int64,size(statestogrid_U,1),1),statestogrid_U,ones(Int64,size(statestogrid_U,1),1))]

    @eval @everywhere statestogrid_E=$statestogrid_E
    @eval @everywhere statestogrid_U=$statestogrid_U
    @eval @everywhere statestogrid=$statestogrid

    nstates0=n_i*n_s*n_a0*n_μ+n_a0+n_i*n_a0
    nstates_E0=n_i*n_s*n_a0*n_μ
    nstates_U0=n_a0+n_i*n_a0
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
        statestogrid_U0[n_a0+(i_i-1)*n_a0+1:n_a0+i_i*n_a0,:]=hcat(i_i*ones(n_a0,1),2*ones(n_a0,1),1:n_a0)
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

    I=zeros(n_i,T)
    E=zeros(n_i,T)
    U_I=zeros(T)
    U_E=zeros(n_i,T)

    @everywhere u(c)=if c>0 (c^(1-σ)-1)/(1-σ) else -Inf end
    @everywhere p(θ)=min(m*θ^(1-ξ),1)
    @everywhere q_inv(y)=if y>1 0.0 elseif y<0 0.0 else (y/m)^(-1/ξ) end

    error=1000

    p0=plot(1:T,Iold')
    plot!(1:T,Eold')
    display(p0)
    ΦT=zeros(length(Φss))
    iter=0

    while error>ϵ
        iter+=1

        wt=zeros(n_i,n_s,T)
        for t in 1:T
            Y=0
            y=zeros(n_i)
            for i_i in 1:n_i
                if t<shockdur
                    z=zt[i_i,t]
                else
                    z=zt[i_i,end]
                end
                y[i_i]=z*Iold[i_i,t]^γ*Eold[i_i,t]^(1-γ)
                Y+=(1/n_i)*y[i_i]^((ν-1)/ν)
            end
            Y=Y^(ν/(ν-1))

            for i_i in 1:n_i
                for s_i in 1:n_s
                    if t<shockdur
                        z=zt[i_i,t]
                    else
                        z=zt[i_i,end]
                    end
                    wt[i_i,s_i,t]=wages(Y,z,Iold[i_i,t],Eold[i_i,t])[s_i]
                end
            end
        end

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

                wage=wt[i_i,s_i,t]

                a_guess=[grid_a0[a_i]+1e-2]

                if s_i==1
                    interp_V_U=LinearInterpolation(grid_a_aux,V_U_old[1:n_a_aux];extrapolation_bc=Line())
                    ind1_ei=[i_i-1,1-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                    interp_W_Ei=LinearInterpolation(grid_a_aux,W_E_old[ind1_ei:n_μ:(ind1_ei-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                    ind1_ee=[i_i-1,2-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                    interp_W_Ee=LinearInterpolation(grid_a_aux,W_E_old[ind1_ee:n_μ:(ind1_ee-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())

                    Veval_i(a1)= -(u((1+r)*grid_a0[a_i]+grid_μ[μ_i]*wage-a1[1])+β*(ρ*interp_V_U(a1[1])+(1-ρ)*((1-α)*interp_W_Ei(a1[1])+α*interp_W_Ee(a1[1]))))
                    if Veval_i(a_min)<Veval_i(a_min+1e-12)
                        pol_a_Eaux[ind,t]=a_min
                        V_Eaux[ind,t]=-Veval_i(a_min)
                    else
                        opt=optimize(Veval_i,a_guess,BFGS())
                        pol_a_Eaux[ind,t]=opt.minimizer[1]
                        V_Eaux[ind,t]=-opt.minimum
                    end
                elseif s_i==2
                    interp_V_Ui=LinearInterpolation(grid_a_aux,V_U_old[1:n_a_aux];extrapolation_bc=Line())
                    interp_V_Ue=LinearInterpolation(grid_a_aux,V_U_old[i_i*n_a_aux+1:(i_i+1)*n_a_aux];extrapolation_bc=Line())
                    ind1_e=[i_i-1,2-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                    interp_W_E=LinearInterpolation(grid_a_aux,W_E_old[ind1_e:n_μ:(ind1_e-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())

                    Veval_e(a1)= -(u((1+r)*grid_a0[a_i]+grid_μ[μ_i]*wage-a1[1])+β*((ρ-δ)*interp_V_Ue(a1[1])+δ*interp_V_Ui(a1[1])+(1-ρ)*interp_W_E(a1[1])))
                    if Veval_e(a_min)<Veval_e(a_min+1e-12)
                        pol_a_Eaux[ind,t]=a_min
                        V_Eaux[ind,t]=-Veval_e(a_min)
                    else
                        opt=optimize(Veval_e,a_guess,BFGS())
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

                ub=b*wt[2,s_i,end]

                if s_i==1
                    interp_W_U=LinearInterpolation(grid_a_aux,W_U_old[1:n_a_aux];extrapolation_bc=Line())

                    Veval_i(a1)=-(u((1+r)*grid_a0[a_i]+ub-a1[1])+β*interp_W_U(a1[1]))

                    a_guess=[grid_a0[a_i]+1e-2]

                    if Veval_i(a_min)<Veval_i(a_min+1e-12)
                        pol_a_Uaux[ind,t]=a_min
                        V_Uaux[ind,t]=-Veval_i(a_min)
                    else
                        opt=optimize(Veval_i,a_guess,BFGS())
                        pol_a_Uaux[ind,t]=opt.minimizer[1]
                        V_Uaux[ind,t]=-opt.minimum
                    end
                elseif s_i==2
                    interp_W_U=LinearInterpolation(grid_a_aux,W_U_old[i_i*n_a_aux+1:(i_i+1)*n_a_aux];extrapolation_bc=Line())
                    interp_W_Ui=LinearInterpolation(grid_a_aux,W_U_old[1:n_a_aux];extrapolation_bc=Line())


                    Veval_e(a1)=-(u((1+r)*grid_a0[a_i]+ub-a1[1])+β*((1-δ)*interp_W_U(a1[1])+δ*interp_W_Ui(a1[1])))

                    a_guess=[grid_a0[a_i]+1e-2]

                    if Veval_e(a_min)<Veval_e(a_min+1e-12)
                        pol_a_Uaux[ind,t]=a_min
                        V_Uaux[ind,t]=-Veval_e(a_min)
                    else
                        opt=optimize(Veval_e,a_guess,BFGS())
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

                wage=wt[i_i,s_i,t]
                a1=pol_a_Eaux[ind,t]

                if s_i==1
                    ind1_i=[i_i-1,s_i-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                    ind1_e=[i_i-1,2-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                    interp_pol_σ_Ei=LinearInterpolation(grid_a_aux,pol_σ_E_old[ind1_i:n_μ:(ind1_i-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                    interp_pol_σ_Ee=LinearInterpolation(grid_a_aux,pol_σ_E_old[ind1_e:n_μ:(ind1_e-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                    σ_probi=interp_pol_σ_Ei(a1)
                    σ_probe=interp_pol_σ_Ee(a1)
                    interp_Ji=LinearInterpolation(grid_a_aux,J_old[ind1_i:n_μ:(ind1_i-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                    interp_Je=LinearInterpolation(grid_a_aux,J_old[ind1_e:n_μ:(ind1_e-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                    Ji=interp_Ji(a1)
                    Je=interp_Je(a1)
                    Jaux[ind,t]=(1-grid_μ[μ_i])*wage+β*(1-ρ)*((1-α)*σ_probi*Ji+α*σ_probe*Je)
                elseif s_i==2
                    ind1=[i_i-1,s_i-1,1-1,μ_i]'*[n_s*n_a_aux*n_μ,n_a_aux*n_μ,n_μ,1]
                    interp_pol_σ_E=LinearInterpolation(grid_a_aux,pol_σ_E_old[ind1:n_μ:(ind1-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                    σ_prob=interp_pol_σ_E(a1)
                    interp_J=LinearInterpolation(grid_a_aux,J_old[ind1:n_μ:(ind1-μ_i)+n_a_aux*n_μ];extrapolation_bc=Line())
                    Je=interp_J(a1)
                    Jaux[ind,t]=(1-grid_μ[μ_i])*wage+β*(1-ρ)*σ_prob*Je
                end
            end

            @sync @distributed for ind in eachindex(θaux[:,t])
                μ_i=statestogrid_E0[ind,4]
                a_i=statestogrid_E0[ind,3]
                s_i=statestogrid_E0[ind,2]
                i_i=statestogrid_E0[ind,1]
                if s_i==1
                    ind1_i=[i_i-1,s_i-1,1-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    interp_Ji=LinearInterpolation(grid_a,Jss[ind1_i:n_μ:(ind1_i-μ_i)+n_a*n_μ];extrapolation_bc=Line())
                    Ji=interp_Ji(grid_a0[a_i])
                    κ=κ_i
                    F=max(F_i*Ji,0.8)
                else
                    ind1=[i_i-1,s_i-1,1-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    interp_J=LinearInterpolation(grid_a,Jss[ind1:n_μ:(ind1-μ_i)+n_a*n_μ];extrapolation_bc=Line())
                    Je=interp_J(grid_a0[a_i])
                    κ=κ_e
                    F=max(F_e*Je,0.8)
                end
                θaux[ind,t]=q_inv(κ/(Jaux[ind,t]-F))
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
                elseif s_i==2
                    W_Eaux[ind,t]=σ_ϵ*((V_Eaux[ind,t]/σ_ϵ)+log(1+exp((V_Uaux[i_i*n_a0+a_i,t]-V_Eaux[ind,t])/σ_ϵ)))
                    pol_σ_Eaux[ind,t]=(1+exp((V_Uaux[i_i*n_a0+a_i,t]-V_Eaux[ind,t])/σ_ϵ))^(-1)
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
                elseif s_i==2
                    for i1_i in 1:n_i
                        for μ1_i in 1:n_μ
                            if i1_i==i_i
                                s1_i=2
                                stu=i_i*n_a0+a_i
                            else
                                s1_i=1
                                stu=i_i*n_a0+a_i #a_i
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


        Φ[:,1]=Φ0

        pol_a_Ei=SharedArray{Int64}(nstates_E,T)
        pol_a_Ui=SharedArray{Int64}(nstates_U,T)

        @sync @distributed for t in 1:T
            pol_val_functionsaux=(V_Eaux[:,t],V_Uaux[:,t],W_Eaux[:,t],W_Uaux[:,t],pol_a_Eaux[:,t],pol_a_Uaux[:,t],pol_μ_Uaux[:,:,t],pol_σ_Eaux[:,t],pol_σ_Uaux[:,:,t],Jaux[:,t],θaux[:,t])
            pol_val_functions=transformVPolFunctions(pol_val_functionsaux,grids0,grid_a)
            (V_E[:,t],V_U[:,t],W_E[:,t],W_U[:,t],pol_a_E[:,t],pol_a_U[:,t],pol_μ_U[:,:,t],pol_σ_E[:,t],pol_σ_U[:,:,t],J[:,t],θ[:,t])=pol_val_functions
            pol_a_Ei[:,t],pol_a_Ui[:,t]=transformPola(pol_a_E[:,t],pol_a_U[:,t],grids)
        end



        for t in 1:T
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
                        ind1_e0=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        ind1_e1=[i_i-1,2-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        ind1_u0=n_i*n_s*n_a*n_μ+a1_i
                        ind1_u1=n_i*n_s*n_a*n_μ+i_i*n_a+a1_i

                        push!(i,ind)
                        push!(j,ind1_e0)
                        push!(k,(1-ρ)*(1-α)*pol_σ_E[ind1_e0,t])

                        push!(i,ind)
                        push!(j,ind1_u0)
                        push!(k,(1-ρ)*(1-α)*(1-pol_σ_E[ind1_e0,t])+ρ)

                        push!(i,ind)
                        push!(j,ind1_e1)
                        push!(k,(1-ρ)*α*pol_σ_E[ind1_e1,t])

                        push!(i,ind)
                        push!(j,ind1_u1)
                        push!(k,(1-ρ)*α*(1-pol_σ_E[ind1_e1,t]))


                    elseif s_i==2
                        ind1_e=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        ind1_u1=n_i*n_s*n_a*n_μ+i_i*n_a+a1_i
                        ind1_u0=n_i*n_s*n_a*n_μ+a_i

                        push!(i,ind)
                        push!(j,ind1_e)
                        push!(k,(1-ρ)*pol_σ_E[ind1_e,t])

                        push!(i,ind)
                        push!(j,ind1_u1)
                        push!(k,(1-ρ)*(1-pol_σ_E[ind1_e,t])+ρ-δ)

                        push!(i,ind)
                        push!(j,ind1_u0)
                        push!(k,δ)

                    end

                elseif e_i==2

                    if s_i==1
                        indU=a_i
                    elseif s_i==2
                        indU=i_i*n_a+a_i
                    end

                    ind1_u=zeros(Int64,n_i)
                    ind1_e=zeros(Int64,n_i)

                    a1_i=pol_a_Ui[indU,t]

                    if s_i==1

                        for i1_i in 1:n_i
                            ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                            ind1_e[i1_i]=[i1_i-1,s_i-1,a1_i-1,pol_μ_U[a1_i,i1_i,t]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]

                            push!(i,ind)
                            push!(j,ind1_e[i1_i])
                            push!(k,pol_σ_U[a1_i,i1_i,t]*p(θ[ind1_e[i1_i],t]))

                            push!(i,ind)
                            push!(j,ind1_u[i1_i])
                            push!(k,pol_σ_U[a1_i,i1_i,t]*(1-p(θ[ind1_e[i1_i],t])))

                        end

                    elseif s_i==2

                        for i1_i in 1:n_i
                            if i1_i==i_i
                                ind1_u[i1_i]=n_i*n_s*n_a*n_μ+i1_i*n_a+a1_i
                                ind1_e[i1_i]=[i1_i-1,s_i-1,a1_i-1,pol_μ_U[i1_i*n_a+a1_i,i1_i,t]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                            else
                                ind1_u[i1_i]=n_i*n_s*n_a*n_μ+i_i*n_a+a1_i #a1_i
                                ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i,t]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                            end

                            push!(i,ind)
                            push!(j,ind1_e[i1_i])
                            push!(k,(1-δ)*pol_σ_U[i_i*n_a+a1_i,i1_i,t]*p(θ[ind1_e[i1_i],t]))

                            push!(i,ind)
                            push!(j,ind1_u[i1_i])
                            push!(k,(1-δ)*pol_σ_U[i_i*n_a+a1_i,i1_i,t]*(1-p(θ[ind1_e[i1_i],t])))


                            ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                            ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i,t]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]

                            push!(i,ind)
                            push!(j,ind1_e[i1_i])
                            push!(k,δ*pol_σ_U[a1_i,i1_i,t]*p(θ[ind1_e[i1_i],t]))

                            push!(i,ind)
                            push!(j,ind1_u[i1_i])
                            push!(k,δ*pol_σ_U[a1_i,i1_i,t]*(1-p(θ[ind1_e[i1_i],t])))

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
                    I[i_i,t]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+1:(i_i-1)*n_s*n_a*n_μ+n_μ*n_a,t])
                    E[i_i,t]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+n_μ*n_a+1:(i_i-1)*n_s*n_a*n_μ+2*n_μ*n_a,t])
                    U_E[i_i,t]=sum(Φ[n_i*n_s*n_a*n_μ+i_i*n_a+1:n_i*n_s*n_a*n_μ+i_i*n_a+n_a,t])
                end
                U_I[t]=sum(Φ[n_i*n_s*n_a*n_μ+1:n_i*n_s*n_a*n_μ+n_a,t])
            end
            println("Distribution computed for t=",t)
        end

        errors=[maximum(abs.(I-Iold));maximum(abs.(E-Eold))]
        if maximum(errors)>error
            println("Solution was diverging")
            pol_val_results=(V_E,V_U,W_E,W_U,pol_a_Ei,pol_a_Ui,pol_μ_U,pol_σ_E,pol_σ_U,J,θ,Φ)
            return Iold,Eold,U_I,U_E,pol_val_results
        end

        error=maximum(errors)

        println("error: ",error)

        Iold=I*wupdate+Iold*(1-wupdate)
        Eold=E*wupdate+Eold*(1-wupdate)

        Iold[:,T]=Iss
        Eold[:,T]=Ess

        p1=plot(1:T,Iold',legend=false)
        plot!(1:T,Eold')
        display(p1)
    end

    errorT=maximum(abs.(ΦT-Φss))

    println("error last period: ",errorT)

    pol_val_results=(V_E,V_U,W_E,W_U,pol_a_Ei,pol_a_Ui,pol_μ_U,pol_σ_E,pol_σ_U,J,θ,Φ)

    return I,E,U_I,U_E,pol_val_results
end

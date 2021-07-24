# Computing a transition

function Transition(grids,StatEq,zt;Guess=false,permanent=0,i_shock=1,n_periods=200)

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

    shockdur=size(zt,2)
    ϵ=2e-4
    wupdate=0.002

    if permanent==1
        (StatEq1,StatEq2)=StatEq
        (V_E0,V_U0,W_E0,W_U0,pol_a_Ei0,pol_a_Ui0,pol_σ_E0,pol_σ_U0,J0,θ0,Φ0,Y0,E0,U0)=StatEq
        (V_Ess,V_Uss,W_Ess,W_Uss,pol_a_Eiss,pol_a_Uiss,pol_σ_Ess,pol_σ_Uss,Jss,θss,Φss,Yss,Ess,Uss)=StatEq
    elseif permanent==0
        (V_Ess,V_Uss,W_Ess,W_Uss,pol_a_Eiss,pol_a_Uiss,pol_σ_Ess,pol_σ_Uss,Jss,θss,Φss,Yss,Ess,Uss)=StatEq
        V_E0,V_U0,W_E0,W_U0=V_Ess,V_Uss,W_Ess,W_Uss
        pol_a_Ei0,pol_a_Ui0,pol_σ_E0,pol_σ_U0=pol_a_Eiss,pol_a_Uiss,pol_σ_Ess,pol_σ_Uss
        J0,θ0,Φ0,Y0,E0,U0=Jss,θss,Φss,Yss,Ess,Uss
    end


    if Guess==false
        T=n_periods
        Eold=[zeros(n_o,n_e) for t=1:T]
        Eold[1]=E0

        if permanent==1
            step_aux=range(0,1,length=T+1)
            for t in 1:T
                Eold[t]=E0*(1-step_aux[t])+Ess*step_aux[t]
            end
        elseif permanent==0

            Eold[2][i_shock,:]=1.0*E0[i_shock,:]
            Eold[3][i_shock,:]=1.0*E0[i_shock,:]
            Eold[4][i_shock,:]=1.0*E0[i_shock,:]
            Eold[2][1:end .!=i_shock,:]=E0[1:end .!=i_shock,:]
            Eold[3][1:end .!=i_shock,:]=E0[1:end .!=i_shock,:]
            Eold[4][1:end .!=i_shock,:]=E0[1:end .!=i_shock,:]
            step_aux=range(0,1,length=min(ceil(Int64,shockdur*1.5),T))
            for t in 5:T
                if t<=length(step_aux)
                    Eold[t][i_shock,:]=Eold[4][i_shock,:]*(1-step_aux[t-1])+Ess[i_shock,:]*step_aux[t-1]
                    Eold[t][1:end .!=i_shock,:]=E0[1:end .!=i_shock,:]
                else
                    Eold[t]=E0
                end
            end
        end

    else
        Eold=Guess
        T=size(Eold)[1]
    end

    prices(ϕ,Y,z,E,e_i)=ϕ*γ[e_i]*Y^(1/ν)*z^(1-(1/ν))*E[e_i]^(ω-1)*sum(γ.*(E.^ω))^((1/ω)*(1-(1/ν))-1)

    z0=zt[:,end]
    p0=zeros(n_o,n_e)
    for o_i in 1:n_o
        for e_i in 1:n_e
            p0[o_i,e_i]=prices(ϕ[o_i],Y0,z0[o_i],E0[o_i,:],e_i)
        end
    end

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

    error=1000

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
    iter=0

    while error>ϵ
        iter+=1

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
                a_i=statestogrid_E0[ind,3]
                e_i=statestogrid_E0[ind,2]
                o_i=statestogrid_E0[ind,1]

                p_eo=p0[o_i,e_i]

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

                ub=b*p0[o_i,e_i]

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

                p_eo=p0[o_i,e_i]
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


        Φ[:,1]=Φ0

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

            # Compute the evolution of the distribution
            if t==0
                Φ[:,t+1]=Tr'*Φ0
            else
                Φ[:,t+1]=Tr'*Φ[:,t]
            end

            # Compute unemployed and employed of each type
            for o_i in 1:n_o
                for e_i in 1:n_e
                    E[t+1][o_i,e_i]=sum(Φ[(o_i-1)*n_e*n_a+(e_i-1)*n_a+1:(o_i-1)*n_e*n_a+e_i*n_a,t+1])
                    if e_i==1
                        if o_i==1
                            U[t+1][o_i,e_i]=sum(Φ[n_o*n_e*n_a+1:n_o*n_e*n_a+n_a,t+1])
                        end
                    else
                        U[t+1][o_i,e_i]=sum(Φ[n_o*n_e*n_a+n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+1:n_o*n_e*n_a+n_a+(o_i-1)*(n_e-1)*n_a+(e_i-1)*n_a,t+1])
                    end
                end
            end

            println("Distribution computed for t=",t)
        end

        errors=zeros(T)
        for t in 1:T
            errors[t]=maximum(abs.(E[t]-Eold[t]))
        end
        if maximum(errors)>error || iter>300
            println("Solution was diverging")
            pol_val_results=(V_E,V_U,W_E,W_U,pol_a_Ei,pol_a_Ui,pol_σ_E,pol_σ_U,J,θ,Φ)
            return Eold,U,pol_val_results
        end

        error=maximum(errors)

        println("iter ",iter," error: ",error)
        for t in 1:T
            Eold[t]=E[t]*wupdate+Eold[t]*(1-wupdate)
        end
        Eold[T]=Ess

        for t in 1:T
            for o_i in 1:n_o
                for e_i in 1:n_e
                    Eplot[(o_i-1)*n_e+e_i,t]=Eold[t][o_i,e_i]
                end
            end
        end
        plot1=plot(1:T,Eplot',legend=false)
        display(plot1)
    end

    pol_val_results=(V_E,V_U,W_E,W_U,pol_a_Ei,pol_a_Ui,pol_σ_E,pol_σ_U,J,θ,Φ)

    return E,U,pol_val_results
end

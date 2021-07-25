# Writing the dynamic programming with small number of grid points, but linear interpolation
# Workers' value function iteration

function ValueFunctions(grids,p;Guess=false,tol=false)

    (grid_o,grid_e,grid_a)=grids
    @eval @everywhere grid_a=$grid_a
    @eval @everywhere p=$p

    if tol==false
        ϵ=1e-6
    else
        ϵ=tol
    end

    n_o,n_e,n_a=length(grid_o),length(grid_e),length(grid_a)

    @eval @everywhere n_o=$n_o
    @eval @everywhere n_e=$n_e
    @eval @everywhere n_a=$n_a

    nstates_E=n_a*n_e*n_o
    ngrids_vars_E=[n_o,n_e,n_a]
    nstates_U=n_a+n_o*(n_e-1)*n_a
    ngrids_vars_U=[n_o,n_e,n_a]

    nsvars_E=3
    nsvars_U=3

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


    if Guess==false
        V_E_old=zeros(nstates_E)
        V_U_old=zeros(nstates_U)
        W_E_old=zeros(nstates_E)
        W_U_old=zeros(nstates_U)
        J_old=zeros(nstates_E)
        pol_σ_E_old=ones(nstates_E)
    else
        (V_E_old,V_U_old,W_E_old,W_U_old,pol_a_E_old,pol_a_U_old,pol_σ_E_old,pol_σ_U_old,J_old,θ_old)=Guess
    end

    V_E=SharedArray{Float64}(nstates_E)
    V_U=SharedArray{Float64}(nstates_U)
    W_E=SharedArray{Float64}(nstates_E)
    W_U=SharedArray{Float64}(nstates_U)

    pol_a_E=SharedArray{Float64}(nstates_E)
    pol_a_U=SharedArray{Float64}(nstates_U)

    pol_σ_E=SharedArray{Float64}(nstates_E)
    pol_σ_U=SharedArray{Float64}(nstates_U,n_o)

    J=SharedArray{Float64}(nstates_E)
    θ=SharedArray{Float64}(nstates_E)


    @everywhere u(c)=if c>0 (c^(1-σ)-1)/(1-σ) else -Inf end
    @everywhere P(θ)=min(m*θ^(1-ξ),1)
    @everywhere q_inv(y)=if y>1 0.0 elseif y<0 0.0 else (y/m)^(-1/ξ) end
    error=1000

    iter=0
    while error>ϵ && iter<3000
        iter+=1

        V_E=SharedArray{Float64}(nstates_E)
        V_U=SharedArray{Float64}(nstates_U)
        W_E=SharedArray{Float64}(nstates_E)
        W_U=SharedArray{Float64}(nstates_U)
        pol_a_E=SharedArray{Float64}(nstates_E)
        pol_a_U=SharedArray{Float64}(nstates_U)
        pol_σ_E=SharedArray{Float64}(nstates_E)
        pol_σ_U=SharedArray{Float64}(nstates_U,n_o)
        @eval @everywhere V_E_old=$V_E_old
        @eval @everywhere V_U_old=$V_U_old
        @eval @everywhere W_E_old=$W_E_old
        @eval @everywhere W_U_old=$W_U_old
        @eval @everywhere J_old=$J_old
        @eval @everywhere pol_σ_E_old=$pol_σ_E_old

        # 1) V_E
        @sync @distributed for ind in eachindex(V_E)
            a_i=statestogrid_E[ind,3]
            e_i=statestogrid_E[ind,2]
            o_i=statestogrid_E[ind,1]

            p_eo=p[o_i,e_i]

            a_guess=[grid_a[a_i]+1e-2]

            if e_i==1
                interp_V_U=LinearInterpolation(grid_a,V_U_old[1:n_a];extrapolation_bc=Line())
                ind1_en=[o_i-1,1-1,1]'*[n_e*n_a,n_a,1]
                interp_W_En=LinearInterpolation(grid_a,W_E_old[ind1_en:(ind1_en-1)+n_a];extrapolation_bc=Line())
                ind1_ep=[o_i-1,2-1,1]'*[n_e*n_a,n_a,1]
                interp_W_Ep=LinearInterpolation(grid_a,W_E_old[ind1_ep:(ind1_ep-1)+n_a];extrapolation_bc=Line())

                Veval_0(a1)= -(u((1+r)*grid_a[a_i]+φ*p_eo-a1[1])+β*(1-δ)*(ρ*interp_V_U(a1[1])+(1-ρ)*((1-α[e_i])*interp_W_En(a1[1])+α[e_i]*interp_W_Ep(a1[1]))))
                if Veval_0(a_min)<Veval_0(a_min+1e-12)
                    pol_a_E[ind]=a_min
                    V_E[ind]=-Veval_0(a_min)
                else
                    opt=optimize(Veval_0,a_guess,BFGS())
                    pol_a_E[ind]=opt.minimizer[1]
                    V_E[ind]=-opt.minimum
                end
            elseif e_i<n_e
                interp_V_Ud=LinearInterpolation(grid_a,V_U_old[1:n_a];extrapolation_bc=Line())
                interp_V_U=LinearInterpolation(grid_a,V_U_old[n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+1:n_a+(o_i-1)*(n_e-1)*n_a+(e_i-1)*n_a];extrapolation_bc=Line())
                ind1_en=[o_i-1,e_i-1,1]'*[n_e*n_a,n_a,1]
                interp_W_En=LinearInterpolation(grid_a,W_E_old[ind1_en:(ind1_en-1)+n_a];extrapolation_bc=Line())
                ind1_ep=[o_i-1,(e_i+1)-1,1]'*[n_e*n_a,n_a,1]
                interp_W_Ep=LinearInterpolation(grid_a,W_E_old[ind1_ep:(ind1_ep-1)+n_a];extrapolation_bc=Line())

                Veval_1(a1)= -(u((1+r)*grid_a[a_i]+φ*p_eo-a1[1])+β*(1-δ)*(ρ*interp_V_U(a1[1])+(1-ρ)*((1-α[e_i])*interp_W_En(a1[1])+α[e_i]*interp_W_Ep(a1[1]))))
                if Veval_1(a_min)<Veval_1(a_min+1e-12)
                    pol_a_E[ind]=a_min
                    V_E[ind]=-Veval_1(a_min)
                else
                    opt=optimize(Veval_1,a_guess,BFGS())
                    pol_a_E[ind]=opt.minimizer[1]
                    V_E[ind]=-opt.minimum
                end
            elseif e_i==n_e
                interp_V_Ud=LinearInterpolation(grid_a,V_U_old[1:n_a];extrapolation_bc=Line())
                interp_V_U=LinearInterpolation(grid_a,V_U_old[n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+1:n_a+(o_i-1)*(n_e-1)*n_a+(e_i-1)*n_a];extrapolation_bc=Line())
                ind1_e=[o_i-1,e_i-1,1]'*[n_e*n_a,n_a,1]
                interp_W_E=LinearInterpolation(grid_a,W_E_old[ind1_e:(ind1_e-1)+n_a];extrapolation_bc=Line())

                Veval_2(a1)= -(u((1+r)*grid_a[a_i]+φ*p_eo-a1[1])+β*(1-δ)*(ρ*interp_V_U(a1[1])+(1-ρ)*interp_W_E(a1[1])))
                if Veval_2(a_min)<Veval_2(a_min+1e-12)
                    pol_a_E[ind]=a_min
                    V_E[ind]=-Veval_2(a_min)
                else
                    opt=optimize(Veval_2,a_guess,BFGS())
                    pol_a_E[ind]=opt.minimizer[1]
                    V_E[ind]=-opt.minimum
                end
            end

        end

        # 2) V_U
        @sync @distributed for ind in eachindex(V_U)
            a_i=statestogrid_U[ind,3]
            e_i=statestogrid_U[ind,2]
            o_i=statestogrid_U[ind,1]

            ub=b*p[o_i,e_i]

            if e_i==1
                interp_W_U=LinearInterpolation(grid_a,W_U_old[1:n_a];extrapolation_bc=Line())

                Veval_0(a1)=-(u((1+r)*grid_a[a_i]+ub-a1[1])+β*(1-δ)*interp_W_U(a1[1]))

                a_guess=[grid_a[a_i]+1e-2]

                if Veval_0(a_min)<Veval_0(a_min+1e-12)
                    pol_a_U[ind]=a_min
                    V_U[ind]=-Veval_0(a_min)
                else
                    opt=optimize(Veval_0,a_guess,BFGS())
                    pol_a_U[ind]=opt.minimizer[1]
                    V_U[ind]=-opt.minimum
                end
            else
                interp_W_Ud=LinearInterpolation(grid_a,W_U_old[1:n_a];extrapolation_bc=Line())
                if e_i==2
                    interp_W_Ul=LinearInterpolation(grid_a,W_U_old[1:n_a];extrapolation_bc=Line())
                elseif e_i>2
                    interp_W_Ul=LinearInterpolation(grid_a,W_U_old[n_a+(o_i-1)*(n_e-1)*n_a+(e_i-3)*n_a+1:n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a];extrapolation_bc=Line())
                end
                interp_W_U=LinearInterpolation(grid_a,W_U_old[n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+1:n_a+(o_i-1)*(n_e-1)*n_a+(e_i-1)*n_a];extrapolation_bc=Line())


                Veval_1(a1)=-(u((1+r)*grid_a[a_i]+ub-a1[1])+β*(1-δ)*((1-χ[e_i])*interp_W_U(a1[1])+χ[e_i]*interp_W_Ul(a1[1])))

                a_guess=[grid_a[a_i]+1e-2]

                if Veval_1(a_min)<Veval_1(a_min+1e-12)
                    pol_a_U[ind]=a_min
                    V_U[ind]=-Veval_1(a_min)
                else
                    opt=optimize(Veval_1,a_guess,BFGS())
                    pol_a_U[ind]=opt.minimizer[1]
                    V_U[ind]=-opt.minimum
                end
            end
        end


        @sync @distributed for ind in eachindex(J)
            a_i=statestogrid_E[ind,3]
            e_i=statestogrid_E[ind,2]
            o_i=statestogrid_E[ind,1]

            p_eo=p[o_i,e_i]
            a1=pol_a_E[ind]

            if e_i<n_e
                ind1_n=[o_i-1,e_i-1,1]'*[n_e*n_a,n_a,1]
                ind1_p=[o_i-1,(e_i+1)-1,1]'*[n_e*n_a,n_a,1]
                interp_pol_σ_En=LinearInterpolation(grid_a,pol_σ_E_old[ind1_n:(ind1_n-1)+n_a];extrapolation_bc=Line())
                interp_pol_σ_Ep=LinearInterpolation(grid_a,pol_σ_E_old[ind1_p:(ind1_p-1)+n_a];extrapolation_bc=Line())
                σ_probn=interp_pol_σ_En(a1)
                σ_probp=interp_pol_σ_Ep(a1)
                interp_Jn=LinearInterpolation(grid_a,J_old[ind1_n:(ind1_n-1)+n_a];extrapolation_bc=Line())
                interp_Jp=LinearInterpolation(grid_a,J_old[ind1_p:(ind1_p-1)+n_a];extrapolation_bc=Line())
                Jn=interp_Jn(a1)
                Jp=interp_Jp(a1)
                J[ind]=(1-φ)*p_eo+β*(1-δ)*(1-ρ)*((1-α[e_i])*σ_probn*Jn+α[e_i]*σ_probp*Jp)
            elseif e_i==n_e
                ind1=[o_i-1,e_i-1,1]'*[n_e*n_a,n_a,1]
                interp_pol_σ_E=LinearInterpolation(grid_a,pol_σ_E_old[ind1:(ind1-1)+n_a];extrapolation_bc=Line())
                σ_prob=interp_pol_σ_E(a1)
                interp_J=LinearInterpolation(grid_a,J_old[ind1:(ind1-1)+n_a];extrapolation_bc=Line())
                Je=interp_J(a1)
                J[ind]=(1-φ)*p_eo+β*(1-δ)*(1-ρ)*σ_prob*Je
            end
        end

        @sync @distributed for ind in eachindex(θ)
            e_i=statestogrid_E[ind,2]
            θ[ind]=q_inv(κ/J[ind])
        end

        # Now compute Ws

        # 1) W_E
        @sync @distributed for ind in eachindex(W_E)
            a_i=statestogrid_E[ind,3]
            e_i=statestogrid_E[ind,2]
            o_i=statestogrid_E[ind,1]

            if e_i==1
                W_E[ind]=σ_ϵ*((V_E[ind]/σ_ϵ)+log(1+exp((V_U[a_i]-V_E[ind])/σ_ϵ)))
                pol_σ_E[ind]=(1+exp((V_U[a_i]-V_E[ind])/σ_ϵ))^(-1)
            else
                W_E[ind]=σ_ϵ*((V_E[ind]/σ_ϵ)+log(1+exp((V_U[n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+a_i]-V_E[ind])/σ_ϵ)))
                pol_σ_E[ind]=(1+exp((V_U[n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+a_i]-V_E[ind])/σ_ϵ))^(-1)
            end
        end

        # 2) W_U
        V_job_s=SharedArray{Float64}(length(W_U),n_o)
        @sync @distributed for ind in eachindex(W_U)
            a_i=statestogrid_U[ind,3]
            e_i=statestogrid_U[ind,2]
            o_i=statestogrid_U[ind,1]

            if e_i==1
                for o1_i in 1:n_o
                    ste=[o1_i-1,1-1,a_i]'*[n_e*n_a,n_a,1]
                    prob=P(θ[ste])
                    V_job_s[ind,o1_i]=prob*V_E[ste]+(1-prob)*V_U[ind]
                end

                W_U[ind]=σ_ϵ*((V_job_s[ind,1]/σ_ϵ)+log(1+(O-1)*exp((V_job_s[ind,2]-V_job_s[ind,1])/σ_ϵ)))
                for o1_i in 1:n_o
                    if o1_i==1
                        pol_σ_U[ind,o1_i]=(1+(O-1)*exp((V_job_s[ind,2]-V_job_s[ind,o1_i])/σ_ϵ))^(-1)
                    else
                        pol_σ_U[ind,o1_i]=(1+(1/(O-1))*exp((V_job_s[ind,1]-V_job_s[ind,o1_i])/σ_ϵ))^(-1)
                    end
                end
            else
                if o_i==1
                    for o1_i in 1:n_o
                        if o1_i==o_i
                            e1_i=e_i
                            stu=ind
                        else
                            e1_i=1
                            stu=ind
                        end
                        ste=[o1_i-1,e1_i-1,a_i]'*[n_e*n_a,n_a,1]

                        prob=P(θ[ste])
                        V_job_s[ind,o1_i]=prob*V_E[ste]+(1-prob)*V_U[stu]
                    end

                    W_U[ind]=σ_ϵ*((V_job_s[ind,1]/σ_ϵ)+log(1+(O-1)*exp((V_job_s[ind,2]-V_job_s[ind,1])/σ_ϵ)))
                    pol_σ_U[ind,1]=(1+(O-1)*exp((V_job_s[ind,2]-V_job_s[ind,1])/σ_ϵ))^(-1)
                    pol_σ_U[ind,2]=(1+(1/(O-1))*exp((V_job_s[ind,1]-V_job_s[ind,2])/σ_ϵ))^(-1)
                elseif o_i==2
                    for o1_i in 1:n_o
                        if o1_i==o_i
                            e1_i=e_i
                            stu=ind
                        else
                            e1_i=1
                            stu=ind
                        end
                        ste=[o1_i-1,e1_i-1,a_i]'*[n_e*n_a,n_a,1]

                        prob=P(θ[ste])
                        V_job_s[ind,o1_i]=prob*V_E[ste]+(1-prob)*V_U[stu]
                    end

                    ste=[2-1,1-1,a_i]'*[n_e*n_a,n_a,1]
                    stu=ind

                    prob=P(θ[ste])
                    V_job2inexperienced=prob*V_E[ste]+(1-prob)*V_U[stu]

                    W_U[ind]=σ_ϵ*((V_job_s[ind,1]/σ_ϵ)+log(1+exp((V_job_s[ind,2]-V_job_s[ind,1])/σ_ϵ)+(O-2)*exp((V_job2inexperienced-V_job_s[ind,1])/σ_ϵ)))
                    pol_σ_U[ind,1]=(1+exp((V_job_s[ind,2]-V_job_s[ind,1])/σ_ϵ)+(O-2)*exp((V_job2inexperienced-V_job_s[ind,1])/σ_ϵ))^(-1)
                    pol_σ_U[ind,2]=(1+exp((V_job_s[ind,1]-V_job_s[ind,2])/σ_ϵ)+(O-2)*exp((V_job2inexperienced-V_job_s[ind,2])/σ_ϵ))^(-1)
                end
            end
        end

        error=maximum(abs.(W_E-W_E_old))+maximum(abs.(W_U-W_U_old))+maximum(abs.(V_E-V_E_old))+maximum(abs.(V_U-V_U_old))+maximum(abs.(J-J_old))
        println("iter ",iter,": error=",error)
        W_E,W_E_old=W_E_old,W_E
        W_U,W_U_old=W_U_old,W_U
        V_E,V_E_old=V_E_old,V_E
        V_U,V_U_old=V_U_old,V_U
        J_old=J_old*0.75+J*0.25
        pol_σ_E,pol_σ_E_old=pol_σ_E_old,pol_σ_E
    end

    W_E,W_E_old=W_E_old,W_E
    W_U,W_U_old=W_U_old,W_U
    V_E,V_E_old=V_E_old,V_E
    V_U,V_U_old=V_U_old,V_U
    J,J_old=J_old,J
    pol_σ_E,pol_σ_E_old=pol_σ_E_old,pol_σ_E

    return V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ
end


@everywhere function transformVPolFunctions(pol_val_functions,grids0,grid_a1)

    (grid_o,grid_e,grid_a)=grids0
    n_o,n_e,n_a=length(grid_o),length(grid_e),length(grid_a)
    n_anew=length(grid_a1)


    (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions


    nstates_E=n_anew*n_e*n_o
    ngrids_vars_E=[n_o,n_e,n_anew]
    nstates_U=n_anew+n_o*(n_e-1)*n_anew
    ngrids_vars_U=[n_o,n_e,n_anew]

    nsvars_E=3
    nsvars_U=3

    statestogrid_E=zeros(Int64,nstates_E,nsvars_E)
    for v in 1:nsvars_E
        statestogrid_E[:,v]=kron(ones(prod(ngrids_vars_E[1:v-1]),1),kron(1:ngrids_vars_E[v],ones(prod(ngrids_vars_E[v+1:nsvars_E]),1)))
    end

    statestogrid_U=zeros(Int64,nstates_U,nsvars_U)
    statestogrid_U[1:n_anew,:]=hcat(ones(n_anew,2),1:n_anew)
    for o_i in 1:n_o
        statestogrid_U[n_anew+(o_i-1)*(n_e-1)*n_anew+1:n_anew+o_i*(n_e-1)*n_anew,:]=hcat(o_i*ones((n_e-1)*n_anew,1),kron(2:n_e,ones(n_anew)),kron(ones(n_e-1),1:n_anew))
    end

    V_Enew=zeros(nstates_E)
    V_Unew=zeros(nstates_U)
    W_Enew=zeros(nstates_E)
    W_Unew=zeros(nstates_U)

    pol_a_Enew=zeros(nstates_E)
    pol_a_Unew=zeros(nstates_U)
    pol_σ_Enew=zeros(nstates_E)
    pol_σ_Unew=zeros(nstates_U,n_o)

    Jnew=zeros(nstates_E)
    θnew=zeros(nstates_E)

    for ind in eachindex(V_Enew)
        a_i=statestogrid_E[ind,3]
        e_i=statestogrid_E[ind,2]
        o_i=statestogrid_E[ind,1]

        ind1=[o_i-1,e_i-1,1]'*[n_e*n_a,n_a,1]

        interp_V_E=LinearInterpolation(grid_a,V_E[ind1:(ind1-1)+n_a])
        V_Enew[ind]=interp_V_E(grid_a1[a_i])
        interp_pol_a_E=LinearInterpolation(grid_a,pol_a_E[ind1:(ind1-1)+n_a])
        pol_a_Enew[ind]=interp_pol_a_E(grid_a1[a_i])
        interp_W_E=LinearInterpolation(grid_a,W_E[ind1:(ind1-1)+n_a])
        W_Enew[ind]=interp_W_E(grid_a1[a_i])
        interp_pol_σ_E=LinearInterpolation(grid_a,pol_σ_E[ind1:(ind1-1)+n_a])
        pol_σ_Enew[ind]=interp_pol_σ_E(grid_a1[a_i])

        interp_J=LinearInterpolation(grid_a,J[ind1:(ind1-1)+n_a])
        Jnew[ind]=interp_J(grid_a1[a_i])
        interp_θ=LinearInterpolation(grid_a,θ[ind1:(ind1-1)+n_a])
        θnew[ind]=interp_θ(grid_a1[a_i])
    end

    for ind in eachindex(V_Unew)
        a_i=statestogrid_U[ind,3]
        e_i=statestogrid_U[ind,2]
        o_i=statestogrid_U[ind,1]

        if e_i==1
            ind1=1
        else
            ind1=n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+1
        end

        interp_V_U=LinearInterpolation(grid_a,V_U[ind1:(ind1-1)+n_a])
        V_Unew[ind]=interp_V_U(grid_a1[a_i])
        interp_pol_a_U=LinearInterpolation(grid_a,pol_a_U[ind1:(ind1-1)+n_a])
        pol_a_Unew[ind]=interp_pol_a_U(grid_a1[a_i])
        interp_W_U=LinearInterpolation(grid_a,W_U[ind1:(ind1-1)+n_a])
        W_Unew[ind]=interp_W_U(grid_a1[a_i])
        for o_i in 1:n_o
            interp_pol_σ_U=LinearInterpolation(grid_a,pol_σ_U[ind1:(ind1-1)+n_a,o_i])
            pol_σ_Unew[ind,o_i]=interp_pol_σ_U(grid_a1[a_i])
        end
    end
    pol_val_functions_new=(V_Enew,V_Unew,W_Enew,W_Unew,pol_a_Enew,pol_a_Unew,pol_σ_Enew,pol_σ_Unew,Jnew,θnew)
    return pol_val_functions_new
end


function multigrid(nGrids_a,p)

    pol_val_functions=false
    pol_val_functione_int=false

    for j in 1:length(nGrids_a)
        grid_a=LinRange(a_min,a_max,nGrids_a[j])
        grids=(grid_o,grid_e,grid_a)

        V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ=ValueFunctions(grids,p; Guess=pol_val_functione_int)
        pol_val_functions=(V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ)

        if j<length(nGrids_a)
            n_anew=nGrids_a[j+1]
            grid_anew=LinRange(a_min,a_max,n_anew)
            pol_val_functione_int=transformVPolFunctions(pol_val_functions,grids,grid_anew)
        end
        println("Problem solved for a grid of size ",nGrids_a[j],". ",length(nGrids_a)-j," more steps for final solution")
    end
    return pol_val_functions
end

@everywhere function transformPola(pol_a_E,pol_a_U,grids)

    (grid_o,grid_e,grid_a)=grids
    n_o,n_e,n_a=length(grid_o),length(grid_e),length(grid_a)

    nstates_E=n_a*n_e*n_o
    ngrids_vars_E=[n_o,n_e,n_a]
    nstates_U=n_a+n_o*(n_e-1)*n_a
    ngrids_vars_U=[n_o,n_e,n_a]

    nsvars_E=3
    nsvars_U=3

    statestogrid_E=zeros(Int64,nstates_E,nsvars_E)
    for v in 1:nsvars_E
        statestogrid_E[:,v]=kron(ones(prod(ngrids_vars_E[1:v-1]),1),kron(1:ngrids_vars_E[v],ones(prod(ngrids_vars_E[v+1:nsvars_E]),1)))
    end

    statestogrid_U=zeros(Int64,nstates_U,nsvars_U)
    statestogrid_U[1:n_a,:]=hcat(ones(n_a,2),1:n_a)
    for o_i in 1:n_o
        statestogrid_U[n_a+(o_i-1)*(n_e-1)*n_a+1:n_a+o_i*(n_e-1)*n_a,:]=hcat(o_i*ones((n_e-1)*n_a,1),kron(2:n_e,ones(n_a)),kron(ones(n_e-1),1:n_a))
    end

    pol_a_Ei=zeros(Int64,nstates_E)
    pol_a_Ui=zeros(Int64,nstates_U)

    for ind in eachindex(pol_a_Ei)
        a_i=statestogrid_E[ind,3]
        e_i=statestogrid_E[ind,2]
        o_i=statestogrid_E[ind,1]

        ind1=[o_i-1,e_i-1,1]'*[n_e*n_a,n_a,1]

        interp_aux=LinearInterpolation(grid_a,1:n_a)
        if pol_a_E[ind]<=a_max && pol_a_E[ind]>=a_min
            pol_a_Ei[ind]=ceil(interp_aux(pol_a_E[ind]))
        elseif pol_a_E[ind]>a_max
            pol_a_Ei[ind]=n_a
        elseif pol_a_E[ind]<a_min
            pol_a_Ei[ind]=1
        end

    end

    for ind in eachindex(pol_a_Ui)
        a_i=statestogrid_U[ind,3]
        e_i=statestogrid_U[ind,2]
        o_i=statestogrid_U[ind,1]

        if e_i==1
            ind1=1
        else
            ind1=n_a+(o_i-1)*(n_e-1)*n_a+1
        end

        interp_aux=LinearInterpolation(grid_a,1:n_a)
        if pol_a_U[ind]<=a_max && pol_a_U[ind]>=a_min
            pol_a_Ui[ind]=ceil(interp_aux(pol_a_U[ind]))
        elseif pol_a_U[ind]>a_max
            pol_a_Ui[ind]=n_a
        elseif pol_a_U[ind]<a_min
            pol_a_Ui[ind]=1
        end

    end
    return pol_a_Ei,pol_a_Ui
end

# Writing the dynamic programming with small number of grid points, but linear interpolation
# Workers' value function iteration

function VFunctionIterEq(grids,w,θ;Vguess=false,tol=false)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    @eval @everywhere grid_a=$grid_a
    @eval @everywhere grid_μ=$grid_μ
    @eval @everywhere w=$w
    @eval @everywhere θ=$θ

    if tol==false
        ϵ=1e-5
    else
        ϵ=tol
    end

    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)

    @eval @everywhere n_i=$n_i
    @eval @everywhere n_s=$n_s
    @eval @everywhere n_a=$n_a
    @eval @everywhere n_μ=$n_μ

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

    @eval @everywhere statestogrid_E=$statestogrid_E
    @eval @everywhere statestogrid_U=$statestogrid_U

    if Vguess==false
        V_E_old=zeros(nstates_E)
        V_U_old=zeros(nstates_U)
        W_E_old=zeros(nstates_E)
        W_U_old=zeros(nstates_U)
    else
        V_E_old,V_U_old,W_E_old,W_U_old=Vguess
    end

    V_E=SharedArray{Float64}(nstates_E)
    V_U=SharedArray{Float64}(nstates_U)
    W_E=SharedArray{Float64}(nstates_E)
    W_U=SharedArray{Float64}(nstates_U)

    pol_a_E=SharedArray{Float64}(nstates_E)
    pol_a_U=SharedArray{Float64}(nstates_U)

    pol_μ_U=SharedArray{Int64}(nstates_U,n_i)

    pol_σ_E=SharedArray{Float64}(nstates_E)
    pol_σ_U=SharedArray{Float64}(nstates_U,n_i)


    @everywhere u(c)=if c>0 (c^(1-σ)-1)/(1-σ) else -Inf end
    @everywhere p(θ)=min(m*θ^(1-ξ),1)
    error=1000

    iter=0
    while error>ϵ && iter<1000
        iter+=1

        V_E=SharedArray{Float64}(nstates_E)
        V_U=SharedArray{Float64}(nstates_U)
        W_E=SharedArray{Float64}(nstates_E)
        W_U=SharedArray{Float64}(nstates_U)
        pol_a_E=SharedArray{Float64}(nstates_E)
        pol_a_U=SharedArray{Float64}(nstates_U)
        pol_μ_U=SharedArray{Int64}(nstates_U,n_i)
        pol_σ_E=SharedArray{Float64}(nstates_E)
        pol_σ_U=SharedArray{Float64}(nstates_U,n_i)
        @eval @everywhere V_E_old=$V_E_old
        @eval @everywhere V_U_old=$V_U_old
        @eval @everywhere W_E_old=$W_E_old
        @eval @everywhere W_U_old=$W_U_old

        # 1) V_E
        @sync @distributed for ind in eachindex(V_E)
            μ_i=statestogrid_E[ind,4]
            a_i=statestogrid_E[ind,3]
            s_i=statestogrid_E[ind,2]
            i_i=statestogrid_E[ind,1]

            wage=w[i_i,s_i]

            a_guess=[grid_a[a_i]+1e-2]

            if s_i==1
                interp_V_U=LinearInterpolation(grid_a,V_U_old[1:n_a];extrapolation_bc=Line())
                ind1_ei=[i_i-1,1-1,1-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                interp_W_Ei=LinearInterpolation(grid_a,W_E_old[ind1_ei:n_μ:(ind1_ei-μ_i)+n_a*n_μ];extrapolation_bc=Line())
                ind1_ee=[i_i-1,2-1,1-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                interp_W_Ee=LinearInterpolation(grid_a,W_E_old[ind1_ee:n_μ:(ind1_ee-μ_i)+n_a*n_μ];extrapolation_bc=Line())

                Veval_i(a1)= -(u((1+r)*grid_a[a_i]+grid_μ[μ_i]*wage-a1[1])+β*(ρ*interp_V_U(a1[1])+(1-ρ)*((1-α)*interp_W_Ei(a1[1])+α*interp_W_Ee(a1[1]))))
                if Veval_i(a_min)<Veval_i(a_min+1e-12)
                    pol_a_E[ind]=a_min
                    V_E[ind]=-Veval_i(a_min)
                else
                    opt=optimize(Veval_i,a_guess,BFGS())
                    pol_a_E[ind]=opt.minimizer[1]
                    V_E[ind]=-opt.minimum
                end
            elseif s_i==2
                interp_V_Ui=LinearInterpolation(grid_a,V_U_old[1:n_a];extrapolation_bc=Line())
                interp_V_Ue=LinearInterpolation(grid_a,V_U_old[i_i*n_a+1:(i_i+1)*n_a];extrapolation_bc=Line())
                ind1_e=[i_i-1,2-1,1-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                interp_W_E=LinearInterpolation(grid_a,W_E_old[ind1_e:n_μ:(ind1_e-μ_i)+n_a*n_μ];extrapolation_bc=Line())

                Veval_e(a1)= -(u((1+r)*grid_a[a_i]+grid_μ[μ_i]*wage-a1[1])+β*((ρ-δ)*interp_V_Ue(a1[1])+δ*interp_V_Ui(a1[1])+(1-ρ)*interp_W_E(a1[1])))
                if Veval_e(a_min)<Veval_e(a_min+1e-12)
                    pol_a_E[ind]=a_min
                    V_E[ind]=-Veval_e(a_min)
                else
                    opt=optimize(Veval_e,a_guess,BFGS())
                    pol_a_E[ind]=opt.minimizer[1]
                    V_E[ind]=-opt.minimum
                end
            end

        end

        # 2) V_U
        @sync @distributed for ind in eachindex(V_U)
            a_i=statestogrid_U[ind,3]
            s_i=statestogrid_U[ind,2]
            i_i=statestogrid_U[ind,1]

            if s_i==1
                interp_W_U=LinearInterpolation(grid_a,W_U_old[1:n_a];extrapolation_bc=Line())

                Veval_i(a1)=-(u((1+r)*grid_a[a_i]+b-a1[1])+β*interp_W_U(a1[1]))

                a_guess=[grid_a[a_i]+1e-2]

                if Veval_i(a_min)<Veval_i(a_min+1e-12)
                    pol_a_U[ind]=a_min
                    V_U[ind]=-Veval_i(a_min)
                else
                    opt=optimize(Veval_i,a_guess,BFGS())
                    pol_a_U[ind]=opt.minimizer[1]
                    V_U[ind]=-opt.minimum
                end
            elseif s_i==2
                interp_W_U=LinearInterpolation(grid_a,W_U_old[i_i*n_a+1:(i_i+1)*n_a];extrapolation_bc=Line())
                interp_W_Ui=LinearInterpolation(grid_a,W_U_old[1:n_a];extrapolation_bc=Line())

                Veval_e(a1)=-(u((1+r)*grid_a[a_i]+b-a1[1])+β*((1-δ)*interp_W_U(a1[1])+δ*interp_W_Ui(a1[1])))

                a_guess=[grid_a[a_i]+1e-2]

                if Veval_e(a_min)<Veval_e(a_min+1e-12)
                    pol_a_U[ind]=a_min
                    V_U[ind]=-Veval_e(a_min)
                else
                    opt=optimize(Veval_e,a_guess,BFGS())
                    pol_a_U[ind]=opt.minimizer[1]
                    V_U[ind]=-opt.minimum
                end
            end
        end


        # Now compute Ws

        # 1) W_E
        @sync @distributed for ind in eachindex(W_E)
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
        V_job_s=SharedArray{Float64}(length(W_U),n_i)
        @sync @distributed for ind in eachindex(W_U)
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
                            stu=i1_i*n_a+a_i
                        else
                            s1_i=1
                            stu=a_i
                        end
                        ste=[i1_i-1,s1_i-1,a_i-1,μ1_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]

                        prob=p(θ[ste])
                        Val_temp[μ1_i]=prob*V_E[ste]+(1-prob)*V_U[stu]
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
        V_E,V_E_old=V_E_old,V_E
        V_U,V_U_old=V_U_old,V_U
    end

    W_E,W_E_old=W_E_old,W_E
    W_U,W_U_old=W_U_old,W_U
    V_E,V_E_old=V_E_old,V_E
    V_U,V_U_old=V_U_old,V_U

    return V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U
end

# Firms' value function iteration

function JFunctionIter(grids,w,policyfunctions_W; Jguess=false,tol=false)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    @eval @everywhere grid_a=$grid_a
    @eval @everywhere grid_μ=$grid_μ

    (pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U)=policyfunctions_W

    @eval @everywhere pol_a_E=$pol_a_E

    if tol==false
        ϵ=1e-6
    else
        ϵ=tol
    end

    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)

    @eval @everywhere n_i=$n_i
    @eval @everywhere n_s=$n_s
    @eval @everywhere n_a=$n_a
    @eval @everywhere n_μ=$n_μ

    nstates=n_μ*n_a*n_s*n_i
    nsvars=4
    ngrids_vars=[n_i,n_s,n_a,n_μ]

    statestogrid=zeros(Int64,nstates,nsvars)
    for v in 1:nsvars
        statestogrid[:,v]=kron(ones(prod(ngrids_vars[1:v-1]),1),kron(1:ngrids_vars[v],ones(prod(ngrids_vars[v+1:nsvars]),1)))
    end

    @eval @everywhere statestogrid=$statestogrid

    if Jguess==false
        J_old=zeros(nstates)
    else
        J_old=zeros(nstates)
        J_old[:]=Jguess
    end

    J=SharedArray{Float64}(nstates)

    error=1000
    iter=0
    while error>ϵ
        iter+=1

        J=SharedArray{Float64}(nstates)
        @eval @everywhere J_old=$J_old

        @sync @distributed for ind in eachindex(J)
            μ_i=statestogrid[ind,4]
            a_i=statestogrid[ind,3]
            s_i=statestogrid[ind,2]
            i_i=statestogrid[ind,1]

            wage=w[i_i,s_i]
            a1=pol_a_E[ind]

            if s_i==1
                ind1_i=[i_i-1,s_i-1,1-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_e=[i_i-1,2-1,1-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                interp_pol_σ_Ei=LinearInterpolation(grid_a,pol_σ_E[ind1_i:n_μ:(ind1_i-μ_i)+n_a*n_μ];extrapolation_bc=Line())
                interp_pol_σ_Ee=LinearInterpolation(grid_a,pol_σ_E[ind1_e:n_μ:(ind1_e-μ_i)+n_a*n_μ];extrapolation_bc=Line())
                σ_probi=interp_pol_σ_Ei(a1)
                σ_probe=interp_pol_σ_Ee(a1)
                interp_Ji=LinearInterpolation(grid_a,J_old[ind1_i:n_μ:(ind1_i-μ_i)+n_a*n_μ];extrapolation_bc=Line())
                interp_Je=LinearInterpolation(grid_a,J_old[ind1_e:n_μ:(ind1_e-μ_i)+n_a*n_μ];extrapolation_bc=Line())
                Ji=interp_Ji(a1)
                Je=interp_Je(a1)
                J[ind]=(1-grid_μ[μ_i])*wage+β*(1-ρ)*((1-α)*σ_probi*Ji+α*σ_probe*Je)
            elseif s_i==2
                ind1=[i_i-1,s_i-1,1-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                interp_pol_σ_E=LinearInterpolation(grid_a,pol_σ_E[ind1:n_μ:(ind1-μ_i)+n_a*n_μ];extrapolation_bc=Line())
                σ_prob=interp_pol_σ_E(a1)
                interp_J=LinearInterpolation(grid_a,J_old[ind1:n_μ:(ind1-μ_i)+n_a*n_μ];extrapolation_bc=Line())
                Jaux=interp_J(a1)
                J[ind]=(1-grid_μ[μ_i])*wage+β*(1-ρ)*σ_prob*Jaux
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

    q_inv(y)=if y>1 0.0 elseif y<0 0.0 else (y/m)^(-1/ξ) end

    θ=zeros(n_μ*n_a*n_s*n_i)
    if Guess==false
        J_old=2.0*ones(n_μ*n_a*n_s*n_i)
        Vfunctions=false
        policyfunctions=false
    else
        (V_E_old,V_U_old,W_E_old,W_U_old,pol_a_E_old,pol_a_U_old,pol_μ_U_old,pol_σ_E_old,pol_σ_U_old,J_old,θ_old)=Guess
        Vfunctions=(V_E_old,V_U_old,W_E_old,W_U_old)
        policyfunctions=(pol_a_E_old,pol_a_U_old,pol_μ_U_old,pol_σ_E_old,pol_σ_U_old)
    end

    dampening=0.5


    ϵ=1e-6
    error=1000

    iter=0
    while error>1e-6 && iter<50
        iter+=1
        for ind in eachindex(θ)
            θ[ind]=q_inv(κ/(J_old[ind]-F))
        end

        V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U=VFunctionIterEq(grids,w,θ,Vguess=Vfunctions,tol=1e-6)

        policyfunctions=(pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U)
        Vfunctions=(V_E,V_U,W_E,W_U)

        J=JFunctionIter(grids,w,policyfunctions,Jguess=J_old,tol=1e-9)

        error=maximum((J-J_old).^2)
        println("iter ",iter," in outward loop, error of ",error)
        if dampening==0.0
            J,J_old=J_old,J
        else
            J_old=dampening*J_old+(1-dampening)*J
        end
    end
    J,J_old=J_old,J

    for ind in eachindex(θ)
        θ[ind]=q_inv(κ/(J[ind]-F))
    end
    (V_E,V_U,W_E,W_U)=Vfunctions
    (pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U)=policyfunctions

    return V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ
end


@everywhere function transformVPolFunctions(pol_val_functions,grids0,grid_a1)

    (grid_i,grid_s,grid_a,grid_μ)=grids0
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    n_anew=length(grid_a1)


    (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions



    nstates_E=n_μ*n_anew*n_s*n_i
    ngrids_vars_E=[n_i,n_s,n_anew,n_μ]
    nstates_U=n_anew+n_i*n_anew
    ngrids_vars_U=[n_i,n_s,n_anew]

    nsvars_E=4
    nsvars_U=3

    statestogrid_E=zeros(Int64,nstates_E,nsvars_E)
    for v in 1:nsvars_E
        statestogrid_E[:,v]=kron(ones(prod(ngrids_vars_E[1:v-1]),1),kron(1:ngrids_vars_E[v],ones(prod(ngrids_vars_E[v+1:nsvars_E]),1)))
    end

    statestogrid_U=zeros(Int64,nstates_U,nsvars_U)
    statestogrid_U[1:n_anew,:]=hcat(ones(n_anew,2),1:n_anew)
    for i_i in 1:n_i
        statestogrid_U[n_anew+(i_i-1)*n_anew+1:n_anew+i_i*n_anew,:]=hcat(i_i*ones(n_anew,1),2*ones(n_anew,1),1:n_anew)
    end

    V_Enew=zeros(nstates_E)
    V_Unew=zeros(nstates_U)
    W_Enew=zeros(nstates_E)
    W_Unew=zeros(nstates_U)

    pol_a_Enew=zeros(nstates_E)
    pol_a_Unew=zeros(nstates_U)
    pol_μ_Unew=zeros(Int64,nstates_U,n_i)
    pol_σ_Enew=zeros(nstates_E)
    pol_σ_Unew=zeros(nstates_U,n_i)

    Jnew=zeros(nstates_E)
    θnew=zeros(nstates_E)

    for ind in eachindex(V_Enew)
        μ_i=statestogrid_E[ind,4]
        a_i=statestogrid_E[ind,3]
        s_i=statestogrid_E[ind,2]
        i_i=statestogrid_E[ind,1]

        ind1=[i_i-1,s_i-1,1-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_a,1]

        interp_V_E=LinearInterpolation(grid_a,V_E[ind1:n_μ:(ind1-μ_i)+n_a*n_μ])
        V_Enew[ind]=interp_V_E(grid_a1[a_i])
        interp_pol_a_E=LinearInterpolation(grid_a,pol_a_E[ind1:n_μ:(ind1-μ_i)+n_a*n_μ])
        pol_a_Enew[ind]=interp_pol_a_E(grid_a1[a_i])
        interp_W_E=LinearInterpolation(grid_a,W_E[ind1:n_μ:(ind1-μ_i)+n_a*n_μ])
        W_Enew[ind]=interp_W_E(grid_a1[a_i])
        interp_pol_σ_E=LinearInterpolation(grid_a,pol_σ_E[ind1:n_μ:(ind1-μ_i)+n_a*n_μ])
        pol_σ_Enew[ind]=interp_pol_σ_E(grid_a1[a_i])

        interp_J=LinearInterpolation(grid_a,J[ind1:n_μ:(ind1-μ_i)+n_a*n_μ])
        Jnew[ind]=interp_J(grid_a1[a_i])
        interp_θ=LinearInterpolation(grid_a,θ[ind1:n_μ:(ind1-μ_i)+n_a*n_μ])
        θnew[ind]=interp_θ(grid_a1[a_i])
    end

    for ind in eachindex(V_Unew)
        a_i=statestogrid_U[ind,3]
        s_i=statestogrid_U[ind,2]
        i_i=statestogrid_U[ind,1]

        if s_i==1
            ind1=1
        elseif s_i==2
            ind1=i_i*n_a+1
        end

        interp_V_U=LinearInterpolation(grid_a,V_U[ind1:(ind1-1)+n_a])
        V_Unew[ind]=interp_V_U(grid_a1[a_i])
        interp_pol_a_U=LinearInterpolation(grid_a,pol_a_U[ind1:(ind1-1)+n_a])
        pol_a_Unew[ind]=interp_pol_a_U(grid_a1[a_i])
        interp_W_U=LinearInterpolation(grid_a,W_U[ind1:(ind1-1)+n_a])
        W_Unew[ind]=interp_W_U(grid_a1[a_i])
        for i_i in 1:n_i
            interp_pol_σ_U=LinearInterpolation(grid_a,pol_σ_U[ind1:(ind1-1)+n_a,i_i])
            pol_σ_Unew[ind,i_i]=interp_pol_σ_U(grid_a1[a_i])
            interp_pol_μ_U=LinearInterpolation(grid_a,pol_μ_U[ind1:(ind1-1)+n_a,i_i])
            pol_μ_Unew[ind,i_i]=ceil(interp_pol_μ_U(grid_a1[a_i]))
        end
    end
    pol_val_functions_new=(V_Enew,V_Unew,W_Enew,W_Unew,pol_a_Enew,pol_a_Unew,pol_μ_Unew,pol_σ_Enew,pol_σ_Unew,Jnew,θnew)
    return pol_val_functions_new
end


function multigrid(nGrids_a,w)

    pol_val_functions=false
    pol_val_functions_int=false

    for j in 1:length(nGrids_a)
        grid_a=LinRange(a_min,a_max,nGrids_a[j])
        grids=(grid_i,grid_s,grid_a,grid_μ)

        V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ=ValueFunctions(grids,w; Guess=pol_val_functions_int)
        pol_val_functions=(V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)

        if j<length(nGrids_a)
            n_anew=nGrids_a[j+1]
            grid_anew=LinRange(a_min,a_max,n_anew)
            pol_val_functions_int=transformVPolFunctions(pol_val_functions,grids,grid_anew)
        end
        println("Problem solved for a grid of size ",nGrids_a[j],". ",length(nGrids_a)-j," more steps for final solution")
    end
    return pol_val_functions
end

@everywhere function transformPola(pol_a_E,pol_a_U,grids)

    (grid_i,grid_s,grid_a,grid_μ)=grids
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

    pol_a_Ei=zeros(Int64,nstates_E)
    pol_a_Ui=zeros(Int64,nstates_U)

    for ind in eachindex(pol_a_Ei)
        μ_i=statestogrid_E[ind,4]
        a_i=statestogrid_E[ind,3]
        s_i=statestogrid_E[ind,2]
        i_i=statestogrid_E[ind,1]

        ind1=[i_i-1,s_i-1,1-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_a,1]

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
        s_i=statestogrid_U[ind,2]
        i_i=statestogrid_U[ind,1]

        if s_i==1
            ind1=1
        elseif s_i==2
            ind1=i_i*n_a+1
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

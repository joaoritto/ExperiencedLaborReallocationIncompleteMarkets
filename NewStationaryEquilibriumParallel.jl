
# Computing the stationary distribution

function ComputeDistribution(grids,pol_functions)

    p(θ)=min(m*θ^(1-ξ),1)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    (pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,θ)=pol_functions

    nstates=n_i*n_s*n_a*n_μ+n_a+n_i*(n_s-1)*n_a
    nsvars=5
    ngrids_vars=[2,n_i,n_s,n_a,n_μ]

    nstates_E=n_i*n_s*n_a*n_μ
    statestogrid_E=ones(Int64,nstates_E,nsvars)
    nstates_U=n_a+n_i*(n_s-1)*n_a
    statestogrid_U=ones(Int64,nstates_U,nsvars)
    for v in 1:nsvars
        if v==1
            statestogrid_E[:,v]=1*ones(nstates_E,1)
        else
            statestogrid_E[:,v]=kron(ones(prod(ngrids_vars[2:v-1]),1),kron(1:ngrids_vars[v],ones(prod(ngrids_vars[v+1:nsvars]),1)))
        end
    end
    statestogrid_U[:,1]=2*ones(nstates_U)
    statestogrid_U[1:n_a,2:nsvars-1]=hcat(ones(n_a,2),1:n_a)
    for i_i in 1:n_i
        statestogrid_U[n_a+(i_i-1)*(n_s-1)*n_a+1:n_a+i_i*(n_s-1)*n_a,2:nsvars-1]=hcat(i_i*ones((n_s-1)*n_a,1),kron(2:n_s,ones(n_a)),kron(ones(n_s-1),1:n_a))
    end

    statestogrid=[statestogrid_E;statestogrid_U]

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
            a1_i=pol_a_E[indE]
            if s_i==1
                ind1_en=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_ep=[i_i-1,(s_i+1)-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_un=n_i*n_s*n_a*n_μ+a1_i
                ind1_up=n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-1)*n_a+a1_i

                push!(i,ind)
                push!(j,ind1_en)
                push!(k,(1-ρ)*(1-α[s_i])*pol_σ_E[ind1_en])

                push!(i,ind)
                push!(j,ind1_un)
                push!(k,(1-ρ)*(1-α[s_i])*(1-pol_σ_E[ind1_en])+ρ)

                push!(i,ind)
                push!(j,ind1_ep)
                push!(k,(1-ρ)*α[s_i]*pol_σ_E[ind1_ep])

                push!(i,ind)
                push!(j,ind1_up)
                push!(k,(1-ρ)*α[s_i]*(1-pol_σ_E[ind1_ep]))

            elseif s_i<n_s
                ind1_en=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_ep=[i_i-1,(s_i+1)-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_un=n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+a1_i
                ind1_up=n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-1)*n_a+a1_i
                ind1_ud=n_i*n_s*n_a*n_μ+a1_i

                push!(i,ind)
                push!(j,ind1_en)
                push!(k,(1-ρ)*(1-α[s_i])*pol_σ_E[ind1_en])

                push!(i,ind)
                push!(j,ind1_un)
                push!(k,(1-ρ)*(1-α[s_i])*(1-pol_σ_E[ind1_en])+ρ-δ)

                push!(i,ind)
                push!(j,ind1_ep)
                push!(k,(1-ρ)*α[s_i]*pol_σ_E[ind1_ep])

                push!(i,ind)
                push!(j,ind1_up)
                push!(k,(1-ρ)*α[s_i]*(1-pol_σ_E[ind1_ep]))

                push!(i,ind)
                push!(j,ind1_ud)
                push!(k,δ)

            elseif s_i==n_s
                ind1_e=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_un=n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+a1_i
                ind1_ud=n_i*n_s*n_a*n_μ+a1_i

                push!(i,ind)
                push!(j,ind1_e)
                push!(k,(1-ρ)*pol_σ_E[ind1_e])

                push!(i,ind)
                push!(j,ind1_un)
                push!(k,(1-ρ)*(1-pol_σ_E[ind1_e])+ρ-δ)

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

            a1_i=pol_a_U[indU]

            if s_i==1

                for i1_i in 1:n_i
                    ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                    ind1_e[i1_i]=[i1_i-1,s_i-1,a1_i-1,pol_μ_U[a1_i,i1_i]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]

                    push!(i,ind)
                    push!(j,ind1_e[i1_i])
                    push!(k,pol_σ_U[a1_i,i1_i]*p(θ[ind1_e[i1_i]]))

                    push!(i,ind)
                    push!(j,ind1_u[i1_i])
                    push!(k,pol_σ_U[a1_i,i1_i]*(1-p(θ[ind1_e[i1_i]])))
                end

            else

                for i1_i in 1:n_i
                    if i1_i==i_i
                        ind1_u[i1_i]=ind
                        ind1_e[i1_i]=[i1_i-1,s_i-1,a1_i-1,pol_μ_U[n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+a1_i,i1_i]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    else
                        ind1_u[i1_i]=ind
                        ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    end

                    push!(i,ind)
                    push!(j,ind1_e[i1_i])
                    push!(k,(1-δ-χ[s_i])*pol_σ_U[n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+a1_i,i1_i]*p(θ[ind1_e[i1_i]]))

                    push!(i,ind)
                    push!(j,ind1_u[i1_i])
                    push!(k,(1-δ-χ[s_i])*pol_σ_U[n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+a1_i,i1_i]*(1-p(θ[ind1_e[i1_i]])))

                    if s_i==2

                        ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                        ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]

                        push!(i,ind)
                        push!(j,ind1_e[i1_i])
                        push!(k,(δ+χ[s_i])*pol_σ_U[a1_i,i1_i]*p(θ[ind1_e[i1_i]]))

                        push!(i,ind)
                        push!(j,ind1_u[i1_i])
                        push!(k,(δ+χ[s_i])*pol_σ_U[a1_i,i1_i]*(1-p(θ[ind1_e[i1_i]])))

                    else
                        ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                        ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]

                        push!(i,ind)
                        push!(j,ind1_e[i1_i])
                        push!(k,δ*pol_σ_U[a1_i,i1_i]*p(θ[ind1_e[i1_i]]))

                        push!(i,ind)
                        push!(j,ind1_u[i1_i])
                        push!(k,δ*pol_σ_U[a1_i,i1_i]*(1-p(θ[ind1_e[i1_i]])))

                        if i1_i==i_i
                            ind1_u[i1_i]=n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-3)*n_a+a1_i
                            ind1_e[i1_i]=[i1_i-1,(s_i-1)-1,a1_i-1,pol_μ_U[n_a+(i_i-1)*(n_s-1)*n_a+(s_i-3)*n_a+a1_i,i1_i]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        else
                            ind1_u[i1_i]=n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-3)*n_a+a1_i
                            ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        end

                        push!(i,ind)
                        push!(j,ind1_e[i1_i])
                        push!(k,χ[s_i]*pol_σ_U[n_a+(i_i-1)*(n_s-1)*n_a+(s_i-3)*n_a+a1_i,i1_i]*p(θ[ind1_e[i1_i]]))

                        push!(i,ind)
                        push!(j,ind1_u[i1_i])
                        push!(k,χ[s_i]*pol_σ_U[n_a+(i_i-1)*(n_s-1)*n_a+(s_i-3)*n_a+a1_i,i1_i]*(1-p(θ[ind1_e[i1_i]])))
                    end

                end
            end
        end
    end

    push!(i,nstates)
    push!(j,nstates)
    push!(k,0.0)

    T=sparse(i,j,k)

    Φ=zeros(nstates)
    ind0=200
    Φ[ind0]=1.0

    for j in 1:7000
        Φ=T'*Φ
    end

    return Φ,T
end

# Computing aggregates in the stationary equilibrium

function ComputeAggregates(grids,pol_functions,Φ,z)

    p(θ)=min(m*θ^(1-ξ),1)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)

    nstates=n_i*n_s*n_a*n_μ+n_i*n_a

    # Compute unemployed and employed of each type
    E=zeros(n_i,n_s)
    U=zeros(n_i,n_s)
    for i_i in 1:n_i
        for s_i in 1:n_s
            E[i_i,s_i]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+(s_i-1)*n_a*n_μ+1:(i_i-1)*n_s*n_a*n_μ+s_i*n_μ*n_a])
            if s_i==1
                if i_i==1
                    U[i_i,s_i]=sum(Φ[n_i*n_s*n_a*n_μ+1:n_i*n_s*n_a*n_μ+n_a])
                end
            else
                U[i_i,s_i]=sum(Φ[n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-2)*n_a+1:n_i*n_s*n_a*n_μ+n_a+(i_i-1)*(n_s-1)*n_a+(s_i-1)*n_a])
            end
        end
    end

    Y=0
    y=zeros(n_i)
    for i_i in 1:n_i
        y[i_i]=z[i_i]*prod(E[i_i,:].^γ)
        Y+=(1/n_i)*y[i_i]^((ν-1)/ν)
    end
    Y=Y^(ν/(ν-1))

    return Y,E,U
end

function GeneralEquilibrium(z)

    n_anew=1500 #nGrids_a[end]

    if maximum(z[2:end]-z[1:end-1])!=0
        diffz=true
        println("message: different productivities")
        return
    else
        diffz=false
    end

    ϵ=1e-6

    function Y_CES(z,E)
        Y=0
        for i_i in 1:n_i
            Y+=(1/n_i)*(z[i_i]*prod(E[i_i,:].^γ))^((ν-1)/ν)
        end
        Y=Y^(ν/(ν-1))
        return Y
    end

    pol_val_functions=false

    El=[0.029737932381096246 - 5e-6;
        0.05017787352590754 - 5e-6;
        0.07629272108164575 - 5e-6;
        0.3137706336240057 - 5e-6]
    Eu=[0.029737932381096246 + 5e-6;
         0.05017787352590754 + 5e-6;
        0.07629272108164575 + 5e-6;
        0.3137706336240057 + 5e-6]

    Es=zeros(n_s)
    U=zeros(n_s)
    Y=false
    Φ=false
    Tr=false
    Eerr=1000*ones(n_s)

    Ed=zeros(n_i,n_s)
    k=1
    Eiter=0
    ready=false

    while ready==false
        Eiter+=1
        Ed[:,k]=ones(n_i,1)*(El[k]*0.5+Eu[k]*0.5)

        if k<n_s
            k=k+1
            continue
        end

        Y=Y_CES(z,Ed)

        wages(Y,z,E,s_i)=(1/n_i)*γ[s_i]*Y^(1/ν)*z^(1-(1/ν))*E[s_i]^(γ[s_i]*(1-(1/ν))-1)*prod(E[1:end .!=s_i].^((γ[1:end .!=s_i])*(1-(1/ν))))

        w=zeros(n_i,n_s)
        for i_i in 1:n_i
            for s_i in 1:n_s
                w[i_i,s_i]=wages(Y,z[i_i],Ed[i_i,:],s_i)
            end
        end

        display(w)
        println("E=",Ed[1,:])

        pol_val_functions=multigrid(nGrids_a,w)
        (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions

        grid_a=LinRange(a_min,a_max,nGrids_a[end])
        grids=(grid_i,grid_s,grid_a,grid_μ)
        if n_anew!=nGrids_a[end]
            grid_a=LinRange(a_min,a_max,n_anew)
            pol_val_functions=transformVPolFunctions(pol_val_functions,grids,grid_a)
            (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions
            grids=(grid_i,grid_s,grid_a,grid_μ)
        end

        pol_a_Ei,pol_a_Ui=transformPola(pol_a_E,pol_a_U,grids)

        pol_functions=(pol_a_Ei,pol_a_Ui,pol_μ_U,pol_σ_E,pol_σ_U,θ)
        pol_val_functions=(V_E,V_U,W_E,W_U,pol_a_Ei,pol_a_Ui,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)

        Φ,Tr=ComputeDistribution(grids,pol_functions)
        Y,Es,U=ComputeAggregates(grids,pol_functions,Φ,z)
        display(Es)

        while k>0
            Eerr[k]=Es[1,k]-Ed[1,k]
            println("E error for s_i=",k," is ",Eerr[k])
            if abs(Eerr[k])>ϵ && Eiter<12
                if Eerr[k]>0
                    El[k]=Ed[1,k]
                else
                    Eu[k]=Ed[1,k]
                end
                break
            else
                k=k-1
                Eiter=0
            end
        end
        if k==0
            ready=true
        end
    end

    return pol_val_functions,Φ,Tr,Y,Es,U
end

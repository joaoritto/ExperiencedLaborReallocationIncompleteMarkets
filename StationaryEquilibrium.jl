
# Computing the stationary distribution

function ComputeDistribution(para,grids,pol_val_functions)

    p(θ)=min(m*θ^(1-ξ),1)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    (V_E,V_U,V_S,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions

    nstates=n_i*n_s*n_a*n_μ+n_i*n_s*n_a
    nsvars=5
    ngrids_vars=[3,n_i,n_s,n_a,n_μ]

    nstates_E=n_i*n_s*n_a*n_μ
    statestogrid_E=ones(Int64,nstates_E,nsvars)
    nstates_U=n_i*n_s*n_a
    statestogrid_U=ones(Int64,nstates_U,nsvars)
    for v in 1:nsvars
        if v==1
            statestogrid_E[:,v]=1*ones(nstates_E,1)
            statestogrid_U[:,v]=2*ones(nstates_U,1)
        else
            statestogrid_E[:,v]=kron(ones(prod(ngrids_vars[2:v-1]),1),kron(1:ngrids_vars[v],ones(prod(ngrids_vars[v+1:nsvars]),1)))
            if v==nsvars
            else
                statestogrid_U[:,v]=kron(ones(prod(ngrids_vars[2:v-1]),1),kron(1:ngrids_vars[v],ones(prod(ngrids_vars[v+1:nsvars-1]),1)))
            end
        end
    end
    statestogrid=[statestogrid_E;statestogrid_U]

    # Construct Transition matrix
    T1=spzeros(nstates,nstates)
    T2=spzeros(nstates,nstates)

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
                ind1_e0=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_e1=[i_i-1,2-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_u=n_i*n_s*n_a*n_μ+[i_i-1,s_i-1,a1_i]'*[n_s*n_a,n_a,1]


                T1[ind,ind1_e0]=(1-ρ)*(1-α)
                T1[ind,ind1_e1]=(1-ρ)*α
                T1[ind,ind1_u]=ρ-δ

                ind1_s=zeros(Int64,n_i)
                for i1_i in 1:n_i
                    ind1_s[i1_i]=n_i*n_s*n_a*n_μ+[i1_i-1,1-1,a1_i]'*[n_s*n_a,n_a,1]
                end

                for i1_i in 1:n_i
                    T1[ind,ind1_s[i1_i]]+=δ*(1/n_i)
                end

            elseif s_i==2
                ind1_e=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_u=n_i*n_s*n_a*n_μ+[i_i-1,s_i-1,a1_i]'*[n_s*n_a,n_a,1]
                ind1_s=n_i*n_s*n_a*n_μ+n_i*n_s*n_a+a1_i

                T1[ind,ind1_e]=(1-ρ)
                T1[ind,ind1_u]=ρ-δ

                ind1_s=zeros(Int64,n_i)
                for i1_i in 1:n_i
                    ind1_s[i1_i]=n_i*n_s*n_a*n_μ+[i1_i-1,1-1,a1_i]'*[n_s*n_a,n_a,1]
                end

                for i1_i in 1:n_i
                    T1[ind,ind1_s[i1_i]]+=δ*(1/n_i)
                end
            end


            ind1_s_u=zeros(Int64,n_i)
            ind1_s_e=zeros(Int64,n_i)
            for i1_i in 1:n_i
                ind1_s_u[i1_i]=n_i*n_s*n_a*n_μ+[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]
                ind1_s_e[i1_i]=[i1_i-1,1-1,a_i-1,pol_μ_U[[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
            end

            T2[ind,ind]=pol_σ_E[indE]
            for i1_i in 1:n_i
                T2[ind,ind1_s_u[i1_i]]+=(1-pol_σ_E[indE])*(1/n_i)*(1-p(θ[[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]]))
                T2[ind,ind1_s_e[i1_i]]+=(1-pol_σ_E[indE])*(1/n_i)*p(θ[[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]])
            end

        elseif e_i==2

            indU=[i_i-1,s_i-1,a_i]'*[n_s*n_a,n_a,1]
            a1_i=pol_a_U[indU]

            ind1_u=n_i*n_s*n_a*n_μ+[i_i-1,s_i-1,a1_i]'*[n_s*n_a,n_a,1]

            T1[ind,ind1_u]=1

            ind1_s_u=zeros(Int64,n_i)
            ind1_s_e=zeros(Int64,n_i)
            for i1_i in 1:n_i
                ind1_s_u[i1_i]=n_i*n_s*n_a*n_μ+[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]
                ind1_s_e[i1_i]=[i1_i-1,1-1,a_i-1,pol_μ_U[[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
            end
            ind1_e=[i_i-1,s_i-1,a_i-1,pol_μ_U[indU]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]


            T2[ind,ind1_u]=pol_σ_U[indU]*(1-p(θ[indU]))
            T2[ind,ind1_e]=pol_σ_U[indU]*p(θ[indU])
            for i1_i in 1:n_i
                T2[ind,ind1_s_u[i1_i]]+=(1-pol_σ_U[indU])*(1/n_i)*(1-p(θ[[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]]))
                T2[ind,ind1_s_e[i1_i]]+=(1-pol_σ_U[indU])*(1/n_i)*p(θ[[i1_i-1,1-1,a_i]'*[n_s*n_a,n_a,1]])
            end

        end
    end
    T=spzeros(nstates,nstates)
    T=T1'*T2'

    Φ=zeros(nstates)
    ind0=200
    Φ[ind0]=1.0

    for i in 1:300
        Φ=T*Φ
    end

    return Φ
end

# Computing aggregates in the stationary equilibrium

function ComputeAggregates(para,grids,pol_val_functions,Φ,z)

    p(θ)=min(m*θ^(1-ξ),1)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    (V_E,V_U,V_S,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions

    nstates=n_i*n_s*n_a*n_μ+n_i*n_s*n_a+n_a
    nsvars=5
    ngrids_vars=[3,n_i,n_s,n_a,n_μ]

    # Compute unemployed and employed of each type
    E_I=zeros(n_i)
    E_E=zeros(n_i)
    U_I=zeros(n_i)
    U_E=zeros(n_i)
    for i_i in 1:n_i
        E_I[i_i]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+1:(i_i-1)*n_s*n_a*n_μ+n_μ*n_a])
        E_E[i_i]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+n_μ*n_a+1:(i_i-1)*n_s*n_a*n_μ+2*n_μ*n_a])
        U_I[i_i]=sum(Φ[n_i*n_s*n_a*n_μ+(i_i-1)*n_s*n_a+1:n_i*n_s*n_a*n_μ+(i_i-1)*n_s*n_a+n_a])
        U_E[i_i]=sum(Φ[n_i*n_s*n_a*n_μ+(i_i-1)*n_s*n_a+n_a+1:n_i*n_s*n_a*n_μ+(i_i-1)*n_s*n_a+2*n_a])
    end

    Y=0
    y=zeros(n_i)
    for i_i in 1:n_i
        y[i_i]=z[i_i]*E_I[i_i]^γ*E_E[i_i]^(1-γ)
        Y+=(1/n_i)*y[i_i]^((ν-1)/ν)
    end
    Y=Y^(ν/(ν-1))

    return Y,E_I,E_E,U_I,U_E
end

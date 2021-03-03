
# Plotting the results at the stationary equilibrium

function PlotResults(para,grids,pol_val_functions,Φ)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    (V_E,V_U,V_S,W_E,W_U,W_S,pol_a_E,pol_a_U,pol_a_S,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions

    nstates=n_i*n_s*n_a*n_μ+n_i*n_s*n_a+n_a
    nsvars=5
    ngrids_vars=[3,n_i,n_s,n_a,n_μ]

    statestogrid=zeros(nstates)

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

    #plot(pol_μ_U[1:n_a])
    #for state in 2:n_s*n_i
    #    plot!(pol_μ_U[(state-1)*n_a+1:state*n_a])
    #end

    plot(pol_μ_U[1:n_a])
    plot!(pol_μ_U[n_a+1:2*n_a])
    plot!(pol_μ_U[2*n_a+1:3*n_a])
    plot!(pol_μ_U[3*n_a+1:4*n_a])

    μ_a=zeros(n_a)
    for a_i in 1:n_a
        μ_a[a_i]=Φ[(a_i-1)*n_μ+1:a_i*n_μ]'*grid_μ[statestogrid[(a_i-1)*n_μ+1:a_i*n_μ,5]]/sum(Φ[(a_i-1)*n_μ+1:a_i*n_μ])
    end

    plot(μ_a)

    μ_a2=zeros(n_a)
    for a_i in 1:n_a
        μ_a2[a_i]=Φ[n_a*n_μ+(a_i-1)*n_μ+1:n_a*n_μ+a_i*n_μ]'*grid_μ[statestogrid[n_a*n_μ+(a_i-1)*n_μ+1:n_a*n_μ+a_i*n_μ,5]]/sum(Φ[n_a*n_μ+(a_i-1)*n_μ+1:n_a*n_μ+a_i*n_μ])
    end

    plot!(μ_a2)
    yaxis!([0,1])

    d_a=zeros(n_a)
    for a_i in 1:n_a
        d_a[a_i]=sum(Φ[(a_i-1)*n_μ+1:a_i*n_μ])
    end
    plot(d_a)

end


# Plotting the results at the stationary equilibrium

function PlotResults(grids,pol_val_functions,Φ)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions

    nstates=n_i*n_s*n_a*n_μ+n_a+n_i*n_a
    nsvars=5
    ngrids_vars=[2,n_i,n_s,n_a,n_μ]

    nstates_E=n_i*n_s*n_a*n_μ
    statestogrid_E=ones(Int64,nstates_E,nsvars)
    nstates_U=n_a+n_i*n_a
    statestogrid_U=ones(Int64,nstates_U,nsvars)
    for v in 1:nsvars
        if v==1
            statestogrid_E[:,v]=1*ones(nstates_E,1)
        else
            statestogrid_E[:,v]=kron(ones(prod(ngrids_vars[2:v-1]),1),kron(1:ngrids_vars[v],ones(prod(ngrids_vars[v+1:nsvars]),1)))
        end
    end
    statestogrid_U[:,1]=2*ones(nstates_U,1)
    statestogrid_U[:,2]=vcat(ones(n_a,1),kron(1:n_i,ones(n_a,1)))
    statestogrid_U[:,4]=kron(ones(n_i+1,1),1:n_a)

    statestogrid=[statestogrid_E;statestogrid_U]

    #---------------------------------------------------------------------------
    # 1) Compute the policy functions of submarket choice of the unemployed
    #---------------------------------------------------------------------------

    myplot=plot(pol_μ_U[1:n_a])
    for state in 2:1+n_i
        myplot=plot!(pol_μ_U[(state-1)*n_a+1:state*n_a])
    end

    #---------------------------------------------------------------------------
    # 2) Compute the distribution of people of different wealth by submarket
    #---------------------------------------------------------------------------

    μ_dist=zeros(n_a,n_i*n_s)
    for i_i in 1:n_i
        for s_i in 1:n_s
            for a_i in 1:n_a
                ind1=[i_i-1,s_i-1,a_i-1,1]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                indend=[i_i-1,s_i-1,a_i-1,n_μ]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                μ_dist[a_i,(i_i-1)*n_s+s_i]=Φ[ind1:indend]'*grid_μ[statestogrid[ind1:indend,5]]/sum(Φ[ind1:indend])
            end
        end
    end


    myplot2=plot(μ_dist[:,1])
    for state in 2:n_s*n_i
        myplot2=plot!(μ_dist[:,state])
    end
    yaxis!([0,1])

    # For each wealth level what is the mu distribution?
    dist_μa=zeros(n_μ,n_a,n_i,n_s)
    Dist_μa=zeros(n_μ,n_a,n_i,n_s)
    for a_i in 1:n_a
        for i_i in 1:n_i
            for s_i in 1:n_s
                ind0=[i_i-1,s_i-1,a_i-1,1]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                for μ_i in 1:n_μ
                    ind=[i_i-1,s_i-1,a_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    Dist_μa[μ_i,a_i,i_i,s_i]=sum(Φ[ind0:ind])
                end
            end
        end
    end

    #---------------------------------------------------------------------------
    # 3) Compute assortment in stationary equilibrium
    #---------------------------------------------------------------------------

    # a) mean wealth of workers in each sector
    a_bar=zeros(n_i,n_s)
    aux=zeros(Int64,n_μ*n_a,1)
    for i_i in 1:n_i
        for s_i in 1:n_s
            ind1=[i_i-1,s_i-1,1-1,1]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
            indend=[i_i-1,s_i-1,n_a-1,n_μ]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
            aux[:]=kron(1:n_a,ones(n_μ,1))
            a_bar[i_i,s_i]=sum(Φ[ind1:indend].*grid_a[aux])/sum(Φ[ind1:indend])
        end
    end

    # b) Percentage of people with given wealth on each sector/skill level
    prob_is_a=zeros(n_a,n_i*n_s)
    for a_i in 1:n_a
        Total=0
        for i_i in 1:n_i
            for s_i in 1:n_s
                ind1=[i_i-1,s_i-1,a_i-1,1]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                indend=[i_i-1,s_i-1,a_i-1,n_μ]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                prob_is_a[a_i,(i_i-1)*n_s+s_i]=sum(Φ[ind1:indend])
                Total+=sum(Φ[ind1:indend])
            end
        end
       prob_is_a[a_i,:]=prob_is_a[a_i,:]/Total
    end



    #---------------------------------------------------------------------------
    # 4) Construct a wealth histogram
    #---------------------------------------------------------------------------

    dist_a=zeros(n_a)
    Dist_a=zeros(n_a)
    for a_i in 1:n_a
        for i_i in 1:n_i
            for s_i in 1:n_s
                ind1=[i_i-1,s_i-1,a_i-1,1]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                indend=[i_i-1,s_i-1,a_i-1,n_μ]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                dist_a[a_i]+=sum(Φ[ind1:indend])
            end
        end

        for s_i in 1:n_s
            if s_i==1
                ind=n_i*n_s*n_a*n_μ+a_i
                dist_a[a_i]+=sum(Φ[ind])
            elseif s_i==2
                for i_i in 1:n_i
                    ind=n_i*n_s*n_a*n_μ+i_i*n_a+a_i
                    dist_a[a_i]+=sum(Φ[ind])
                end
            end
        end
        Dist_a[a_i]=sum(dist_a[1:a_i])
    end

    #---------------------------------------------------------------------------
    # 5) Simulations
    #---------------------------------------------------------------------------
    # i) Worker with little income unemployed


    # ii) Worker with little wealth who is experienced in the better sector


    # iii) Worker with little wealth who is experienced in the worse sector

    return myplot,myplot2,a_bar,prob_is_a
end

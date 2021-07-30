

#function PlotResultsStatEq(grids,StatEq)

    (grid_beq,grid_o,grid_e,grid_a)=grids
    n_beq,n_o,n_e,n_a=length(grid_beq),length(grid_o),length(grid_e),length(grid_a)
    (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ,Φ,Y,E,U)=StatEq

    nstates=n_beq*(n_o*n_e*n_a+n_a+n_o*(n_e-1)*n_a)
    nsvars=5
    ngrids_vars=[2,n_beq,n_o,n_e,n_a]

    nstates_E=n_a*n_e*n_o*n_beq
    ngrids_vars_E=[n_beq,n_o,n_e,n_a]
    nstates_U=n_beq*(n_a+n_o*(n_e-1)*n_a)
    ngrids_vars_U=[n_beq,n_o,n_e-1,n_a]

    nsvars_E=4
    nsvars_U=4

    statestogrid_E=zeros(Int64,nstates_E,nsvars_E)
    for v in 1:nsvars_E
        statestogrid_E[:,v]=kron(ones(prod(ngrids_vars_E[1:v-1]),1),kron(1:ngrids_vars_E[v],ones(prod(ngrids_vars_E[v+1:nsvars_E]),1)))
    end

    statestogrid_U=zeros(Int64,nstates_U,nsvars_U)
    for beq_i in 1:n_beq
        statestogrid_U[(beq_i-1)*n_a+1:beq_i*n_a,:]=hcat(beq_i*ones(n_a),ones(n_a,2),1:n_a)
    end
    for v in 1:nsvars_E
        if v==3
            statestogrid_U[n_beq*n_a+1:end,v]=kron(ones(prod(ngrids_vars_U[1:v-1]),1),kron(2:ngrids_vars_U[v]+1,ones(prod(ngrids_vars_U[v+1:nsvars_U]),1)))
        else
            statestogrid_U[n_beq*n_a+1:end,v]=kron(ones(prod(ngrids_vars_U[1:v-1]),1),kron(1:ngrids_vars_U[v],ones(prod(ngrids_vars_U[v+1:nsvars_U]),1)))
        end
    end

    statestogrid=[hcat(ones(Int64,size(statestogrid_E,1),1),statestogrid_E);hcat(2*ones(Int64,size(statestogrid_U,1),1),statestogrid_U)]


    #---------------------------------------------------------------------------
    #  Mean wealth of workers in stationary equilibrium
    #---------------------------------------------------------------------------
#=
    a_bar=zeros(n_o,n_e)
    aux=zeros(Int64,n_beq*n_a,1)
    for o_i in 1:n_o
        for e_i in 1:n_e
            A=zeros(Bool,nstates)
            for beq_i in 1:n_beq
                ind1=[beq_i-1,o_i-1,e_i-1,1]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
                A[ind1:ind1-1+n_a].=1
            end
            a_bar[o_i,e_i]=sum(Φ[A].*kron(ones(n_beq),grid_a))/sum(Φ[A])
        end
    end
    display(a_bar)
=#

    #---------------------------------------------------------------------------
    #  Construct a wealth histogram in stationary equilibrium
    #---------------------------------------------------------------------------

    dist_a=zeros(n_a)
    Dist_a=zeros(n_a)
    for a_i in 1:n_a
        for beq_i in 1:n_beq
            for o_i in 1:n_o
                for e_i in 1:n_e
                    ind=[beq_i-1,o_i-1,e_i-1,a_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
                    dist_a[a_i]+=Φ[ind]
                end
            end
        end

        for beq_i in 1:n_beq
            for e_i in 1:n_e
                if e_i==1
                    ind=n_beq*n_o*n_e*n_a+(beq_i-1)*n_a+a_i
                    dist_a[a_i]+=Φ[ind]
                else
                    for o_i in 1:n_o
                        ind=n_beq*n_o*n_e*n_a+n_beq*n_a+[beq_i-1,o_i-1,(e_i-1)-1,a_i]'*[n_o*(n_e-1)*n_a,(n_e-1)*n_a,n_a,1]
                        dist_a[a_i]+=Φ[ind]
                    end
                end
            end
        end
        Dist_a[a_i]=sum(dist_a[1:a_i])
    end

    plot2=plot(grid_a/(φ*p[1,1]),Dist_a,title="Distribution of wealth",titlefontsize=7)
    display(plot2)

#end

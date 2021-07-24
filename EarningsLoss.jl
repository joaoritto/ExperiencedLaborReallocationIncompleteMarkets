# Computing the earnings loss from job loss

function EarningsLoss(grids,pol_functions,Tr,p)

    (grid_i,grid_s,grid_a)=grids
    n_o,n_e,n_a=length(grid_i),length(grid_s),length(grid_a)
    (pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,θ)=pol_functions

    nstates=n_o*n_e*n_a+n_a+n_o*(n_e-1)*n_a
    nsvars=4
    ngrids_vars=[2,n_o,n_e,n_a]

    nstates_E=n_o*n_e*n_a
    statestogrid_E=ones(Int64,nstates_E,nsvars)
    nstates_U=n_a+n_o*(n_e-1)*n_a
    statestogrid_U=ones(Int64,nstates_U,nsvars)
    for v in 1:nsvars
        if v==1
            statestogrid_E[:,v]=1*ones(nstates_E,1)
        else
            statestogrid_E[:,v]=kron(ones(prod(ngrids_vars[2:v-1]),1),kron(1:ngrids_vars[v],ones(prod(ngrids_vars[v+1:nsvars]),1)))
        end
    end
    statestogrid_U[:,1]=2*ones(nstates_U)
    statestogrid_U[1:n_a,2:nsvars]=hcat(ones(n_a,2),1:n_a)
    for o_i in 1:n_o
        statestogrid_U[n_a+(o_i-1)*(n_e-1)*n_a+1:n_a+o_i*(n_e-1)*n_a,2:nsvars]=hcat(o_i*ones((n_e-1)*n_a,1),kron(2:n_e,ones(n_a)),kron(ones(n_e-1),1:n_a))
    end

    statestogrid=[statestogrid_E;statestogrid_U]

    T=30
    earningsdist=zeros(nstates)
    for ind in 1:nstates
        a_i=statestogrid[ind,4]
        e_i=statestogrid[ind,3]
        o_i=statestogrid[ind,2]
        s_i=statestogrid[ind,1]


        if s_i==1
            earningsdist[ind]=p[o_i,e_i]
        elseif s_i==2
            earningsdist[ind]=p[o_i,e_i]*b
        end
    end

    earningsloss0=zeros(nstates_E,T)
    earningsloss1=zeros(nstates_E,T)
    for ind in eachindex(V_E)
        a_i=statestogrid_E[ind,4]
        e_i=statestogrid_E[ind,3]
        o_i=statestogrid_E[ind,2]

        a1_i=pol_a_E[ind]
        indE=ind
        if e_i==1
            indU=n_o*n_e*n_a+a1_i
        else
            indU=n_o*n_e*n_a+n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+a1_i
        end

        distE=zeros(nstates)
        distE[indE]=1.0
        distU=zeros(nstates)
        distU[indU]=1.0
        for t in 1:T
            distE=Tr'*distE
            distU=Tr'*distU
            earningsE0=distE'*earningsdist
            earningsU0=distU'*earningsdist
            earningsE1=(distE[1:nstates_E]'*earningsdist[1:nstates_E])/sum(distE[1:nstates_E])
            earningsU1=(distU[1:nstates_E]'*earningsdist[1:nstates_E])/sum(distU[1:nstates_E])
            earningsloss0[ind,t]=1-earningsU0/earningsE0
            earningsloss1[ind,t]=1-earningsU1/earningsE1
        end
    end
    return earningsloss0,earningsloss1
end

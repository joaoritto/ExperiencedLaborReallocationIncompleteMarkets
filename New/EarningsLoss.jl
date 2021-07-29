# Computing the earnings loss from job loss

function EarningsLoss(grids,pol_functions,Tr,p)

    (grid_beq,grid_o,grid_e,grid_a)=grids
    n_beq,n_o,n_e,n_a=length(grid_beq),length(grid_o),length(grid_e),length(grid_a)
    (pol_a_Ei,pol_a_Ui,pol_σ_E,pol_σ_U,θ)=pol_functions

    @eval @everywhere pol_a_Ei=$pol_a_Ei
    @eval @everywhere Tr=$Tr
    @eval @everywhere p=$p

    @eval @everywhere grid_beq=$grid_beq
    @eval @everywhere grid_o=$grid_o
    @eval @everywhere grid_e=$grid_e
    @eval @everywhere grid_a=$grid_a


    @eval @everywhere n_beq=$n_beq
    @eval @everywhere n_o=$n_o
    @eval @everywhere n_e=$n_e
    @eval @everywhere n_a=$n_a

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

    @eval @everywhere statestogrid_E=$statestogrid_E
    @eval @everywhere statestogrid_U=$statestogrid_U
    @eval @everywhere statestogrid=$statestogrid

    T=30
    earningsdist=SharedArray{Float64}(nstates)
    @sync @distributed for ind in 1:nstates
        a_i=statestogrid[ind,5]
        e_i=statestogrid[ind,4]
        o_i=statestogrid[ind,3]
        beq_i=statestogrid[ind,2]
        s_i=statestogrid[ind,1]


        if s_i==1
            earningsdist[ind]=φ*p[o_i,e_i]
        elseif s_i==2
            earningsdist[ind]=p[o_i,e_i]*b
        end
    end

    earningsloss0=SharedArray{Float64}(nstates_E,T)
    earningsloss1=SharedArray{Float64}(nstates_E,T)
    @sync @distributed for ind in 1:1500:nstates_E
        a_i=statestogrid_E[ind,4]
        e_i=statestogrid_E[ind,3]
        o_i=statestogrid_E[ind,2]
        beq_i=statestogrid_E[ind,1]

        println(ind)

        a1_i=pol_a_Ei[ind]
        indE=ind
        if e_i==1
            indU=n_beq*n_o*n_e*n_a+(beq_i-1)*n_a+a1_i
        else
            indU=n_beq*n_o*n_e*n_a+n_beq*n_a+[beq_i-1,o_i-1,(e_i-1)-1,a1_i]'*[n_o*(n_e-1)*n_a,(n_e-1)*n_a,n_a,1]
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

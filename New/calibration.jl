function tenure_returns(E)

    γnew=zeros(4)

    R=[1; 1.089; 1.20; 1.279]

    Lratio=[1; E[1,2]/E[1,1]; E[1,3]/E[1,1]; E[1,4]/E[1,1]].^(1-ω)

    γnew[1]=(Lratio'*R)^(-1)
    γnew=(R.*Lratio)*γnew[1]
    return γnew
end

function quintile_w_to_y(grids,E,Dist_a,p)

    (grid_beq,grid_o,grid_e,grid_a)=grids

    quintiles_w=zeros(5)
    quintiles_y=zeros(5)

    Ecumul=zeros(n_e)
    for e_i in 1:n_e
        Ecumul[e_i]=sum(E[:,1:e_i])
    end

    aux=zeros(Int64,5)
    for quintile in 1:5
        aux[quintile]=findmin(abs.(Dist_a.-0.2*quintile))[2]
        if quintile==1
            quintiles_w[quintile]=sum(grid_a[1:aux[quintile]].*Dist_a[1:aux[quintile]])/sum(Dist_a[1:aux[quintile]])
        else
            quintiles_w[quintile]=sum(grid_a[aux[quintile-1]+1:aux[quintile]].*Dist_a[aux[quintile-1]+1:aux[quintile]])/sum(Dist_a[aux[quintile-1]+1:aux[quintile]])
        end
    end

    quintiles_y[1]=([sum(E[:,1]);sum(E[:,2]);0.2-Ecumul[2]]'*p[1,1:3]*φ)/0.2
    quintiles_y[2]=([sum(E[:,3])-(0.2-Ecumul[2]);0.4-(sum(E[:,3])-(0.2-Ecumul[2]))]'*p[1,3:4]*φ)/0.2
    quintiles_y[3]=p[1,4]*φ
    quintiles_y[4]=p[1,4]*φ
    quintiles_y[5]=p[1,4]*φ

    w_to_y=quintiles_w./quintiles_y
    #(sum(p.*E)*φ)

    return w_to_y
end


function NewDistribution(grids,Φ)

    (grid_beq,grid_o,grid_e,grid_a)=grids
    nstates=length(Φ)

    n_beq,n_o,n_e,n_a=length(grid_beq),length(grid_o),length(grid_e),length(grid_a)

    @eval @everywhere n_beq=$n_beq
    @eval @everywhere n_o=$n_o
    @eval @everywhere n_e=$n_e
    @eval @everywhere n_a=$n_a

    nstates_E=n_beq*n_o*n_e*n_a
    nstates_U=n_beq*(n_a+n_o*(n_e-1)*n_a)
    nsvars_E=4
    ngrids_vars_E=[n_beq,n_o,n_e,n_a]
    nsvars_U=4
    ngrids_vars_U=[n_beq,n_o,n_e-1,n_a]

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

    @eval @everywhere statestogrid=$statestogrid
    @eval @everywhere Φ=$Φ

    Φnew=SharedArray{Float64}(nstates)

    A=Φ'*kron(ones(div(nstates,n_a)),grid_a)

    A_i=findmin(abs.(grid_a.-A))[2]

    @sync @distributed for ind in 1:nstates
        a_i=statestogrid[ind,5]
        e_i=statestogrid[ind,4]
        o_i=statestogrid[ind,3]
        beq_i=statestogrid[ind,2]
        s_i=statestogrid[ind,1]

        if s_i==1
            indnew=[beq_i-1,o_i-1,e_i-1,A_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
        elseif s_i==2
            if e_i==1
                indnew=n_beq*n_o*n_e*n_a+(beq_i-1)*n_a+A_i
            else
                indnew=n_beq*n_o*n_e*n_a+n_beq*n_a+[beq_i-1,o_i-1,(e_i-1)-1,A_i]'*[n_o*(n_e-1)*n_a,(n_e-1)*n_a,n_a,1]
            end
        end

        Φnew[indnew]+=Φ[ind]
    end
    return Φnew
end

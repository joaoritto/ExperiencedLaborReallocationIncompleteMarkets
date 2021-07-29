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

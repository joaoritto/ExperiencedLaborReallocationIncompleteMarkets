
# Computing the stationary distribution

function ComputeDistribution(grids,pol_functions)

    p(θ)=min(m*θ^(1-ξ),1)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    (pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,θ)=pol_functions

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

    # Construct Transition matrix
    T=spzeros(nstates,nstates)

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
                ind1_u0=n_i*n_s*n_a*n_μ+a1_i
                ind1_u1=n_i*n_s*n_a*n_μ+i_i*n_a+a1_i


                T[ind,ind1_e0]=(1-ρ)*(1-α)*pol_σ_E[ind1_e0]
                T[ind,ind1_u0]=(1-ρ)*(1-α)*(1-pol_σ_E[ind1_e0])+ρ
                T[ind,ind1_e1]=(1-ρ)*α*pol_σ_E[ind1_e1]
                T[ind,ind1_u1]=(1-ρ)*α*(1-pol_σ_E[ind1_e1])


            elseif s_i==2
                ind1_e=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                ind1_u1=n_i*n_s*n_a*n_μ+i_i*n_a+a1_i
                ind1_u0=n_i*n_s*n_a*n_μ+a_i

                T[ind,ind1_e]=(1-ρ)*pol_σ_E[ind1_e]
                T[ind,ind1_u1]=(1-ρ)*(1-pol_σ_E[ind1_e])+ρ-δ
                T[ind,ind1_u0]=δ

            end


        elseif e_i==2

            if s_i==1
                indU=a_i
            elseif s_i==2
                indU=i_i*n_a+a_i
            end

            ind1_u=zeros(Int64,n_i)
            ind1_e=zeros(Int64,n_i)

            a1_i=pol_a_U[indU]

            if s_i==1

                for i1_i in 1:n_i
                    ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                    ind1_e[i1_i]=[i1_i-1,s_i-1,a1_i-1,pol_μ_U[a1_i,i1_i]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]

                    T[ind,ind1_e[i1_i]]=pol_σ_U[a1_i,i1_i]*p(θ[ind1_e[i1_i]])
                    T[ind,ind1_u[i1_i]]+=pol_σ_U[a1_i,i1_i]*(1-p(θ[ind1_e[i1_i]]))
                end



            elseif s_i==2

                for i1_i in 1:n_i
                    if i1_i==i_i
                        ind1_u[i1_i]=n_i*n_s*n_a*n_μ+i1_i*n_a+a1_i
                        ind1_e[i1_i]=[i1_i-1,s_i-1,a1_i-1,pol_μ_U[i1_i*n_a+a1_i,i1_i]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    else
                        ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                        ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    end

                    T[ind,ind1_e[i1_i]]=pol_σ_U[i_i*n_a+a1_i,i1_i]*p(θ[ind1_e[i1_i]])
                    T[ind,ind1_u[i1_i]]+=pol_σ_U[i_i*n_a+a1_i,i1_i]*(1-p(θ[ind1_e[i1_i]]))
                end
            end
        end
    end

    Φ=zeros(nstates)
    ind0=200
    Φ[ind0]=1.0

    for j in 1:1000
        Φ=T'*Φ
    end

    return Φ
end

# Computing aggregates in the stationary equilibrium

function ComputeAggregates(grids,pol_functions,Φ,z)

    p(θ)=min(m*θ^(1-ξ),1)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,θ=pol_functions

    nstates=n_i*n_s*n_a*n_μ+n_i*n_a

    # Compute unemployed and employed of each type
    E_I=zeros(n_i)
    E_E=zeros(n_i)
    U_E=zeros(n_i)
    for i_i in 1:n_i
        E_I[i_i]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+1:(i_i-1)*n_s*n_a*n_μ+n_μ*n_a])
        E_E[i_i]=sum(Φ[(i_i-1)*n_s*n_a*n_μ+n_μ*n_a+1:(i_i-1)*n_s*n_a*n_μ+2*n_μ*n_a])
        U_E[i_i]=sum(Φ[n_i*n_s*n_a*n_μ+i_i*n_a+1:n_i*n_s*n_a*n_μ+i_i*n_a+n_a])
    end
    U_I=sum(Φ[n_i*n_s*n_a*n_μ+1:n_i*n_s*n_a*n_μ+n_a])

    Y=0
    y=zeros(n_i)
    for i_i in 1:n_i
        y[i_i]=z[i_i]*E_I[i_i]^γ*E_E[i_i]^(1-γ)
        Y+=(1/n_i)*y[i_i]^((ν-1)/ν)
    end
    Y=Y^(ν/(ν-1))

    return Y,E_I,E_E,U_I,U_E
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

    ϵ=1e-4

    function Y_CES(z,I,E)
        Y=0
        for i_i in 1:n_i
            Y+=(1/n_i)*(z[i_i]*I[i_i]^γ*E[i_i]^(1-γ))^((ν-1)/ν)
        end
        Y=Y^(ν/(ν-1))
        return Y
    end

    pol_val_functions=false

    Iiter=0

    IL=0.411486102666538-0.01 #0.2
    IU=0.411486102666538+0.01 #0.35
    Ierr=1000

    Is=false
    Es=false
    U_E=false
    U_I=false
    Y=false
    Φ=false


    while Ierr>ϵ
        Iiter+=1
        println("I iter ",Iiter)

        Id=ones(n_i,1)*(IL*0.5+IU*0.5)

        Eiter=0

        EL=0.07609338860042157-0.01 # 0.05
        EU=min(0.07609338860042157+0.01,1/n_i-Id[1])  #min(1/n_i-Id[1],0.25)
        Eerr=1000

        while Eerr>ϵ && Eiter<10
            Eiter+=1
            println("E iter ",Eiter)

            Ed=ones(n_i,1)*(EL*0.5+EU*0.5)

            Y=Y_CES(z,Id,Ed)

            wages(Y,z,i,e)=((1/n_i)*γ*Y^(1/ν)*z^(1-(1/ν))*i^(γ*(1-(1/ν))-1)*e^((1-γ)*(1-(1/ν))),
                            (1/n_i)*(1-γ)*Y^(1/ν)*z^(1-(1/ν))*i^(γ*(1-(1/ν)))*e^((1-γ)*(1-(1/ν))-1))

            w=zeros(n_i,n_s)
            for i_i in 1:n_i
                for s_i in 1:n_s
                    w[i_i,s_i]=wages(Y,z[i_i],Id[i_i],Ed[i_i])[s_i]
                end
            end

            display(w)
            println("E=",Ed[1])
            println("I=",Id[1])

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

            Φ=ComputeDistribution(grids,pol_functions)
            Y,Is,Es,U_I,U_E=ComputeAggregates(grids,pol_functions,Φ,z)

            Eerr=Es[1]-Ed[1]

            if Eerr>0
                EL=Ed[1]
            else
                EU=Ed[1]
            end
            Eerr=abs(Eerr)
            println("E error is: ",Eerr)
        end

        Ierr=Is[1]-Id[1]

        if Ierr>0
            IL=Id[1]
        else
            IU=Id[1]
        end
        Ierr=abs(Ierr)
        println("I error is: ",Ierr)

    end

    return pol_val_functions,Φ,Y,Is,Es,U_I,U_E
end


# Computing the stationary distribution

function ComputeDistribution(grids,pol_functions)

    P(θ)=min(m*θ^(1-ξ),1)

    (grid_beq,grid_o,grid_e,grid_a)=grids
    n_beq,n_o,n_e,n_a=length(grid_beq),length(grid_o),length(grid_e),length(grid_a)
    (pol_a_Ei,pol_a_Ui,pol_σ_E,pol_σ_U,θ)=pol_functions

    nstates=n_beq*(n_o*n_e*n_a+n_a+n_o*(n_e-1)*n_a)
    nsvars=5
    ngrids_vars=[2,n_beq,n_o,n_e,n_a]

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


    # Construct Transition matrix

    i=Int64[]
    j=Int64[]
    k=Float64[]

    for ind in 1:nstates
        a_i=statestogrid[ind,5]
        e_i=statestogrid[ind,4]
        o_i=statestogrid[ind,3]
        beq_i=statestogrid[ind,2]
        s_i=statestogrid[ind,1]


        if s_i==1

            indE=[beq_i-1,o_i-1,e_i-1,a_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
            a1_i=pol_a_Ei[indE]

            ind1_ud=zeros(Int64,n_beq)

            if e_i==1
                ind1_en=[beq_i-1,o_i-1,e_i-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
                ind1_ep=[beq_i-1,o_i-1,(e_i+1)-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
                ind1_un=n_beq*n_o*n_e*n_a+(beq_i-1)*n_a+a1_i
                ind1_up=n_beq*n_o*n_e*n_a+n_beq*n_a+[beq_i-1,o_i-1,(e_i)-1,a1_i]'*[n_o*(n_e-1)*n_a,(n_e-1)*n_a,n_a,1]

                push!(i,ind)
                push!(j,ind1_en)
                push!(k,(1-δ)*(1-ρ)*(1-α[e_i])*pol_σ_E[ind1_en])

                push!(i,ind)
                push!(j,ind1_un)
                push!(k,(1-δ)*(1-ρ)*(1-α[e_i])*(1-pol_σ_E[ind1_en])+(1-δ)*ρ)

                push!(i,ind)
                push!(j,ind1_ep)
                push!(k,(1-δ)*(1-ρ)*α[e_i]*pol_σ_E[ind1_ep])

                push!(i,ind)
                push!(j,ind1_up)
                push!(k,(1-δ)*(1-ρ)*α[e_i]*(1-pol_σ_E[ind1_ep]))

                for beq1_i in 1:n_beq
                    ind1_ud[beq1_i]=n_beq*n_o*n_e*n_a+(beq1_i-1)*n_a+a1_i

                    push!(i,ind)
                    push!(j,ind1_ud[beq1_i])
                    push!(k,δ*weight_beq[beq1_i])
                end

            elseif e_i<n_e
                ind1_en=[beq_i-1,o_i-1,e_i-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
                ind1_ep=[beq_i-1,o_i-1,(e_i+1)-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
                ind1_un=n_beq*n_o*n_e*n_a+n_beq*n_a+[beq_i-1,o_i-1,(e_i-1)-1,a1_i]'*[n_o*(n_e-1)*n_a,(n_e-1)*n_a,n_a,1]
                ind1_up=n_beq*n_o*n_e*n_a+n_beq*n_a+[beq_i-1,o_i-1,(e_i)-1,a1_i]'*[n_o*(n_e-1)*n_a,(n_e-1)*n_a,n_a,1]

                push!(i,ind)
                push!(j,ind1_en)
                push!(k,(1-δ)*(1-ρ)*(1-α[e_i])*pol_σ_E[ind1_en])

                push!(i,ind)
                push!(j,ind1_un)
                push!(k,(1-δ)*(1-ρ)*(1-α[e_i])*(1-pol_σ_E[ind1_en])+(1-δ)*ρ)

                push!(i,ind)
                push!(j,ind1_ep)
                push!(k,(1-δ)*(1-ρ)*α[e_i]*pol_σ_E[ind1_ep])

                push!(i,ind)
                push!(j,ind1_up)
                push!(k,(1-δ)*(1-ρ)*α[e_i]*(1-pol_σ_E[ind1_ep]))

                for beq1_i in 1:n_beq
                    ind1_ud[beq1_i]=n_beq*n_o*n_e*n_a+(beq1_i-1)*n_a+a1_i

                    push!(i,ind)
                    push!(j,ind1_ud[beq1_i])
                    push!(k,δ*weight_beq[beq1_i])
                end

            elseif e_i==n_e
                ind1_e=[beq_i-1,o_i-1,e_i-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
                ind1_un=n_beq*n_o*n_e*n_a+n_beq*n_a+[beq_i-1,o_i-1,(e_i-1)-1,a1_i]'*[n_o*(n_e-1)*n_a,(n_e-1)*n_a,n_a,1]

                push!(i,ind)
                push!(j,ind1_e)
                push!(k,(1-δ)*(1-ρ)*pol_σ_E[ind1_e])

                push!(i,ind)
                push!(j,ind1_un)
                push!(k,(1-δ)*(1-ρ)*(1-pol_σ_E[ind1_e])+(1-δ)*ρ)

                for beq1_i in 1:n_beq
                    ind1_ud[beq1_i]=n_beq*n_o*n_e*n_a+(beq1_i-1)*n_a+a1_i

                    push!(i,ind)
                    push!(j,ind1_ud[beq1_i])
                    push!(k,δ*weight_beq[beq1_i])
                end
            end


        elseif s_i==2

            if e_i==1
                indU=(beq_i-1)*n_a+a_i
            else
                indU=n_beq*n_a+[beq_i-1,o_i-1,(e_i-1)-1,a_i]'*[n_o*(n_e-1)*n_a,(n_e-1)*n_a,n_a,1]
            end

            ind1_u=zeros(Int64,n_o)
            ind1_ud=zeros(Int64,n_beq)
            ind1_e=zeros(Int64,n_o)

            a1_i=pol_a_Ui[indU]

            if e_i==1

                for o1_i in 1:n_o
                    ind1_u[o1_i]=n_beq*n_o*n_e*n_a+(beq_i-1)*n_a+a1_i
                    ind1_e[o1_i]=[beq_i-1,o1_i-1,e_i-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]

                    push!(i,ind)
                    push!(j,ind1_e[o1_i])
                    push!(k,(1-δ)*pol_σ_U[ind1_u[o1_i]-nstates_E,o1_i]*P(θ[ind1_e[o1_i]]))

                    push!(i,ind)
                    push!(j,ind1_u[o1_i])
                    push!(k,(1-δ)*pol_σ_U[ind1_u[o1_i]-nstates_E,o1_i]*(1-P(θ[ind1_e[o1_i]])))
                end

            else

                for o1_i in 1:n_o
                    if o1_i==o_i
                        ind1_u[o1_i]=ind
                        ind1_e[o1_i]=[beq_i-1,o1_i-1,e_i-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
                    else
                        ind1_u[o1_i]=ind
                        ind1_e[o1_i]=[beq_i-1,o1_i-1,1-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
                    end

                    push!(i,ind)
                    push!(j,ind1_e[o1_i])
                    push!(k,(1-δ)*(1-χ[e_i])*pol_σ_U[ind1_u[o1_i]-nstates_E,o1_i]*P(θ[ind1_e[o1_i]]))

                    push!(i,ind)
                    push!(j,ind1_u[o1_i])
                    push!(k,(1-δ)*(1-χ[e_i])*pol_σ_U[ind1_u[o1_i]-nstates_E,o1_i]*(1-P(θ[ind1_e[o1_i]])))

                    if e_i==2

                        ind1_u[o1_i]=n_beq*n_o*n_e*n_a+(beq_i-1)*n_a+a1_i
                        ind1_e[o1_i]=[beq_i-1,o1_i-1,1-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]

                        push!(i,ind)
                        push!(j,ind1_e[o1_i])
                        push!(k,(1-δ)*χ[e_i]*pol_σ_U[ind1_u[o1_i]-nstates_E,o1_i]*P(θ[ind1_e[o1_i]]))

                        push!(i,ind)
                        push!(j,ind1_u[o1_i])
                        push!(k,(1-δ)*χ[e_i]*pol_σ_U[ind1_u[o1_i]-nstates_E,o1_i]*(1-P(θ[ind1_e[o1_i]])))

                    else
                        if o1_i==o_i
                            ind1_u[o1_i]=n_beq*n_o*n_e*n_a+n_beq*n_a+[beq_i-1,o_i-1,(e_i-2)-1,a1_i]'*[n_o*(n_e-1)*n_a,(n_e-1)*n_a,n_a,1]
                            ind1_e[o1_i]=[beq_i-1,o1_i-1,(e_i-1)-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
                        else
                            ind1_u[o1_i]=n_beq*n_o*n_e*n_a+n_beq*n_a+[beq_i-1,o_i-1,(e_i-2)-1,a1_i]'*[n_o*(n_e-1)*n_a,(n_e-1)*n_a,n_a,1]
                            ind1_e[o1_i]=[beq_i-1,o1_i-1,1-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]
                        end

                        push!(i,ind)
                        push!(j,ind1_e[o1_i])
                        push!(k,(1-δ)*χ[e_i]*pol_σ_U[ind1_u[o1_i]-nstates_E,o1_i]*P(θ[ind1_e[o1_i]]))

                        push!(i,ind)
                        push!(j,ind1_u[o1_i])
                        push!(k,(1-δ)*χ[e_i]*pol_σ_U[ind1_u[o1_i]-nstates_E,o1_i]*(1-P(θ[ind1_e[o1_i]])))
                    end
                end

                if o_i==2

                    ind1_u=ind
                    ind1_e=[beq_i-1,2-1,1-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]

                    push!(i,ind)
                    push!(j,ind1_e)
                    push!(k,(1-δ)*(1-χ[e_i])*(1-sum(pol_σ_U[ind1_u-nstates_E,:]))*P(θ[ind1_e]))

                    push!(i,ind)
                    push!(j,ind1_u)
                    push!(k,(1-δ)*(1-χ[e_i])*(1-sum(pol_σ_U[ind1_u-nstates_E,:]))*(1-P(θ[ind1_e])))

                    if e_i==2

                        ind1_u=n_beq*n_o*n_e*n_a+(beq_i-1)*n_a+a1_i
                        ind1_e=[beq_i-1,2-1,1-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]

                        push!(i,ind)
                        push!(j,ind1_e)
                        push!(k,(1-δ)*χ[e_i]*(1-sum(pol_σ_U[ind1_u-nstates_E,:]))*P(θ[ind1_e]))

                        push!(i,ind)
                        push!(j,ind1_u)
                        push!(k,(1-δ)*χ[e_i]*(1-sum(pol_σ_U[ind1_u-nstates_E,:]))*(1-P(θ[ind1_e])))

                    else
                        ind1_u=n_beq*n_o*n_e*n_a+n_beq*n_a+[beq_i-1,o_i-1,(e_i-2)-1,a1_i]'*[n_o*(n_e-1)*n_a,(n_e-1)*n_a,n_a,1]
                        ind1_e=[beq_i-1,2-1,1-1,a1_i]'*[n_o*n_e*n_a,n_e*n_a,n_a,1]

                        push!(i,ind)
                        push!(j,ind1_e)
                        push!(k,(1-δ)*χ[e_i]*(1-sum(pol_σ_U[ind1_u-nstates_E,:]))*P(θ[ind1_e]))

                        push!(i,ind)
                        push!(j,ind1_u)
                        push!(k,(1-δ)*χ[e_i]*(1-sum(pol_σ_U[ind1_u-nstates_E,:]))*(1-P(θ[ind1_e])))
                    end


                end

            end

            for beq1_i in 1:n_beq
                ind1_ud[beq1_i]=n_beq*n_o*n_e*n_a+(beq1_i-1)*n_a+a1_i

                push!(i,ind)
                push!(j,ind1_ud[beq1_i])
                push!(k,δ*weight_beq[beq1_i])
            end

        end
    end

    push!(i,nstates)
    push!(j,nstates)
    push!(k,0.0)

    T=sparse(i,j,k)

    Φ=zeros(nstates)
    ind0=n_beq*n_o*n_e*n_a+1:n_a:n_beq*n_o*n_e*n_a+n_beq*n_a
    Φ[ind0]=weight_beq

    for j in 1:50000
        Φ=T'*Φ
    end

    return Φ,T
end

# Computing aggregates in the stationary equilibrium

function ComputeAggregates(grids,Φ,z)

    (grid_beq,grid_o,grid_e,grid_a)=grids
    n_beq,n_o,n_e,n_a=length(grid_beq),length(grid_o),length(grid_e),length(grid_a)

    # Compute unemployed and employed of each type
    E=zeros(n_o,n_e)
    U=zeros(n_o,n_e)
    for o_i in 1:n_o

        for e_i in 1:n_e
            aux=0
            for beq_i in 1:n_beq
                aux+=sum(Φ[(beq_i-1)*n_o*n_e*n_a+(o_i-1)*n_e*n_a+(e_i-1)*n_a+1:(beq_i-1)*n_o*n_e*n_a+(o_i-1)*n_e*n_a+e_i*n_a])
            end
            E[o_i,e_i]=aux
            aux2=0
            if e_i==1
                if o_i==1
                    U[o_i,e_i]=sum(Φ[n_beq*n_o*n_e*n_a+1:n_beq*n_o*n_e*n_a+n_beq*n_a])
                end
            else
                for beq_i in 1:n_beq
                    aux2+=sum(Φ[n_beq*n_o*n_e*n_a+n_beq*n_a+(beq_i-1)*n_o*(n_e-1)*n_a+(o_i-1)*(n_e-1)*n_a+(e_i-2)*n_a+1:n_beq*n_o*n_e*n_a+n_beq*n_a+(beq_i-1)*n_o*(n_e-1)*n_a+(o_i-1)*(n_e-1)*n_a+(e_i-1)*n_a])
                end
                U[o_i,e_i]=aux2
            end
        end
    end

    Y=0.0
    y=zeros(n_o)
    for o_i in 1:n_o
        y[o_i]=z[o_i]*sum(γ.*(E[o_i,:].^ω))^(1/ω)
        Y+=ϕ[o_i]*y[o_i]^((ν-1)/ν)
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

    ϵ=5e-6

    function Y_CES(z,E)
        Y=0
        for o_i in 1:n_o
            Y+=ϕ[o_i]*(z[o_i]*prod(E[o_i,:].^γ))^((ν-1)/ν)
        end
        Y=Y^(ν/(ν-1))
        return Y
    end

    pol_val_functions=false

    El=[0.000539028 - 5e-6;
        0.000745483 - 5e-6;
        0.00104219 - 5e-6;
        0.00237116 - 5e-6]
    Eu=[0.000539028 + 5e-6;
        0.000745483 + 5e-6;
        0.00104219 + 5e-6;
        0.00237116 + 5e-6]

    Es=zeros(n_e)
    U=zeros(n_e)
    Y=false
    Φ=false
    Tr=false
    Eerr=1000*ones(n_e)

    Ed=zeros(n_o,n_e)
    k=1
    Eiter=0
    ready=false

    while ready==false
        Eiter+=1
        Ed[:,k]=[1;O-1]*(El[k]*0.5+Eu[k]*0.5)

        if k<n_e
            k=k+1
            continue
        end

        Y=Y_CES(z,Ed)

        prices(ϕ,Y,z,E,e_i)=ϕ*γ[e_i]*Y^(1/ν)*z^(1-(1/ν))*E[e_i]^(ω-1)*sum(γ.*(E.^ω))^((1/ω)*(1-(1/ν))-1)

        p=zeros(n_o,n_e)
        for o_i in 1:n_o
            for e_i in 1:n_e
                p[o_i,e_i]=prices(ϕ[o_i],Y,z[o_i],Ed[o_i,:],e_i)
            end
        end

        display(p)
        println("E=",Ed[1,:])

        pol_val_functions=multigrid(nGrids_a,p)
        (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions

        grid_a=LinRange(a_min,a_max,nGrids_a[end])
        grids=(grid_beq,grid_o,grid_e,grid_a)
        if n_anew!=nGrids_a[end]
            grid_a=LinRange(a_min,a_max,n_anew)
            pol_val_functions=transformVPolFunctions(pol_val_functions,grids,grid_a)
            (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_σ_E,pol_σ_U,J,θ)=pol_val_functions
            grids=(grid_beq,grid_o,grid_e,grid_a)
        end

        pol_a_Ei,pol_a_Ui=transformPola(pol_a_E,pol_a_U,grids)

        pol_functions=(pol_a_Ei,pol_a_Ui,pol_σ_E,pol_σ_U,θ)
        pol_val_functions=(V_E,V_U,W_E,W_U,pol_a_Ei,pol_a_Ui,pol_σ_E,pol_σ_U,J,θ)

        Φ,Tr=ComputeDistribution(grids,pol_functions)
        Y,Es,U=ComputeAggregates(grids,Φ,z)
        display(Es)

        while k>0
            Eerr[k]=Es[1,k]-Ed[1,k]
            println("E error for e_i=",k," is ",Eerr[k])
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

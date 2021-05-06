

function PlotResultsStatEq(grids,Stateq)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    (V_E,V_U,W_E,W_U,pol_a_E,pol_a_U,pol_μ_U,pol_σ_E,pol_σ_U,J,θ,Φ,Y,E_I,E_E,U_I,U_E)=StatEq

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
    statestogrid_U[end-2*n_a+1:end,3].=2
    statestogrid_U[:,4]=kron(ones(n_i+1,1),1:n_a)

    statestogrid=[statestogrid_E;statestogrid_U]

    #---------------------------------------------------------------------------
    # 1) Compute the policy functions of submarket choice of the unemployed
    #---------------------------------------------------------------------------

    plot1=plot(grid_μ[pol_μ_U[1:n_a]],title="Search strategy of unemployed",titlefontsize=7)
    for state in 2:1+n_i
        plot1=plot!(grid_μ[pol_μ_U[(state-1)*n_a+1:state*n_a]])
    end
    display(plot1)

    #---------------------------------------------------------------------------
    # 2) Mean wealth of workers in stationary equilibrium
    #---------------------------------------------------------------------------

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
    display(a_bar)


    #---------------------------------------------------------------------------
    # 3) Construct a wealth histogram in stationary equilibrium
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

    plot2=plot(grid_a,Dist_a,title="Distribution of wealth",titlefontsize=7)
    display(plot2)

end

function PlotResultsTransition(grids,zt,pol_val_results,aggregates_transition)

    (grid_i,grid_s,grid_a,grid_μ)=grids
    n_i,n_s,n_a,n_μ=length(grid_i),length(grid_s),length(grid_a),length(grid_μ)
    (V_E,V_U,W_E,W_U,pol_a_Ei,pol_a_Ui,pol_μ_U,pol_σ_E,pol_σ_U,J,θ,Φ)=pol_val_results
    (I,E,U_I,U_E)=aggregates_transition
    T=size(V_E,2)
    newzt=zt[:,end]*ones(1,T)
    newzt[:,1:size(zt,2)]=zt

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
    statestogrid_U[end-2*n_a+1:end,3].=2
    statestogrid_U[:,4]=kron(ones(n_i+1,1),1:n_a)

    statestogrid=[statestogrid_E;statestogrid_U]



    ########################################################
    # 1) Evolution of aggregates in transition
    ########################################################
    x=1:T
    plot1=plot(x,NewI',title="Panel A: Inexperienced employed",titlefontsize=7,legend=false,linestyle=[:solid :dash :dashdot :solid :dash :dashdot :solid :dash :dashdot],lc=[:black :red :blue :black :red :blue :black :red :blue])
    plot2=plot(x,NewE',title="Panel B: Experienced employed",titlefontsize=7,legend=false,linestyle=[:solid :dash :dashdot :solid :dash :dashdot :solid :dash :dashdot],lc=[:black :red :blue :black :red :blue :black :red :blue])
    plot3=plot(x,newzt',title="Panel C: productivity shock",titlefontsize=7,legend=false,linestyle=[:solid :dash :dashdot :solid :dash :dashdot :solid :dash :dashdot],lc=[:black :red :blue :black :red :blue :black :red :blue])

    plot4=plot(plot1,plot2,plot3,layout=(1,3))
    display(plot4)


    ########################################################
    # 2) Evolution of wealth in transition
    ########################################################

    a_bar=zeros(n_i,n_s,T)
    for t in 1:T
        aux=zeros(Int64,n_μ*n_a,1)
        for i_i in 1:n_i
            for s_i in 1:n_s
                ind1=[i_i-1,s_i-1,1-1,1]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                indend=[i_i-1,s_i-1,n_a-1,n_μ]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                aux[:]=kron(1:n_a,ones(n_μ,1))
                a_bar[i_i,s_i,t]=sum(Φ[ind1:indend,t].*grid_a[aux])/sum(Φ[ind1:indend,t])
            end
        end
    end

    plot5=plot(a_bar[1,1,:],title="Average wealth of inexperienced workers",titlefontsize=7)
    plot!(a_bar[2,1,:])
    display(plot5)

    plot6=plot(a_bar[1,2,:],title="Average wealth of experienced workers",titlefontsize=7)
    plot!(a_bar[2,2,:])
    display(plot6)


    ###############################################################
    # 3) Evolution of variables in transition for different agents
    ###############################################################
    N=50000
    i_agent1,s_agent1,a_agent1,μ_agent1=1,2,1,95
    i_agent2,s_agent2,a_agent2,μ_agent2=1,2,50,95

    ind_agent1=zeros(Int64,N,T+1)
    ind_agent2=zeros(Int64,N,T+1)
    moving1=zeros(Int64,N,T)
    moving2=zeros(Int64,N,T)
    prob_moving1=zeros(T)
    prob_moving2=zeros(T)

    ind_agent1[:,1].=[i_agent1-1,s_agent1-1,a_agent1-1,μ_agent1]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
    ind_agent2[:,1].=[i_agent2-1,s_agent2-1,a_agent2-1,μ_agent2]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]


    for t in 1:T
        i=Int64[]
        j=Int64[]
        k=Float64[]
        for ind in 1:nstates
            μ_i=statestogrid[ind,5]
            a_i=statestogrid[ind,4]
            s_i=statestogrid[ind,3]
            i_i=statestogrid[ind,2]
            e_i=statestogrid[ind,1]

            if e_i==1

                indE=[i_i-1,s_i-1,a_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                a1_i=pol_a_Ei[indE,t]
                if s_i==1
                    ind1_e0=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    ind1_e1=[i_i-1,2-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    ind1_u0=n_i*n_s*n_a*n_μ+a1_i
                    ind1_u1=n_i*n_s*n_a*n_μ+i_i*n_a+a1_i

                    push!(i,ind)
                    push!(j,ind1_e0)
                    push!(k,(1-ρ)*(1-α)*pol_σ_E[ind1_e0,t])

                    push!(i,ind)
                    push!(j,ind1_u0)
                    push!(k,(1-ρ)*(1-α)*(1-pol_σ_E[ind1_e0,t])+ρ)

                    push!(i,ind)
                    push!(j,ind1_e1)
                    push!(k,(1-ρ)*α*pol_σ_E[ind1_e1,t])

                    push!(i,ind)
                    push!(j,ind1_u1)
                    push!(k,(1-ρ)*α*(1-pol_σ_E[ind1_e1,t]))

                elseif s_i==2
                    ind1_e=[i_i-1,s_i-1,a1_i-1,μ_i]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                    ind1_u1=n_i*n_s*n_a*n_μ+i_i*n_a+a1_i
                    ind1_u0=n_i*n_s*n_a*n_μ+a_i

                    push!(i,ind)
                    push!(j,ind1_e)
                    push!(k,(1-ρ)*pol_σ_E[ind1_e,t])

                    push!(i,ind)
                    push!(j,ind1_u1)
                    push!(k,(1-ρ)*(1-pol_σ_E[ind1_e,t])+ρ-δ)

                    push!(i,ind)
                    push!(j,ind1_u0)
                    push!(k,δ)
                end


            elseif e_i==2

                if s_i==1
                    indU=a_i
                elseif s_i==2
                    indU=i_i*n_a+a_i
                end

                ind1_u=zeros(Int64,n_i)
                ind1_e=zeros(Int64,n_i)

                a1_i=pol_a_Ui[indU,t]

                if s_i==1

                    for i1_i in 1:n_i
                        ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                        ind1_e[i1_i]=[i1_i-1,s_i-1,a1_i-1,pol_μ_U[a1_i,i1_i,t]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]

                        push!(i,ind)
                        push!(j,ind1_e[i1_i])
                        push!(k,pol_σ_U[a1_i,i1_i,t]*p(θ[ind1_e[i1_i],t]))

                        push!(i,ind)
                        push!(j,ind1_u[i1_i])
                        push!(k,pol_σ_U[a1_i,i1_i,t]*(1-p(θ[ind1_e[i1_i],t])))
                    end

                elseif s_i==2

                    for i1_i in 1:n_i
                        if i1_i==i_i
                            ind1_u[i1_i]=n_i*n_s*n_a*n_μ+i1_i*n_a+a1_i
                            ind1_e[i1_i]=[i1_i-1,s_i-1,a1_i-1,pol_μ_U[i1_i*n_a+a1_i,i1_i,t]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        else
                            ind1_u[i1_i]=n_i*n_s*n_a*n_μ+a1_i
                            ind1_e[i1_i]=[i1_i-1,1-1,a1_i-1,pol_μ_U[a1_i,i1_i,t]]'*[n_s*n_a*n_μ,n_a*n_μ,n_μ,1]
                        end

                        push!(i,ind)
                        push!(j,ind1_e[i1_i])
                        push!(k,pol_σ_U[i_i*n_a+a1_i,i1_i,t]*p(θ[ind1_e[i1_i],t]))

                        push!(i,ind)
                        push!(j,ind1_u[i1_i])
                        push!(k,pol_σ_U[i_i*n_a+a1_i,i1_i,t]*(1-p(θ[ind1_e[i1_i],t])))
                    end
                end
            end
        end

        push!(i,nstates)
        push!(j,nstates)
        push!(k,0.0)

        Tr=sparse(i,j,k)

        for s in 1:N
            weights=pweights(Tr[ind_agent1[s,t],:])
            ind_agent1[s,t+1]=sample(weights)
            s_i=statestogrid[ind_agent1[s,t+1],3]
            if t==1
                if s_i==1
                    moving1[s,t]=1
                end
            elseif t>1
                if s_i==1 && sum(moving1[s,1:t-1])==0
                    moving1[s,t]=1
                end
            end
            weights2=pweights(Tr[ind_agent2[s,t],:])
            ind_agent2[s,t+1]=sample(weights2)
            s_i=statestogrid[ind_agent2[s,t+1],3]
            if t==1
                if s_i==1
                    moving2[s,t]=1
                end
            elseif t>1
                if s_i==1 && sum(moving2[s,1:t-1])==0
                    moving2[s,t]=1
                end
            end
        end
        if t==1
            prob_moving1[t]=sum(moving1[:,t])/N
            prob_moving2[t]=sum(moving2[:,t])/N
        else
            prob_moving1[t]=sum(moving1[:,t])/(N-sum(moving1[:,1:t-1]))
            prob_moving2[t]=sum(moving2[:,t])/(N-sum(moving1[:,1:t-1]))
        end
        println("t=",t," complete")
    end

    plot7=plot(prob_moving1,title="Probability of moving sector",titlefontsize=7)
    plot!(prob_moving2)

end

@everywhere function bisection_derivative(F,x_l,x_u;initialguess=0,weight=0.5,ϵ=1e-6)

    dx=1e-6
    f(x)=(F(x+dx)-F(x))/dx
    iter=0
    error=1000

    # finding if function is upward or downward sloping
    if f(x_l)>0 && f(x_u)<0 #downward sloping
        type=1
    elseif f(x_l)<0 && f(x_u)>0 # upward sloping
        type=2
    elseif isnan(f(x_u))

        if f(x_l)>0
            type=1
            iteru=0
            done=false
            while done==false
                iteru+=1
                if isnan(f(x_u))
                    x_u=x_u-(x_u-x_l)/2
                    done=false
                elseif f(x_u)>0
                    x_u=x_u+(x_u-x_l)/2
                    done=false
                elseif f(x_u)<0
                    done=true
                end
                if iteru>=1000
                    println("error: no guarantee of a 0 in bisection")
                    return
                end
            end
        elseif f(x_l)<0
            type=2
            iteru=0
            done=false
            while done==false
                iteru+=1
                if isnan(f(x_u))
                    x_u=x_u-(x_u-x_l)/2
                    done=false
                elseif f(x_u)<0
                    x_u=x_u+(x_u-x_l)/2
                    done=false
                elseif f(x_u)>0
                    done=true
                end
                if iteru>=1000
                    println("error: no guarantee of a 0 in bisection")
                    return
                end
            end
        end

    else
        println("error: no guarantee of a 0 in bisection")
        return
    end

    if initialguess==0
        x_0=weight*x_l+(1-weight)*x_u
    else
        x_0=initialguess
    end

    x_m=0.0
    while error>ϵ
        iter=iter+1
        if iter==1
            x_m=x_0
        else
            x_m=weight*x_l+(1-weight)*x_u
        end
        error=f(x_m)

        if type==1 # Downward sloping
            if error>0
                x_l=x_m
            else
                x_u=x_m
            end
        elseif type==2 # Upward sloping
            if error>0
                x_u=x_m
            else
                x_l=x_m
            end
        end
        error=abs(error)
    end
    return x_m
end

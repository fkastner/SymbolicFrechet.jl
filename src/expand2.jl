function expand2(x)
    !iscall(x) && return x
    
    op = operation(x)
    args = arguments(x)
    op_expanded = expand2(op)
    args_expanded = map(expand2, args)
    expand2(op_expanded, args_expanded)
end

expand2(op, args) = op(args...)

function expand2(op::typeof(*), args)
    isempty(args) && return 1
    for (i,arg) in enumerate(args)
        iscall(arg) || continue
        inner_op = operation(arg)
        inner_args = arguments(arg)

        if inner_op == (+)
            # rest_prod = prod(args[(1:i-1)∪(i+1:end)])
            rest_prod = expand2(*, args[!=(i).(eachindex(args))])
            return sum(inner_args) do a
                # expand2(a*rest_prod)
                expand2(*, (a,rest_prod))
            end
        end
    end
    return prod(args)
end

function expand2(op::AbstractMultiLinearOperator, args)
    for (i,arg) in enumerate(args)
        iscall(arg) || continue
        inner_op = operation(arg)
        inner_args = arguments(arg)

        if inner_op == (+)
            return sum(inner_args) do a
                new_args = copy(args)
                new_args[i] = a
                expand2(op, new_args)
            end
        elseif inner_op == (*)
            for (j,iarg) in enumerate(inner_args)
                if isscalar(iarg)
                    new_args = copy(args)
                    # new_args[i] = prod(inner_args[(1:j-1)∪(j+1:end)])
                    new_args[i] = prod(inner_args[!=(j).(eachindex(inner_args))])
                    # return expand2(iarg*expand2(op, new_args))
                    # return expand2(iarg*op(new_args...))
                    return expand2(*, (iarg, expand2(op, new_args)))
                end
                # new_iarg = expand2(iarg)
                # if !isequal(new_iarg, iarg)
                #     new_iargs = copy(inner_args)
                #     new_iargs[j] = new_iarg
                #     new_args = copy(args)
                #     new_args[i] = expand2(prod(new_iargs))
                #     return expand2(op, new_args)
                # end
            end
        end
    end

    return op(args...)
end

function expand2(op::FrechetDerivative, args)
    order = op.order
    fun = only(args)
    
    # @show op, args
    # return op(fun)

    isscalar(fun) && return Returns(0)
    !iscall(fun) && return op(fun)

    inner_op = operation(fun)
    inner_args = arguments(fun)

    if inner_op == (+)
        # linearity
        return (args...) -> sum(inner_args) do arg
            # expand2(FrechetDifferential(order, arg), collect(args))
            expand2(FrechetDifferential(order, arg))(args...)
        end
    elseif isa(inner_op, AbstractMultiLinearOperator) || inner_op == (*)
        # chain rule; don't ask
        len = length(inner_args)
        return (args...) -> sum(Iterators.product(ntuple(Returns(0:len),order)...)) do t
            new_op = inner_op
            new_args = similar(inner_args)

            c = count(==(0), t)
            new_op = (c == 0) ? inner_op : expand2(FrechetDifferential(c, inner_op))
            (c>0) && append!(new_args, args[collect(==(0).(t))])

            for i = 1:len
                c = count(==(i), t)
                # new_args[i] = (c == 0) ? inner_args[i] : expand2(expand2(FrechetDifferential(c, inner_args[i])), collect(args[collect(t) .== i]))
                new_args[i] = (c == 0) ? inner_args[i] : expand2(FrechetDifferential(c, inner_args[i]))(args[collect(==(i).(t))]...)
            end

            expand2(new_op, new_args)
        end
    elseif inner_op == (^) && isa(inner_args[2], Number)
        # chain rule, again
        base = inner_args[1]
        len = inner_args[2]
        return (args...) -> sum(Iterators.product(ntuple(Returns(1:len),order)...)) do t
            # new_args = Vector{Any}(undef, len)
            # for i = 1:len
            #     c = count(==(i), t)
            #     # new_args[i] = (c == 0) ? base : expand2(FrechetDifferential(c, base), collect(args[collect(t) .== i]))
            #     new_args[i] = (c == 0) ? base : FrechetDifferential(c, base)(args[collect(==(i).(t))]...)
            # end
            # expand2(prod(new_args))

            expand2(prod(1:len) do i
                c = count(==(i), t)
                (c == 0) ? base : FrechetDifferential(c, base)(args[collect(==(i).(t))]...)
            end)
        end
    end

    return op(fun)
end
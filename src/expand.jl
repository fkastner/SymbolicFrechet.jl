function expand(x)
    !iscall(x) && return x
    
    op = operation(x)
    args = arguments(x)
    op_expanded = expand(op)
    args_expanded = map(expand, args)
    expand(op_expanded, args_expanded)
end

expand(op, args) = op(args...)

function expand(op::typeof(*), args)
    
    # helper(x) = SymbolicUtils.is_operation(+)(x) ? arguments(x) : Any[x]
    # return sum(prod.(Iterators.product(helper.(args)...)))
    
    
    isempty(args) && return 1
    for (i,arg) in enumerate(args)
        iscall(arg) || continue
        inner_op = operation(arg)
        inner_args = arguments(arg)

        if inner_op == (+)
            # rest_prod = prod(args[(1:i-1)∪(i+1:end)])
            rest_prod = expand(*, args[!=(i).(eachindex(args))])
            return sum(inner_args) do a
                # expand(a*rest_prod)
                expand(*, (a,rest_prod))
            end
        end
    end
    return prod(args)
end

function expand(op::typeof(^), args)
    base = args[1]
    exponent = args[2]

    isinteger(exponent) && exponent > 0 || return base^exponent
    iscall(base) && operation(base) == (+) || return base^exponent

    inner_args = arguments(base)
    n = length(inner_args)
    return sum(multiexponents(n, exponent)) do exponents
        multinomial(exponents...)*prod(inner_args.^exponents)
    end
end

function expand(op::AbstractMultiLinearOperator, args)
    for (i,arg) in enumerate(args)
        iscall(arg) || continue
        inner_op = operation(arg)
        inner_args = arguments(arg)

        if inner_op == (+)
            return sum(inner_args) do a
                new_args = copy(args)
                new_args[i] = a
                expand(op, new_args)
            end
        elseif inner_op == (*)
            for (j,iarg) in enumerate(inner_args)
                if isscalar(iarg)
                    new_args = copy(args)
                    # new_args[i] = prod(inner_args[(1:j-1)∪(j+1:end)])
                    new_args[i] = prod(inner_args[!=(j).(eachindex(inner_args))])
                    # return expand(iarg*expand(op, new_args))
                    # return expand(iarg*op(new_args...))
                    return expand(*, (iarg, expand(op, new_args)))
                end
                # new_iarg = expand(iarg)
                # if !isequal(new_iarg, iarg)
                #     new_iargs = copy(inner_args)
                #     new_iargs[j] = new_iarg
                #     new_args = copy(args)
                #     new_args[i] = expand(prod(new_iargs))
                #     return expand(op, new_args)
                # end
            end
        end
    end

    return op(args...)
end

function expand(op::FrechetDifferential, args)
    order = op.order
    fun = op.fun

    isscalar(fun) && return 0
    !iscall(fun) && return @invoke expand(op::AbstractMultiLinearOperator, args)

    inner_op = operation(fun)
    inner_args = arguments(fun)

    if inner_op == (+)
        # linearity
        return sum(inner_args) do arg
            # expand(FrechetDifferential(order, arg)(args...))
            expand(FrechetDifferential(order, arg), args)
        end
    elseif isa(inner_op, AbstractMultiLinearOperator) || inner_op == (*)
        # chain rule; don't ask
        len = length(inner_args)
        return sum(Iterators.product(ntuple(Returns(0:len),order)...)) do t
            new_op = inner_op
            new_args = similar(inner_args)

            c = count(==(0), t)
            new_op = (c == 0) ? inner_op : FrechetDifferential(c, inner_op)
            (c>0) && append!(new_args, args[collect(==(0).(t))])

            for i = 1:len
                c = count(==(i), t)
                # new_args[i] = (c == 0) ? inner_args[i] : expand(expand(FrechetDifferential(c, inner_args[i])), args[collect(t) .== i])
                # new_args[i] = (c == 0) ? inner_args[i] : expand(FrechetDifferential(c, inner_args[i])(args[collect(==(i).(t))]...))
                new_args[i] = (c == 0) ? inner_args[i] : expand(FrechetDifferential(c, inner_args[i]), args[collect(==(i).(t))])
            end

            expand(new_op, new_args)
        end
    elseif inner_op == (^) && isa(inner_args[2], Number)
        # chain rule, again
        base = inner_args[1]
        len = inner_args[2]
        return sum(Iterators.product(ntuple(Returns(1:len),order)...)) do t
            # new_args = Vector{Any}(undef, len)
            # for i = 1:len
            #     c = count(==(i), t)
            #     # new_args[i] = (c == 0) ? base : expand(FrechetDifferential(c, base), args[collect(t) .== i])
            #     new_args[i] = (c == 0) ? base : FrechetDifferential(c, base)(args[collect(==(i).(t))]...)
            # end
            # expand(prod(new_args))

            expand(prod(1:len) do i
                c = count(==(i), t)
                (c == 0) ? base : FrechetDifferential(c, base)(args[collect(==(i).(t))]...)
            end)
        end
    end

    return @invoke expand(op::AbstractMultiLinearOperator, args)
end

# function expand(op::FrechetDerivative, args)
#     order = op.order
#     fun = only(args)
    
#     # @show op, args
#     # return op(fun)

#     isscalar(fun) && return Returns(0)
#     !iscall(fun) && return op(fun)

#     inner_op = operation(fun)
#     inner_args = arguments(fun)

#     if inner_op == (+)
#         # linearity
#         return (args...) -> sum(inner_args) do arg
#             # expand(FrechetDifferential(order, arg), collect(args))
#             expand(FrechetDifferential(order, arg)(args...))
#         end
#     elseif isa(inner_op, AbstractMultiLinearOperator) || inner_op == (*)
#         # chain rule; don't ask
#         len = length(inner_args)
#         return (args...) -> sum(Iterators.product(ntuple(Returns(0:len),order)...)) do t
#             new_op = inner_op
#             new_args = similar(inner_args)

#             c = count(==(0), t)
#             new_op = (c == 0) ? inner_op : expand(FrechetDifferential(c, inner_op))
#             (c>0) && append!(new_args, args[collect(==(0).(t))])

#             for i = 1:len
#                 c = count(==(i), t)
#                 # new_args[i] = (c == 0) ? inner_args[i] : expand(expand(FrechetDifferential(c, inner_args[i])), collect(args[collect(t) .== i]))
#                 new_args[i] = (c == 0) ? inner_args[i] : expand(FrechetDifferential(c, inner_args[i])(args[collect(==(i).(t))]...))
#             end

#             expand(new_op, new_args)
#         end
#     elseif inner_op == (^) && isa(inner_args[2], Number)
#         # chain rule, again
#         base = inner_args[1]
#         len = inner_args[2]
#         return (args...) -> sum(Iterators.product(ntuple(Returns(1:len),order)...)) do t
#             new_args = Vector{Any}(undef, len)
#             for i = 1:len
#                 c = count(==(i), t)
#                 # new_args[i] = (c == 0) ? base : expand(FrechetDifferential(c, base), collect(args[collect(t) .== i]))
#                 new_args[i] = (c == 0) ? base : FrechetDifferential(c, base)(args[collect(==(i).(t))]...)
#             end
#             expand(prod(new_args))
#         end
#     end

#     return op(fun)
# end
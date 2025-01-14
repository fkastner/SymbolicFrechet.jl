module SymbolicFrechet

using SymbolicUtils
import SymbolicUtils: Symbolic

export MultiLinearOperator, FrechetDifferential
export expand_fdiffs, expand_MLOs

abstract type AbstractMultiLinearOperator end

(T::AbstractMultiLinearOperator)(args...) = length(args) == nargs(T) ? term(T, args...) : throw(ArgumentError("need $(nargs(T)) many arguments, got $(length(args))"))
SymbolicUtils.promote_symtype(::AbstractMultiLinearOperator, Ts...) = Real
SymbolicUtils.isbinop(::AbstractMultiLinearOperator) = false

function expand_MLOs(O::Symbolic)
    iscall(O) || return O

    op = operation(O)
    args = arguments(O)

    if isa(op, AbstractMultiLinearOperator)
        return expand_MLO(op, args)
    else
        args_expanded = map(expand_MLOs, args)
        return op(args_expanded...)
    end
end
expand_MLOs(x) = x

function expand_MLO(T, args)
    for (i,arg) in enumerate(args)
        iscall(arg) || continue
        op = operation(arg)

        if op == (+)
            inner_args = arguments(arg)
            return sum(inner_args) do a
                new_args = copy(args)
                new_args[i] = a
                expand_MLO(T,new_args)
            end
        elseif op == (*)
            inner_args = arguments(arg)
            if isa(inner_args[1], Number)
                new_args = copy(args)
                new_args[i] = prod(inner_args[2:end])
                return inner_args[1]*expand_MLO(T,new_args)
            end
        end
    end

    return T(args...)
end

struct MultiLinearOperator <: AbstractMultiLinearOperator
    name
    nargs
end
nargs(T::MultiLinearOperator) = T.nargs
Base.nameof(T::MultiLinearOperator) = T.name
Base.show(io::IO, T::MultiLinearOperator) = print(io, T.name)


struct FrechetDifferential <: AbstractMultiLinearOperator
    order
    fun
end
FrechetDifferential(order=1) = fun -> FrechetDifferential(order, fun)
FrechetDifferential(order,d::FrechetDifferential) = FrechetDifferential(order+d.order,d.fun)

nargs(D::FrechetDifferential) = isa(D.fun, AbstractMultiLinearOperator) ? D.order + nargs(D.fun) : D.order

SymbolicUtils.operation(D::FrechetDifferential) = FrechetDifferential
SymbolicUtils.arguments(D::FrechetDifferential) = D.fun

Base.nameof(D::FrechetDifferential) = Symbol("d($(D.fun))")
Base.show(io::IO, D::FrechetDifferential) = print(io, (D.order == 1 ? "d" : "d^$(D.order)") * "($(D.fun))")

function expand_fdiffs(O::Symbolic)
    iscall(O) || return O

    op = operation(O)
    args = arguments(O)

    if isa(op, FrechetDifferential)
        return expand_fdiff(op, args)
    else
        args_expanded = map(expand_fdiffs, args)
        return op(args_expanded...)
    end
end
expand_fdiffs(x) = x

function expand_fdiff(T, args)
    order = T.order
    fun = T.fun

    isa(fun, Number) && return 0
    !iscall(fun) && return T(args...)

    op = operation(fun)
    inner_args = arguments(fun)

    if op == (+)
        return sum(inner_args) do arg
            expand_fdiff(FrechetDifferential(order, arg), args)
        end
    elseif isa(op, AbstractMultiLinearOperator)
        # chain rule; don't ask
        len = length(inner_args)
        return sum(Iterators.product(ntuple(Returns(0:len),order)...)) do t
            new_op = op
            new_args = similar(inner_args)

            c = count(==(0), t)
            new_op = (c == 0) ? op : FrechetDifferential(c, op)
            (c>0) && append!(new_args, args[collect(t) .== 0])

            for i = 1:len
                c = count(==(i), t)
                new_args[i] = (c == 0) ? inner_args[i] : FrechetDifferential(c, inner_args[i])(args[collect(t) .== i]...)
            end

            new_op(new_args...)
        end

        # if order == 1
        #     return FrechetDifferential(1,op)(inner_args...,args...) + sum(1:len) do i
        #         new_args = copy(inner_args)
        #         new_args[i] = expand_fdiff(FrechetDifferential(1,inner_args[i]),args)
        #         return op(new_args...)
        #     end
        # elseif order == 2
        #     res = 0
        #     for i=0:len
        #         for j=0:len
        #             new_op = op
        #             new_args = copy(inner_args)
        #             copy_args = copy(args)
        #             if i == 0
        #                 new_op = FrechetDifferential(1,new_op)
        #                 push!(new_args, popfirst!(copy_args))
        #             elseif i == j # != 0
        #                 new_args[i] = FrechetDifferential(2,new_args[i])(copy_args...)
        #             else
        #                 # new_args[i] = FrechetDifferential(1,new_args[i])(popfirst!(copy_args))
        #                 new_args[i] = expand_fdiff(FrechetDifferential(1, new_args[i]), [popfirst!(copy_args)])
        #             end
        #             if j == 0
        #                 new_op = FrechetDifferential(1,new_op)
        #                 push!(new_args, popfirst!(copy_args))
        #             elseif i != j
        #                 # new_args[j] = FrechetDifferential(1,new_args[j])(popfirst!(copy_args))
        #                 new_args[j] = expand_fdiff(FrechetDifferential(1, new_args[j]), [popfirst!(copy_args)])
        #             end
        #             res += new_op(new_args...)
        #         end
        #     end
        #     return res
        # else
        #     # too much chain rule :(
        #     # return expand_fdiff(FrechetDifferential(order-1, expand_fdiff(FrechetDifferential(1, fun), args[[1]])), args[2:end]) 
        #     throw(ArgumentError("too much order: $order"))
        # end
    end

    return T(args...)
end

end # module
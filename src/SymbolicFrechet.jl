module SymbolicFrechet

using SymbolicUtils
import SymbolicUtils: Symbolic, symtype

using Combinatorics

export MultiLinearOperator, FrechetDerivative
export expand_fdiffs, expand_MLOs

abstract type AbstractMultiLinearOperator end

function (T::AbstractMultiLinearOperator)(args...)
    length(args) != nargs(T) && throw(ArgumentError("need $(nargs(T)) many arguments, got $(length(args))"))
    issymmetricoperator(T) ? term(T, sort!(collect(args); by=hash)...) : term(T, args...)
end
issymmetricoperator(T::AbstractMultiLinearOperator) = false

SymbolicUtils.promote_symtype(::AbstractMultiLinearOperator, Ts...) = promote_type(Ts...)
SymbolicUtils.isbinop(::AbstractMultiLinearOperator) = false


struct ZeroDeriv <: Number end
isscalar(x) = iscall(x) ? isscalar(operation(x)) && all(isscalar, arguments(x)) : false
isscalar(::typeof(^)) = true
isscalar(::typeof(+)) = true
isscalar(::typeof(*)) = true
isscalar(::Number) = true
isscalar(::SymbolicUtils.BasicSymbolic{ZeroDeriv}) = true

#### distribute * over +
function distribute(O::Symbolic)
    iscall(O) || return O

    op = operation(O)
    args = arguments(O)

    if op == (*)
        for (i,arg) in enumerate(args)
            iscall(arg) || continue
            inner_op = operation(arg)
            inner_args = arguments(arg)
            if inner_op == (+)
                rest_prod = prod(args[(1:i-1)∪(i+1:end)])
                return sum(inner_args) do a
                    distribute(a*rest_prod)
                end
            end
        end
        return O
    else
        args_expanded = map(distribute, args)
        return op(args_expanded...)
    end
end
distribute(x) = x 

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
        inner_args = arguments(arg)

        if op == (+)
            return sum(inner_args) do a
                new_args = copy(args)
                new_args[i] = a
                expand_MLO(T,new_args)
            end
        elseif op == (*)
            for (j,iarg) in enumerate(inner_args)
                if isscalar(iarg)
                    new_args = copy(args)
                    new_args[i] = prod(inner_args[(1:j-1)∪(j+1:end)])
                    return distribute(iarg*expand_MLO(T,new_args))
                end
                new_iarg = expand_MLOs(iarg)
                if !isequal(new_iarg, iarg)
                    new_iargs = copy(inner_args)
                    new_iargs[j] = new_iarg
                    new_args = copy(args)
                    new_args[i] = distribute(prod(new_iargs))
                    return expand_MLO(T, new_args)
                end
            end
        elseif isa(op, AbstractMultiLinearOperator)
            new_args = copy(args)
            new_args[i] = expand_MLO(op, inner_args)
            !isequal(new_args[i], args[i]) && return expand_MLO(T, new_args)
        end
    end

    return T(args...)
end

struct MultiLinearOperator <: AbstractMultiLinearOperator
    name
    nargs::Int
    symmetric::Bool
    MultiLinearOperator(name,nargs,symmetric=false) = new(name,nargs,symmetric)
end
nargs(T::MultiLinearOperator) = T.nargs
issymmetricoperator(T::MultiLinearOperator) = T.symmetric

Base.nameof(T::MultiLinearOperator) = T.name
Base.show(io::IO, T::MultiLinearOperator) = print(io, T.name)

struct FrechetDifferential <: AbstractMultiLinearOperator
    order
    fun
end
FrechetDifferential(order,d::FrechetDifferential) = FrechetDifferential(order+d.order,d.fun)
FrechetDifferential(order, ::typeof(*)) = Returns(0)
FrechetDifferential(order, ::SymbolicUtils.BasicSymbolic{ZeroDeriv}) = Returns(0)

nargs(D::FrechetDifferential) = isa(D.fun, AbstractMultiLinearOperator) ? D.order + nargs(D.fun) : D.order
issymmetricoperator(::FrechetDifferential) = true

# TermInterface.jl interface
SymbolicUtils.isexpr(x::FrechetDifferential) = true
SymbolicUtils.head(x::FrechetDifferential) = FrechetDerivative(D.order)
SymbolicUtils.children(x::FrechetDifferential) = Any[D.fun]
SymbolicUtils.iscall(x::FrechetDifferential) = true
SymbolicUtils.operation(D::FrechetDifferential) = FrechetDerivative(D.order)
SymbolicUtils.arguments(D::FrechetDifferential) = Any[D.fun]
SymbolicUtils.maketerm(::Type{FrechetDifferential}, f, args, metadata) = f(args...)

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

expand_fdiff(T, args) = T(args...)
function expand_fdiff(T::FrechetDifferential, args)
    order = T.order
    fun = T.fun

    isscalar(fun) && return 0
    !iscall(fun) && return T(args...)

    op = operation(fun)
    inner_args = arguments(fun)

    if op == (+)
        # linearity
        return sum(inner_args) do arg
            expand_fdiff(FrechetDifferential(order, arg), args)
        end
    elseif isa(op, AbstractMultiLinearOperator) || op == (*)
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
                new_args[i] = (c == 0) ? inner_args[i] : expand_fdiff(FrechetDifferential(c, inner_args[i]), args[collect(t) .== i])
            end

            new_op(new_args...)
        end
    elseif op == (^) && isa(inner_args[2], Number)
        # chain rule, again
        base = inner_args[1]
        len = inner_args[2]
        return sum(Iterators.product(ntuple(Returns(1:len),order)...)) do t
            new_args = Vector{Any}(undef, len)
            for i = 1:len
                c = count(==(i), t)
                new_args[i] = (c == 0) ? base : expand_fdiff(FrechetDifferential(c, base), args[collect(t) .== i])
            end
            *(new_args...)
        end
    end

    return T(args...)
end

struct FrechetDerivative
    order::Int64
    FrechetDerivative(order=1) = order <= 0 ? error("Order must be a positive natural number.") : new(order)
end
(D::FrechetDerivative)(u) = FrechetDifferential(D.order, u)

Base.:*(D1::FrechetDerivative, D2::FrechetDerivative) = FrechetDerivative(D1.order+D2.order)
Base.:*(D::FrechetDerivative, DU::FrechetDifferential) = FrechetDifferential(D.order+DU.order, DU.fun)
Base.:^(D::FrechetDerivative, n::Integer) = iszero(n) ? identity : FrechetDerivative(D.order * n)

Base.nameof(D::FrechetDerivative) = :d
Base.show(io::IO, D::FrechetDerivative) = print(io, D.order == 1 ? "d" : "d^$(D.order)")


include("expand.jl")

end # module
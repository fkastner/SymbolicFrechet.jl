# SymbolicFréchet.jl
*Compute Fréchet derivatives symbolically*

This package builds on [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl) and provides the types `MultiLinearOperator` and `FrechetDifferential` to represent the corresponding concepts symbolically as well as functions to ''expand'' these operators.
At its core, it is just an implementation of the chain rule of arbitrary order for multilinear operators.

## Example

```julia-repl
julia> using SymbolicFrechet

julia> T = MultiLinearOperator(:T, 2) # provide name and number of arguments
T

julia> T(:a, :b)
T(a, b)

julia> D = FrechetDifferential(1);

julia> D(T(:a, :b))(:c) 
d(T(a,b))(c)

julia> D(T(:a, :b))(:c) |> expand_fdiffs
d(T)(a, b, c) + T(d(a)(c), b) + T(a, d(b)(c))

julia> FrechetDifferential(2)(D(:u)(:x))(:y,:z) |> expand_fdiffs
d(u)(d^2(x)(y, z)) + d^3(u)(x, y, z) + d^2(u)(d(x)(z), y) + d^2(u)(d(x)(y), z)
```

> [!NOTE]
> Here we used `Symbol`s as basic variables but in practice one could/would use symbolic variables from `SymbolicUtils.jl`. 
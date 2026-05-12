# CartesianRuns.jl

A read-only `AbstractVector{CartesianIndex{N}}` view of the positions where a
boolean mask is `true`, using
[SAMURAI](https://github.com/hpc-maths/samurai)-style interval compression for
arbitrary `N`.

## Coordinate spaces

Two coordinate spaces are distinguished throughout:

- **Mask space**: the Cartesian indices of the original array. Contains gaps —
  the `false` cells.
- **Compact space**: the 1-based linear positions `1…count(mask)`, enumerating
  only the `true` cells in column-major order. No gaps by definition.

`CartesianRunIndices` bridges the two spaces via a pair of dual operations:

| Operation   | Direction              | Julia API        |
|-------------|------------------------|------------------|
| `cri[k]`    | compact → mask         | `Base.getindex`  |
| `ci ∈ cri`  | mask → compact (membership) | `Base.in`   |

The representation itself is also dual: `CartesianRunIndices(mask)` compresses
a boolean mask into interval form, and `expand(cri)` reconstructs the mask
exactly. They are exact inverses:

```julia
expand(CartesianRunIndices(mask)) == mask   # always true
```

This duality mirrors the philosophy of
[SAMURAI](https://github.com/hpc-maths/samurai), which separates the
compressed mesh structure from the data it indexes.

## Example

```julia
using CartesianRuns

mask = Bool[0 1 0; 1 1 0; 0 1 1]
cri  = CartesianRunIndices(mask)

collect(cri) == findall(mask)      # true — compact → mask
CartesianIndex(2, 1) ∈ cri        # true  — mask membership
```

## Set operations

`CartesianRunIndices` supports set-algebraic operations that work directly on the
interval representation — no boolean mask is ever materialised.

### Unary operations

```julia
mask = Bool[0 1; 1 1; 0 1]
cri  = CartesianRunIndices(mask)

expand(cri) == mask        # true — exact inverse of constructor
complement(cri)            # positions in domain(cri) not in cri
```

`CartesianRunIndices(mask)` and `expand` are exact inverses.

### Binary operations

Binary operations require both operands to share the same `domain`.

```julia
A = CartesianRunIndices(Bool[0 1; 1 1; 0 1])
B = CartesianRunIndices(Bool[1 1; 1 0; 0 1])

intersect(A, B)  # positions in both A and B
union(A, B)      # positions in A or B (or both)
setdiff(A, B)    # positions in A but not in B
```

## API reference

```@docs
Interval
shift
CartesianRunIndices
domain
Base.size(::CartesianRunIndices)
Base.getindex(::CartesianRunIndices, ::Int)
Base.in(::CartesianIndex{N}, ::CartesianRunIndices{N}) where {N}
expand
complement
Base.intersect(::CartesianRunIndices{N}, ::CartesianRunIndices{N}) where {N}
Base.setdiff(::CartesianRunIndices{N}, ::CartesianRunIndices{N}) where {N}
Base.union(::CartesianRunIndices{N}, ::CartesianRunIndices{N}) where {N}
```

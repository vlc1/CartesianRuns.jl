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
a boolean mask into interval form, and `expand(cri, axes(mask))` reconstructs
the mask exactly. They are exact inverses:

```julia
expand(CartesianRunIndices(mask), axes(mask)) == mask   # always true
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

expand(cri, axes(mask)) == mask    # true — exact inverse of constructor
complement(cri, axes(mask))        # positions in axes(mask) not in cri
```

`CartesianRunIndices(mask)` and `expand(cri, axes(mask))` are exact inverses.

`expand` requires a **1-based domain** (`NTuple{N,Base.OneTo{Int}}`); for
non-1-based domains load `OffsetArrays.jl` first, which returns an
`OffsetArray{Bool,N}` with `axes` matching `domain`.
`complement` accepts **any `AbstractUnitRange`** and always returns a
`CartesianRunIndices` (no `OffsetArray` needed):

```julia
# The constructor and complement accept any AbstractUnitRange.
cri_sub = CartesianRunIndices(mask, (1:2, 2:3))   # restricted to sub-domain
c   = complement(cri_sub, (1:2, 2:3))  # CartesianRunIndices within sub-domain

# expand with a non-1-based domain requires OffsetArrays.
using OffsetArrays
out = expand(cri_sub, (1:2, 2:3))      # OffsetArray{Bool,2} with axes (1:2, 2:3)
```

`expand` returns an `OffsetArray{Bool,N}` whose `axes` match `domain` when
loaded with OffsetArrays.

A second constructor restricts the scan to a sub-domain:

```julia
cri_sub = CartesianRunIndices(mask, (1:2, 1:2))  # only rows 1:2, cols 1:2
```

`domain[d]` must be a subset of `axes(mask, d)` for every dimension;
an `ArgumentError` is thrown otherwise.

### Binary operations

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
Base.size(::CartesianRunIndices)
Base.getindex(::CartesianRunIndices, ::Int)
Base.in(::CartesianIndex{N}, ::CartesianRunIndices{N}) where {N}
expand
complement
Base.intersect(::CartesianRunIndices{N}, ::CartesianRunIndices{N}) where {N}
Base.setdiff(::CartesianRunIndices{N}, ::CartesianRunIndices{N}) where {N}
Base.union(::CartesianRunIndices{N}, ::CartesianRunIndices{N}) where {N}
```

# CartesianRuns.jl

A read-only `AbstractVector{CartesianIndex{N}}` view of the positions where a
boolean mask is `true`, using SAMURAI-style interval compression for arbitrary
`N`.

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

## Example

```julia
using CartesianRuns

mask = Bool[0 1 0; 1 1 0; 0 1 1]
cri  = CartesianRunIndices(mask)

collect(cri) == findall(mask)      # true — compact → mask
CartesianIndex(2, 1) ∈ cri        # true  — mask membership
```

## API reference

```@docs
Interval
shift
CartesianRunIndices
Base.size(::CartesianRunIndices)
Base.getindex(::CartesianRunIndices, ::Int)
Base.in(::CartesianIndex{N}, ::CartesianRunIndices{N}) where {N}
```

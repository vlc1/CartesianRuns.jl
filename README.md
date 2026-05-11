# CartesianRuns.jl

A Julia package providing a read-only `AbstractVector{CartesianIndex{N}}` view
of the positions where a boolean mask is `true`, using SAMURAI-style interval
compression for arbitrary `N`.

Two coordinate spaces are distinguished throughout:

- **Mask space**: the Cartesian indices of the original array (e.g.
  `CartesianIndex(i, j)` for a 2D mask). Contains gaps — the `false` cells.
- **Compact space**: the 1-based linear positions 1…`count(mask)`, enumerating
  only the `true` cells in column-major order. No gaps by definition.

`CartesianRunIndices` bridges the two spaces via a pair of dual operations:

| Operation | Direction | Julia API |
|-----------|-----------|-----------|
| `cri[k]`  | compact → mask | `Base.getindex` |
| `ci ∈ cri` | mask → compact (membership) | `Base.in` |

Both are backed by `Interval` objects — each storing a contiguous mask-space
run (`mask::UnitRange{Int}`) alongside the corresponding compact-space run
(`compact::UnitRange{Int}`), always of equal length.

## Example

```julia
using CartesianRuns

mask = Bool[0 1 0; 1 1 0; 0 1 1]
cri  = CartesianRunIndices(mask)
collect(cri) == findall(mask)           # true — compact → mask
CartesianIndex(2, 1) ∈ cri             # true  — mask → compact membership
```

## Documentation

Build locally with:

```
julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=".")); Pkg.instantiate()'
julia --project=docs/ docs/make.jl
```

Inspired by [SAMURAI](https://github.com/hpc-maths/samurai).

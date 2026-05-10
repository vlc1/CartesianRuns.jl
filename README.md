# CartesianRuns.jl

A Julia package providing a read-only `AbstractVector{CartesianIndex{N}}` view
of the positions where a boolean mask is `true`, using SAMURAI-style interval
compression.

The 1D building block (`find_runs`, `Interval`) is implemented. The N-D
extension and the exported `CartesianRunIndices{N}` type are pending.

## Example

```julia
using CartesianRuns

mask = Bool[0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0]
runs = find_runs(mask)
# 2-element Vector{Interval}:
#  Interval(5:7,  -4)
#  Interval(9:10, -5)
```

For each `run` and any `i ∈ run.range`, `i + run.shift` is the 1-based linear
index of `mask[i]` in a compact data array storing only the `true` cells, in
scan order.

Inspired by [SAMURAI](https://github.com/hpc-maths/samurai).

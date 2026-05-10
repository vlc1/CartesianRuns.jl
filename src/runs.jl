"""
    find_runs(mask::AbstractVector{Bool}, prior::Int=0) -> AbstractVector{Interval}

Find all contiguous runs of `true` values in `mask` and return them as `Interval`s
in scan order.

For each returned `run` and any `i ∈ run.range`, `i + run.shift` is the 1-based
linear position of the corresponding cell in a compact data array that stores
entries only for `true` positions, in scan order.

# Arguments

- `mask::AbstractVector{Bool}`: Boolean mask; `true` values mark cells inside a run.
- `prior::Int=0`: Cumulative count of `true` values from prior scanlines, used as
  the seed for the shift offset in higher-dimensional contexts. With `prior = 0`,
  the first `true` cell maps to compact position `1` (Julia's 1-based convention).

# Returns

A vector of `Interval`s in scan order. The container type follows
`similar(mask, Interval, 0)`, so callers using non-`Vector` mask types
(e.g. `OffsetVector`) get a matching result type.

# Example

```julia
mask = Bool[0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0]
runs = find_runs(mask)
# 2-element Vector{Interval}:
#  Interval(5:7,  -4)   # cells 5,6,7 → compact positions 1,2,3
#  Interval(9:10, -5)   # cells 9,10  → compact positions 4,5
```

# Note

This 1D scan is the building block for an N-D constructor that runs `find_runs`
per scanline along dim 1, propagating the cumulative `true` count via `prior`
and recording per-line offsets externally.
"""
function find_runs(mask::AbstractVector{Bool}, prior::Int = 0)
    runs = similar(mask, Interval, 0)
    prev = false
    start = 0
    n = 0
    for (i, val) in pairs(IndexLinear(), mask)
        if val && !prev
            start = i
        elseif !val && prev
            stop = i - 1
            push!(runs, Interval(start:stop, n - start + prior + 1))
            n += stop - start + 1
        end
        prev = val
    end
    prev && push!(runs, Interval(start:lastindex(mask), n - start + prior + 1))
    runs
end

"""
    find_runs(mask::AbstractArray{Bool,N}, line::NTuple{M,Int}, prior::Int=0)

N-D entry point (planned): find runs along dim 1 of `mask` for the scanline
identified by `line`, where `line` holds the indices in dims 2..N and must have
length `N - 1`.

Currently a stub; implementation pending.
"""
function find_runs(
    mask::AbstractArray{Bool,N}, line::NTuple{M,Int}, prior::Int = 0
) where {N,M}
    M + 1 == N || throw(ArgumentError(
        "`line` must have length ndims(mask) - 1 (got $M for an $N-D mask)"
    ))
    error("N-D find_runs not yet implemented")
end

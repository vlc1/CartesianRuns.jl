"""
    push_runs!(runs::AbstractVector{Interval}, mask::AbstractVector{Bool}, prior::Int) -> Tuple{Int,Int}

Find all contiguous runs of `true` values in `mask`, pushing them as `Interval`s
onto `runs` in scan order.

For each pushed `run` and any `i ∈ run.range`, `i + run.shift` is the 1-based
linear position of the corresponding cell in a compact data array that stores
entries only for `true` positions, in scan order.

# Arguments

- `runs::AbstractVector{Interval}`: Output buffer; new intervals are appended in scan order.
- `mask::AbstractVector{Bool}`: Boolean mask; `true` values mark cells inside a run.
- `prior::Int`: Cumulative count of `true` values from prior scanlines, used as
  the seed for the shift offset in higher-dimensional contexts. Pass `0` for a
  standalone 1D scan so that the first `true` cell maps to compact position `1`
  (Julia's 1-based convention).

# Returns

`(n, m)` where:
- `n`: total `true` values seen in this scan (to be passed as `prior` for the
  next scanline in an N-D context).
- `m`: number of runs pushed onto `runs` during this call.

# Example

```julia
mask = Bool[0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0]
runs = similar(mask, Interval, 0)
(n, m) = push_runs!(runs, mask, 0)
# n = 5, m = 2
# runs[1] = Interval(5:7,  -4)   # cells 5,6,7 → compact positions 1,2,3
# runs[2] = Interval(9:10, -5)   # cells 9,10  → compact positions 4,5
```

# Note

This 1D scan is the building block for an N-D constructor that calls `push_runs!`
per scanline along dim 1, propagating the cumulative `true` count via `prior`
and recording per-line offsets externally.
"""
function push_runs!(runs::AbstractVector{Interval}, mask::AbstractVector{Bool}, prior::Int)
    prev = false
    start = 0
    n = 0
    m = 0
    for (i, val) in pairs(IndexLinear(), mask)
        if val && !prev
            start = i
        elseif !val && prev
            stop = i - 1
            push!(runs, Interval(start:stop, n - start + prior + 1))
            n += stop - start + 1
            m += 1
        end
        prev = val
    end
    if prev
        stop = lastindex(mask)
        push!(runs, Interval(start:stop, n - start + prior + 1))
        n += stop - start + 1
        m += 1
    end
    (n, m)
end

"""
    push_runs!(runs::AbstractVector{Interval}, mask::AbstractArray{Bool,N}, line::NTuple{M,Int}, prior::Int)

N-D entry point (planned): find runs along dim 1 of `mask` for the scanline
identified by `line`, where `line` holds the indices in dims 2..N and must have
length `N - 1`.

Currently a stub; implementation pending.
"""
function push_runs!(
    runs::AbstractVector{Interval}, mask::AbstractArray{Bool,N}, line::NTuple{M,Int}, prior::Int
) where {N,M}
    M + 1 == N || throw(ArgumentError(
        "`line` must have length ndims(mask) - 1 (got $M for an $N-D mask)"
    ))
    error("N-D push_runs! not yet implemented")
end

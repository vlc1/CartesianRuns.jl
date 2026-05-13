"""
    _build_runs!(runs::AbstractVector{Interval}, mask::AbstractVector{Bool},
                 range::AbstractUnitRange{Int}, prior::Int) -> Tuple{Int,Int}

Find all contiguous runs of `true` values within `range` of `mask`, pushing
them as `Interval`s onto `runs` in scan order.

For each pushed `run` and any `i ∈ run.mask`, `i + shift(run)` is the 1-based
linear position of the corresponding cell in a compact data array that stores
entries only for `true` positions, in scan order.

# Arguments

- `runs::AbstractVector{Interval}`: Output buffer; new intervals are appended in scan order.
- `mask::AbstractVector{Bool}`: Boolean mask; `true` values mark cells inside a run.
- `range::AbstractUnitRange{Int}`: Subset of valid indices into `mask` to scan.
  Must be a sub-range of `axes(mask, 1)`. Indices in `range` are used directly
  as mask-space coordinates, so callers that pass a sub-range automatically
  produce intervals in the original coordinate space.
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
(n, m) = _build_runs!(runs, mask, axes(mask, 1), 0)
# n = 5, m = 2
# runs[1] = Interval(5:7,  -4)   # mask 5:7  → compact 1:3
# runs[2] = Interval(9:10, -5)   # mask 9:10 → compact 4:5
```

# Note

This 1D scan is the building block for an N-D constructor that calls `_build_runs!`
per scanline along dim 1, propagating the cumulative `true` count via `prior`
and recording per-line offsets externally.
"""
function _build_runs!(runs::AbstractVector{Interval}, mask::AbstractVector{Bool},
                      range::AbstractUnitRange{Int}, prior::Int)
    prev = false
    start = 0
    n = 0
    m = 0
    for i in range
        val = mask[i]
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
        stop = last(range)
        push!(runs, Interval(start:stop, n - start + prior + 1))
        n += stop - start + 1
        m += 1
    end
    (n, m)
end

# Backward-compatible 3-arg form: scans the full mask axis.
_build_runs!(runs, mask, prior::Int) = _build_runs!(runs, mask, axes(mask, 1), prior)

"""
    _build_fused!(ivs, offs, x_prior, d_prior, mask, domain) -> (x_prior, n_d, m_d)

Internal recursive kernel. Populates the pre-allocated interval and offset
buffers for all dimensions in a single pass over `mask`, without allocating
any intermediate arrays.

# Arguments

- `ivs`: `NTuple{N}` of `AbstractVector{Interval}` — one buffer per dimension.
- `offs`: `NTuple{N-1}` of `AbstractVector{Int}` — CSR offset tables; each
  `offs[d]` connects active dim-`d+1` positions to their dim-`d` intervals:
  `ivs[d][offs[d][p-1]+1 : offs[d][p]]` are the dim-`d` intervals for the
  `p`-th active dim-`d+1` position.
- `x_prior`: cumulative count of `true` dim-1 cells seen across all previously
  processed scanlines (threaded through every recursive call).
- `d_prior`: cumulative count of active dim-`D` positions from previous calls
  at this recursion level (supplied by the dim-`D+1` caller; `0` at the top level).
- `mask`: the slice of the boolean mask to process.
- `domain`: `NTuple{D}` of `AbstractUnitRange{Int}` specifying which index range
  to scan along each dimension. Must be sub-ranges of the corresponding axes of
  `mask`. Using the actual mask-space indices (not 1-based sub-array indices)
  ensures that all stored intervals are in the original coordinate space.

# Returns

`(x_prior, n_d, m_d)` where:
- `x_prior`: updated cumulative dim-1 count (to be forwarded to the next call).
- `n_d`: number of active dim-`D` positions in this slice (caller adds this to
  its own `inner_prior` before the next recursive call).
- `m_d`: number of dim-`D` intervals pushed (caller uses `m_d > 0` to decide
  whether this position is "active" at the next level up).

# Implementation note

All run-detection state (`prev`, `start`, `n_d`, `m_d`, `inner_prior`) lives
in local variables on the call stack — no heap allocation occurs beyond writes
to the pre-allocated `ivs` and `offs` buffers.
"""
function _build_fused! end

# D=1 base: domain is a NTuple{1}; d_prior is accepted (passed from D=2 caller)
# but unused; x_prior is what _build_runs! needs.  Returns (new_x_prior, n_d, m_d).
function _build_fused!(ivs, offs, x_prior::Int, _::Int,
                       mask::AbstractVector{Bool},
                       domain::NTuple{1,<:AbstractUnitRange{Int}})
    (n, m) = _build_runs!(ivs[1], mask, domain[1], x_prior)
    (x_prior + n, n, m)
end

# D≥2: all mutable state lives in local variables on the call stack — no heap
# allocation beyond the pre-allocated ivs/offs vectors.
function _build_fused!(ivs, offs, x_prior::Int, d_prior::Int,
                       mask::AbstractArray{Bool,D},
                       domain::NTuple{D,<:AbstractUnitRange{Int}}) where {D}
    prev        = false
    start       = first(domain[D])
    n_d         = 0   # active D-positions counted so far (updated at run-close)
    m_d         = 0   # D-intervals pushed so far
    inner_prior = 0   # prior for dim D-1, accumulated across slices

    for k in domain[D]
        slice = view(mask, ntuple(_ -> :, Val(D - 1))..., k)
        (x_prior, n_inner, m_inner) =
            _build_fused!(ivs, offs, x_prior, inner_prior, slice, Base.front(domain))
        inner_prior += n_inner

        active = m_inner > 0
        active && push!(offs[D - 1], offs[D - 1][end] + m_inner)

        if active && !prev
            start = k
        elseif !active && prev
            push!(ivs[D], Interval(start:(k - 1), n_d - start + d_prior + 1))
            n_d += (k - 1) - start + 1
            m_d += 1
        end
        prev = active
    end

    if prev
        stop = last(domain[D])
        push!(ivs[D], Interval(start:stop, n_d - start + d_prior + 1))
        n_d += stop - start + 1
        m_d += 1
    end
    (x_prior, n_d, m_d)
end

"""
    _build_cartesian_runs(mask::AbstractArray{Bool,N}) -> (ivs, offs)
    _build_cartesian_runs(mask::AbstractArray{Bool,N}, domain::NTuple{N}) -> (ivs, offs)

Internal constructor helper. Scans `mask` (restricted to `domain` when provided)
and returns the fully populated interval and offset buffers needed to build a
`CartesianRunIndices{N}`.

Returns `(ivs, offs)` where:
- `ivs::NTuple{N, AbstractVector{Interval}}` — dim-`d` intervals in `ivs[d]`.
- `offs::NTuple{N-1, AbstractVector{Int}}` — CSR offset tables; `offs[d]`
  maps active dim-`d+1` positions to their slice of `ivs[d]`.

The 1D method returns `((buf,), ())` (empty offsets tuple). The N-D method
pre-allocates all buffers and delegates to `_build_fused!` for a single-pass,
allocation-free fill.
"""
function _build_cartesian_runs end

# Default: use the full mask axes as domain.
_build_cartesian_runs(mask::AbstractVector{Bool}) =
    _build_cartesian_runs(mask, (axes(mask, 1),))

_build_cartesian_runs(mask::AbstractArray{Bool,N}) where {N} =
    _build_cartesian_runs(mask, axes(mask))

# Base case: 1D mask with explicit range.
function _build_cartesian_runs(mask::AbstractVector{Bool},
                               domain::NTuple{1,<:AbstractUnitRange{Int}})
    buf = similar(mask, Interval, 0)
    _build_runs!(buf, mask, domain[1], 0)
    ((buf,), ())
end

# N-D: pre-allocate all buffers, then fill them in one fused pass.
function _build_cartesian_runs(mask::AbstractArray{Bool,N},
                               domain::NTuple{N,<:AbstractUnitRange{Int}}) where {N}
    ivs  = ntuple(_ -> similar(mask, Interval, 0), Val(N))
    offs = ntuple(Val(N - 1)) do _
        o = similar(mask, Int, 1); o[begin] = 1; o
    end
    _build_fused!(ivs, offs, 0, 0, mask, domain)
    (ivs, offs)
end

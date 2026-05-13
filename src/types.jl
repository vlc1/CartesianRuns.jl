"""
    Interval(mask, shift)

A run of consecutive `true` cells in a boolean mask, paired with the
corresponding run of compact positions.

# Fields

- `mask::UnitRange{Int}`: Indices of the run in the original mask (mask space).
- `compact::UnitRange{Int}`: Corresponding 1-based positions in the compact
  array (compact space). Always satisfies `length(compact) == length(mask)`.

The two ranges always have the same length; the constructor enforces this by
deriving `compact` from `mask` and `shift`:

    Interval(mask::UnitRange{Int}, shift::Int)

# Methods

    shift(iv::Interval) -> Int

Retrieve the integer offset such that `i + shift(iv)` maps any `i ∈ iv.mask`
to its compact position. Equivalent to `iv.compact.start - iv.mask.start`.
"""
struct Interval
    mask::UnitRange{Int}
    compact::UnitRange{Int}
    Interval(mask::UnitRange{Int}, shift::Int) = new(mask, mask .+ shift)
end

"""
    shift(iv::Interval) -> Int

Return the integer offset between mask space and compact space for `iv`:
for any `i ∈ iv.mask`, `i + shift(iv)` is the corresponding compact position.
"""
shift(iv::Interval) = iv.compact.start - iv.mask.start

"""
    CartesianRunIndices{N,M,I,O} <: AbstractVector{CartesianIndex{N}}

A read-only `AbstractVector` view of the positions where an N-dimensional
boolean mask is `true`, using [SAMURAI](https://github.com/hpc-maths/samurai)-style interval compression.

# Type parameters

- `N`: number of spatial dimensions.
- `M = N - 1`: number of offset arrays (enforced by inner constructor).
- `I <: AbstractVector{Interval}`: element type of each dimension's interval vector.
- `O <: AbstractVector{Int}`: element type of each CSR offset array.

# Fields

- `intervals::NTuple{N,I}`: packed interval vectors, one per dimension.
  `intervals[d]` contains all dim-d runs in scan order across all scanlines.
- `offsets::NTuple{M,O}`: 1-based CSR offset arrays connecting active cells
  in dimension `d+1` to their ranges in `intervals[d]`.
  Empty tuple for `N = 1`.

# Construction

    CartesianRunIndices(mask::AbstractArray{Bool,N}) -> CartesianRunIndices

Construct from a boolean mask of any dimension.

    CartesianRunIndices(mask::AbstractArray{Bool,N}, domain::NTuple{N,AbstractUnitRange{Int}}) -> CartesianRunIndices

Construct from a boolean mask, restricting the scan to positions within
`domain`. Each element of `domain` must be a subset of the corresponding
axis of `mask`; an `ArgumentError` is thrown otherwise. The returned
`CartesianRunIndices` contains only the `true` cells of `mask` that lie
within `domain`.
"""
struct CartesianRunIndices{N,M,I<:AbstractVector{Interval},O<:AbstractVector{Int}} <:
        AbstractVector{CartesianIndex{N}}
    intervals::NTuple{N,I}
    offsets::NTuple{M,O}

    # General inner constructor (N ≥ 2, M ≥ 1)
    function CartesianRunIndices(
        intervals::NTuple{N,I}, offsets::NTuple{M,O}
    ) where {N,M,I<:AbstractVector{Interval},O<:AbstractVector{Int}}
        M + 1 == N || throw(ArgumentError(
            "length(offsets) must equal ndims - 1: got $M offsets for N=$N"
        ))
        new{N,M,I,O}(intervals, offsets)
    end

    # Special inner constructor for 1D (empty offsets; O cannot be inferred)
    function CartesianRunIndices(
        intervals::NTuple{1,I}, ::Tuple{}
    ) where {I<:AbstractVector{Interval}}
        new{1,0,I,Vector{Int}}(intervals, ())
    end
end

# --- AbstractVector interface ---

"""
    size(cri::CartesianRunIndices) -> Tuple{Int}

Return `(n,)` where `n = count(mask)` is the total number of `true` cells.
Consistent with the `AbstractVector` interface.
"""
Base.size(cri::CartesianRunIndices) =
    (sum(length(iv.mask) for iv in cri.intervals[1]; init = 0),)

"""
    _search_compact(intervals, k) -> p

Return the index `p` of the unique interval in `intervals` whose compact range
contains `k`. The compact range of interval `p` is `intervals[p].compact`, so
the search key is `intervals[p].compact.stop`.
Intervals must be sorted and non-overlapping (guaranteed by construction).
Uses `searchsortedfirst` with a custom `lt`; `by` cannot be used as Julia
applies it to both array elements and the key.
"""
function _search_compact(intervals::AbstractVector{Interval}, k::Int)
    searchsortedfirst(intervals, k; lt = (iv, key) -> iv.compact.stop < key)
end

"""
    _recover(intervals, offsets, k) -> CartesianIndex

Reconstruct the `CartesianIndex` corresponding to compact position `k`.

Works bottom-up through the dimension hierarchy:

1. **Dim-1** (`intervals[1]`): binary-search via `_search_compact` to find
   interval `p`; the mask index is `k - intervals[1][p].shift`.
2. **Dim-2+**: `offsets[1]` is a CSR table mapping x-interval index `p` to
   the active dim-2 position `q` via `searchsortedfirst(offsets[1], p+1)-1`.
   Recurse with `(Base.tail(intervals), Base.tail(offsets), q)` to recover
   the remaining coordinates, then prepend the dim-1 index.

The 1D base case dispatches on `NTuple{1,...}` / `Tuple{}` and terminates the
recursion.
"""
function _recover end

# 1D base case
function _recover(
    intervals::NTuple{1,<:AbstractVector{Interval}},
    ::Tuple{},
    k::Int,
)
    p = _search_compact(intervals[1], k)
    CartesianIndex(intervals[1][p].mask.start + (k - intervals[1][p].compact.start))
end

# N-D recursive case
function _recover(
    intervals::NTuple{N,<:AbstractVector{Interval}},
    offsets::NTuple{M,<:AbstractVector{Int}},
    k::Int,
) where {N,M}
    p     = _search_compact(intervals[1], k)
    q     = searchsortedfirst(offsets[1], p + 1) - 1
    outer = _recover(Base.tail(intervals), Base.tail(offsets), q)
    CartesianIndex(intervals[1][p].mask.start + (k - intervals[1][p].compact.start), Tuple(outer)...)
end

"""
    getindex(cri::CartesianRunIndices, k::Int) -> CartesianIndex

Return the `k`-th `true` cell of the mask, in column-major (Julia-default)
scan order.  Delegates to `_recover` after a bounds check.
"""
function Base.getindex(cri::CartesianRunIndices, k::Int)
    @boundscheck 1 <= k <= length(cri) || throw(BoundsError(cri, k))
    _recover(cri.intervals, cri.offsets, k)
end

# Check if coord falls in any interval in iv[lo:hi].
@inline function _in_intervals(iv::AbstractVector{Interval}, lo::Int, hi::Int, coord::Int)
    r = lo - 1 + searchsortedfirst(view(iv, lo:hi), coord; lt = (a, b) -> a.mask.stop < b)
    r <= hi && iv[r].mask.start <= coord
end

# Base.in helpers — top-down through dimensions (outermost first).

# 1D base case: check dim-1 coord in iv[lo:hi].
function _in_cri(
    idx::CartesianIndex{1},
    intervals::NTuple{1,<:AbstractVector{Interval}},
    ::Tuple{},
    lo::Int, hi::Int,
)
    _in_intervals(intervals[1], lo, hi, idx[1])
end

# N-D recursive case: check outermost coord, narrow via offsets, recurse.
function _in_cri(
    idx::CartesianIndex{N},
    intervals::NTuple{N,<:AbstractVector{Interval}},
    offsets::NTuple{M,<:AbstractVector{Int}},
    lo::Int, hi::Int,
) where {N,M}
    iv    = last(intervals)
    coord = idx[N]
    r = lo - 1 + searchsortedfirst(view(iv, lo:hi), coord; lt = (a, b) -> a.mask.stop < b)
    r > hi && return false
    iv[r].mask.start <= coord || return false
    q      = iv[r].compact.start + (coord - iv[r].mask.start)
    new_lo = last(offsets)[q]
    new_hi = last(offsets)[q + 1] - 1
    _in_cri(CartesianIndex(Base.front(Tuple(idx))),
            Base.front(intervals), Base.front(offsets), new_lo, new_hi)
end

"""
    in(idx::CartesianIndex{N}, cri::CartesianRunIndices{N}) -> Bool

Return `true` if `idx` corresponds to a `true` cell in the mask used to
construct `cri`. Equivalent to `mask[idx]`, but works directly on the
compressed representation in O(N log R) time (R = maximum run count per
dimension) without reconstructing the mask.
"""
function Base.in(idx::CartesianIndex{N}, cri::CartesianRunIndices{N}) where {N}
    isempty(cri) && return false
    _in_cri(idx, cri.intervals, cri.offsets, 1, length(last(cri.intervals)))
end

# --- Outer constructors ---

function CartesianRunIndices(mask::AbstractVector{Bool})
    buf = similar(mask, Interval, 0)
    _build_runs!(buf, mask, 0)
    CartesianRunIndices((buf,), ())
end

function CartesianRunIndices(mask::AbstractArray{Bool,N}) where {N}
    intervals, offsets = _build_cartesian_runs(mask)
    CartesianRunIndices(intervals, offsets)
end

function CartesianRunIndices(
    mask::AbstractArray{Bool,N},
    domain::NTuple{N,<:AbstractUnitRange{Int}},
) where {N}
    for d in 1:N
        ax = axes(mask, d)
        dm = domain[d]
        (first(dm) >= first(ax) && last(dm) <= last(ax)) || throw(ArgumentError(
            "domain[$d] = $dm is not a subset of axes(mask, $d) = $ax"))
    end
    dom = map(r -> Int(first(r)):Int(last(r)), domain)
    intervals, offsets = _build_cartesian_runs(mask, dom)
    CartesianRunIndices(intervals, offsets)
end

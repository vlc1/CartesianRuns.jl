function _expand_runs!(
    out::AbstractVector{Bool},
    ivs::AbstractVector{Interval}, lo::Int, hi::Int,
    x_offset::Int,
)
    for k in lo:hi
        out[ivs[k].mask .- x_offset] .= true
    end
end

function _expand_fused!(
    out::AbstractVector{Bool},
    ivs::NTuple{1,<:AbstractVector{Interval}}, ::Tuple{},
    lo::Int, hi::Int,
    domain::NTuple{1,<:AbstractUnitRange{Int}},
)
    _expand_runs!(out, ivs[1], lo, hi, first(domain[1]) - 1)
end

function _expand_fused!(
    out,
    ivs::NTuple{N,<:AbstractVector{Interval}},
    offs::NTuple{M,<:AbstractVector{Int}},
    lo::Int, hi::Int,
    domain::NTuple{N,<:AbstractUnitRange{Int}},
) where {N,M}
    a_outer    = last(ivs)
    a_off      = last(offs)
    front_ivs  = Base.front(ivs)
    front_offs = Base.front(offs)
    outer_offset = first(last(domain)) - 1
    for k in lo:hi
        iv = a_outer[k]
        for r in iv.mask
            ia_lo, ia_hi = _inner_slice(iv, a_off, r)
            _expand_fused!(
                view(out, ntuple(_ -> :, Val(N - 1))..., r - outer_offset),
                front_ivs, front_offs, ia_lo, ia_hi, Base.front(domain))
        end
    end
end

"""
    expand(cri::CartesianRunIndices{N}, domain::NTuple{N,Base.OneTo{Int}}) -> Array{Bool,N}

Reconstruct the boolean mask from `cri`: the returned array has `axes` equal to
`domain` and is `true` exactly at the positions in `cri`.

`CartesianRunIndices(mask)` and `expand(cri, axes(mask))` are exact inverses:

```julia
expand(CartesianRunIndices(mask), axes(mask)) == mask   # always true
```

`domain` must be a `NTuple` of `Base.OneTo{Int}` ranges (i.e., 1-based).
For non-1-based domains, load `OffsetArrays.jl` first; `expand` will then
return an `OffsetArray` with `axes` equal to `domain`.
"""
function expand(cri::CartesianRunIndices{N},
                domain::NTuple{N,Base.OneTo{Int}}) where {N}
    out = similar(cri.intervals[1], Bool, length.(domain))
    fill!(out, false)
    _expand_fused!(out, cri.intervals, cri.offsets, 1, length(last(cri.intervals)), domain)
    out
end

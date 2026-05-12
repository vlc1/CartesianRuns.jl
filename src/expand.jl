function _expand_runs!(
    out::AbstractVector{Bool},
    ivs::AbstractVector{Interval}, lo::Int, hi::Int,
)
    for k in lo:hi
        out[ivs[k].mask] .= true
    end
end

function _expand_fused!(
    out::AbstractVector{Bool},
    ivs::NTuple{1,<:AbstractVector{Interval}}, ::Tuple{},
    lo::Int, hi::Int,
)
    _expand_runs!(out, ivs[1], lo, hi)
end

function _expand_fused!(
    out,
    ivs::NTuple{N,<:AbstractVector{Interval}},
    offs::NTuple{M,<:AbstractVector{Int}},
    lo::Int, hi::Int,
) where {N,M}
    a_outer    = last(ivs)
    a_off      = last(offs)
    front_ivs  = Base.front(ivs)
    front_offs = Base.front(offs)
    for k in lo:hi
        iv = a_outer[k]
        for r in iv.mask
            ia_lo, ia_hi = _inner_slice(iv, a_off, r)
            _expand_fused!(
                view(out, ntuple(_ -> :, Val(N - 1))..., r),
                front_ivs, front_offs, ia_lo, ia_hi)
        end
    end
end

"""
    expand(cri::CartesianRunIndices{N}) -> AbstractArray{Bool,N}

Reconstruct the boolean mask from `cri`: the returned array has axes
`domain(cri)` and is `true` exactly at the positions in `cri`.

`CartesianRunIndices(mask)` and `expand` are exact inverses.
"""
function expand(cri::CartesianRunIndices{N}) where {N}
    out = similar(cri.intervals[1], Bool, length.(cri.domain))
    fill!(out, false)
    _expand_fused!(out, cri.intervals, cri.offsets, 1, length(last(cri.intervals)))
    out
end

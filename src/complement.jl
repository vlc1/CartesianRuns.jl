function _complement_runs!(
    out::AbstractVector{Interval},
    a_ivs::AbstractVector{Interval}, a_lo::Int, a_hi::Int,
    x_range::UnitRange{Int}, prior::Int,
)
    pos = x_range.start
    n = 0; m = 0
    for k in a_lo:a_hi
        iv = a_ivs[k]
        if pos < iv.mask.start
            push!(out, Interval(pos:(iv.mask.start - 1), n - pos + prior + 1))
            n += iv.mask.start - pos; m += 1
        end
        pos = iv.mask.stop + 1
    end
    if pos <= x_range.stop
        push!(out, Interval(pos:x_range.stop, n - pos + prior + 1))
        n += x_range.stop - pos + 1; m += 1
    end
    (n, m)
end

function _complement_fused!(
    out_ivs, ::Tuple{},
    a_ivs::NTuple{1,<:AbstractVector{Interval}}, ::Tuple{},
    a_lo::Int, a_hi::Int,
    x_prior::Int, _::Int, domain,
)
    (n, m) = _complement_runs!(out_ivs[1], a_ivs[1], a_lo, a_hi, domain[1], x_prior)
    (x_prior + n, n, m)
end

function _complement_fused!(
    out_ivs, out_offs,
    a_ivs::NTuple{N,<:AbstractVector{Interval}},
    a_offs::NTuple{M,<:AbstractVector{Int}},
    a_lo::Int, a_hi::Int,
    x_prior::Int, d_prior::Int, domain,
) where {N,M}
    front_out  = Base.front(out_ivs);  front_ooffs = Base.front(out_offs)
    front_a    = Base.front(a_ivs);    front_ao    = Base.front(a_offs)
    a_outer    = last(a_ivs)
    a_off      = last(a_offs)
    inner_dom  = Base.front(domain)

    inner_prior = 0; n_d = 0; m_d = 0
    prev = false; start = 0
    i = a_lo

    for r in last(domain)
        while i <= a_hi && a_outer[i].mask.stop < r
            i += 1
        end
        ia_lo = a_hi + 1; ia_hi = a_hi
        if i <= a_hi && a_outer[i].mask.start <= r
            ia_lo, ia_hi = _inner_slice(a_outer[i], a_off, r)
        end

        (x_prior, n_inner, m_inner) = _complement_fused!(
            front_out, front_ooffs,
            front_a, front_ao, ia_lo, ia_hi,
            x_prior, inner_prior, inner_dom)
        inner_prior += n_inner

        active = m_inner > 0
        if active
            push!(last(out_offs), last(out_offs)[end] + m_inner)
            !prev && (start = r)
        elseif prev
            push!(last(out_ivs), Interval(start:(r - 1), n_d - start + d_prior + 1))
            n_d += r - 1 - start + 1; m_d += 1
        end
        prev = active
    end

    if prev
        stop = last(last(domain))
        push!(last(out_ivs), Interval(start:stop, n_d - start + d_prior + 1))
        n_d += stop - start + 1; m_d += 1
    end
    (x_prior, n_d, m_d)
end

"""
    complement(cri::CartesianRunIndices{N}) -> CartesianRunIndices{N}

Return a `CartesianRunIndices` containing every position in `domain(cri)` that
is **not** in `cri`, computed directly from the interval representation.
"""
function complement(cri::CartesianRunIndices{1})
    out = similar(cri.intervals[1], Interval, 0)
    _complement_runs!(out, cri.intervals[1], 1, length(cri.intervals[1]),
                      cri.domain[1], 0)
    CartesianRunIndices((out,), (), cri.domain)
end

function complement(cri::CartesianRunIndices{N}) where {N}
    out_ivs  = ntuple(_ -> similar(cri.intervals[1], Interval, 0), Val(N))
    out_offs = ntuple(Val(N - 1)) do _
        o = similar(cri.intervals[1], Int, 1); o[begin] = 1; o
    end
    _complement_fused!(out_ivs, out_offs,
                       cri.intervals, cri.offsets, 1, length(last(cri.intervals)),
                       0, 0, cri.domain)
    CartesianRunIndices(out_ivs, out_offs, cri.domain)
end

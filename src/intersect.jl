function _intersect_runs!(
    out::AbstractVector{Interval},
    a_ivs::AbstractVector{Interval}, a_lo::Int, a_hi::Int,
    b_ivs::AbstractVector{Interval}, b_lo::Int, b_hi::Int,
    prior::Int,
)
    n = 0; m = 0
    i, j = a_lo, b_lo
    while i <= a_hi && j <= b_hi
        a_iv = a_ivs[i]; b_iv = b_ivs[j]
        lo = max(a_iv.mask.start, b_iv.mask.start)
        hi = min(a_iv.mask.stop,  b_iv.mask.stop)
        if lo <= hi
            push!(out, Interval(lo:hi, n - lo + prior + 1))
            n += hi - lo + 1; m += 1
        end
        a_iv.mask.stop <= b_iv.mask.stop ? (i += 1) : (j += 1)
    end
    (n, m)
end

function _intersect_fused!(
    out_ivs, ::Tuple{},
    a_ivs::NTuple{1,<:AbstractVector{Interval}}, ::Tuple{}, a_lo::Int, a_hi::Int,
    b_ivs::NTuple{1,<:AbstractVector{Interval}}, ::Tuple{}, b_lo::Int, b_hi::Int,
    x_prior::Int, _::Int,
)
    (n, m) = _intersect_runs!(out_ivs[1], a_ivs[1], a_lo, a_hi,
                                           b_ivs[1], b_lo, b_hi, x_prior)
    (x_prior + n, n, m)
end

function _intersect_fused!(
    out_ivs, out_offs,
    a_ivs::NTuple{N,<:AbstractVector{Interval}},
    a_offs::NTuple{M,<:AbstractVector{Int}},
    a_lo::Int, a_hi::Int,
    b_ivs::NTuple{N,<:AbstractVector{Interval}},
    b_offs::NTuple{M,<:AbstractVector{Int}},
    b_lo::Int, b_hi::Int,
    x_prior::Int, d_prior::Int,
) where {N,M}
    front_out  = Base.front(out_ivs);  front_ooffs = Base.front(out_offs)
    front_a    = Base.front(a_ivs);    front_ao    = Base.front(a_offs)
    front_b    = Base.front(b_ivs);    front_bo    = Base.front(b_offs)
    a_outer    = last(a_ivs);          b_outer     = last(b_ivs)
    a_off      = last(a_offs);         b_off       = last(b_offs)

    inner_prior = 0; n_d = 0; m_d = 0
    prev = false; start = 0; last_r = 0
    i, j = a_lo, b_lo
    ar = i <= a_hi ? a_outer[i].mask.start : typemax(Int)
    br = j <= b_hi ? b_outer[j].mask.start : typemax(Int)

    while ar < typemax(Int) || br < typemax(Int)
        r = min(ar, br)
        if prev && r > last_r + 1
            push!(last(out_ivs), Interval(start:last_r, n_d - start + d_prior + 1))
            n_d += last_r - start + 1; m_d += 1; prev = false
        end
        ia_lo = a_hi + 1; ia_hi = a_hi
        ib_lo = b_hi + 1; ib_hi = b_hi

        if ar == r
            ia_lo, ia_hi = _inner_slice(a_outer[i], a_off, r)
            r == a_outer[i].mask.stop ? (i += 1; ar = i <= a_hi ? a_outer[i].mask.start : typemax(Int)) : (ar = r + 1)
        end
        if br == r
            ib_lo, ib_hi = _inner_slice(b_outer[j], b_off, r)
            r == b_outer[j].mask.stop ? (j += 1; br = j <= b_hi ? b_outer[j].mask.start : typemax(Int)) : (br = r + 1)
        end

        if ia_lo <= ia_hi && ib_lo <= ib_hi
            (x_prior, n_inner, m_inner) = _intersect_fused!(
                front_out, front_ooffs,
                front_a, front_ao, ia_lo, ia_hi,
                front_b, front_bo, ib_lo, ib_hi,
                x_prior, inner_prior)
            inner_prior += n_inner
            active = m_inner > 0
            if active
                push!(last(out_offs), last(out_offs)[end] + m_inner)
                !prev && (start = r)
            elseif prev
                push!(last(out_ivs), Interval(start:(r - 1), n_d - start + d_prior + 1))
                n_d += r - 1 - start + 1; m_d += 1
            end
            prev = active; last_r = r
        else
            if prev
                push!(last(out_ivs), Interval(start:(r - 1), n_d - start + d_prior + 1))
                n_d += r - 1 - start + 1; m_d += 1; prev = false
            end
        end
    end
    if prev
        push!(last(out_ivs), Interval(start:last_r, n_d - start + d_prior + 1))
        n_d += last_r - start + 1; m_d += 1
    end
    (x_prior, n_d, m_d)
end

function Base.intersect(a::CartesianRunIndices{1}, b::CartesianRunIndices{1})
    _check_domain(a, b)
    out = similar(a.intervals[1], Interval, 0)
    _intersect_runs!(out, a.intervals[1], 1, length(a.intervals[1]),
                          b.intervals[1], 1, length(b.intervals[1]), 0)
    CartesianRunIndices((out,), (), a.domain)
end

"""
    intersect(a::CartesianRunIndices{N}, b::CartesianRunIndices{N}) -> CartesianRunIndices{N}

Return a `CartesianRunIndices` containing every `CartesianIndex` that is `true`
in **both** `a` and `b`.  Both operands must have the same `domain`.
"""
function Base.intersect(a::CartesianRunIndices{N}, b::CartesianRunIndices{N}) where {N}
    _check_domain(a, b)
    out_ivs  = ntuple(_ -> similar(a.intervals[1], Interval, 0), Val(N))
    out_offs = ntuple(Val(N - 1)) do _
        o = similar(a.intervals[1], Int, 1); o[begin] = 1; o
    end
    _intersect_fused!(out_ivs, out_offs,
                      a.intervals, a.offsets, 1, length(last(a.intervals)),
                      b.intervals, b.offsets, 1, length(last(b.intervals)),
                      0, 0)
    CartesianRunIndices(out_ivs, out_offs, a.domain)
end

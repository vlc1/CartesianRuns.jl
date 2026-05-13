module OffsetArraysExt

using CartesianRuns, OffsetArrays

import CartesianRuns: CartesianRunIndices, _expand_fused!

"""
    expand(cri::CartesianRunIndices{N}, domain::NTuple{N,<:AbstractUnitRange{Int}}) -> OffsetArray{Bool,N}

Reconstruct the boolean mask from `cri` over `domain`. When `domain` is not
1-based the returned array is an `OffsetArray` with `axes(out) == domain`, so
that `out[ci]` is `true` exactly when `ci ∈ cri`.

Requires `using OffsetArrays`.
"""
function CartesianRuns.expand(
        cri::CartesianRunIndices{N},
        domain::NTuple{N,<:AbstractUnitRange{Int}}) where {N}
    oa = OffsetArray(falses(length.(domain)...), map(r -> first(r) - 1, domain)...)
    # The kernel (_expand_fused!) subtracts `first(domain[d]) - 1` from every
    # mask-space coordinate before indexing `out`. For a 1-based domain that
    # offset is 0 and the subtraction is a no-op. For a non-1-based domain
    # (e.g. 3:7) the kernel would subtract 2 — correct for a plain Array, but
    # wrong for `oa`: OffsetArray already maps index 5 → storage[3], so a
    # further -2 would land on storage[1] instead.
    #
    # Fix: pass a *proxy* domain whose starts are all 1. Every offset becomes
    # 0, so the kernel writes `oa[mask_pos]` unchanged, and OffsetArray's own
    # indexing performs the sole coordinate translation. Passing the real
    # domain would double-subtract the offset and corrupt the output.
    proxy = map(r -> Base.OneTo(length(r)), domain)
    _expand_fused!(oa, cri.intervals, cri.offsets,
                   1, length(last(cri.intervals)), proxy)
    oa
end

end # module

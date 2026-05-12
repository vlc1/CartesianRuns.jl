# Shared helpers used by all set operations.

function _check_domain(a::CartesianRunIndices{N}, b::CartesianRunIndices{N}) where {N}
    a.domain == b.domain || throw(ArgumentError(
        "domain mismatch: $(a.domain) vs $(b.domain)"))
end

# Given outer Interval `iv`, its CSR offset table `offs`, and row `r`,
# return the [lo, hi] slice into the inner interval vector.
@inline function _inner_slice(iv::Interval, offs::AbstractVector{Int}, r::Int)
    c = iv.compact.start + (r - iv.mask.start)
    (offs[c], offs[c + 1] - 1)
end

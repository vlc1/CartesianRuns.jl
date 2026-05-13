This new branch (ghost) will add of a new constructor, that takes two arguments:
```julia
function CCrtesianRunIndices(mask::AbstractArray{Bool,N}, domain::NTuple{N,AbstractUnitRange{Int}}) where {N}
```
Element-wise, `domain` is a subset of of `axes(mask)` (the constructor checks this on construction).

Then, the construction proceeds by only adding the indices where `mask` is `true` that lie within ranges provided in `domain`.

Propose a plan detailing this implementation.

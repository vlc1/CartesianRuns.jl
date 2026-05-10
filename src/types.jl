"""
    Interval(range, shift)

A run of consecutive `true` cells in a boolean mask, paired with an offset that
maps each cell index in `range` to its position in a compact data array.

# Fields

- `range::UnitRange{Int}`: The indices of the run in the original mask.
- `shift::Int`: For any `i ∈ range`, `i + shift` is the 1-based linear position
  of the corresponding entry in the compact data array.
"""
struct Interval
    range::UnitRange{Int}
    shift::Int
end

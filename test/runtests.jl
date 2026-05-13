using CartesianRuns
using Test

_to_set(cri) = Set(collect(cri))

# ──────────────────────────────────────────────────────────────────────────────
# CartesianRunIndices — construction and indexing
# ──────────────────────────────────────────────────────────────────────────────

@testset "Interval and shift" begin
    iv = Interval(3:5, 2)        # mask = 3:5, compact = 5:7
    @test iv.mask    == 3:5
    @test iv.compact == 5:7
    @test shift(iv) == 2

    iv2 = Interval(1:4, -3)      # compact = -2:1
    @test shift(iv2) == -3
end

@testset "construction: 1D" begin
    mask = Bool[0, 1, 1, 0, 1]
    cri  = CartesianRunIndices(mask)
    @test size(cri) == (count(mask),)
    @test length(cri.intervals[1]) == 2
    @test cri.intervals[1][1].mask == 2:3
    @test cri.intervals[1][2].mask == 5:5
    @test cri.offsets === ()
    @test collect(cri) == CartesianIndex.(findall(mask))
end

@testset "construction: 2D" begin
    mask = Bool[1 0 0 1; 1 1 0 1; 0 1 0 1]
    cri  = CartesianRunIndices(mask)
    @test length(cri) == count(mask)
    @test length(cri.intervals[1]) == 3
    @test length(cri.intervals[2]) == 2
    @test length(cri.offsets[1])   == 4
    @test collect(cri) == findall(mask)
end

@testset "construction: 3D" begin
    mask = rand(Bool, 4, 5, 3)
    @test collect(CartesianRunIndices(mask)) == findall(mask)
end

@testset "construction: edge cases" begin
    @test length(CartesianRunIndices(trues(3))) == 3
    @test isempty(CartesianRunIndices(falses(3)))
    @test isempty(CartesianRunIndices(falses(3, 4)))
    m = falses(4, 4); m[2, 3] = true
    @test CartesianRunIndices(m)[1] == CartesianIndex(2, 3)
end

@testset "construction: domain" begin
    # 1D: restrict to a subrange
    mask1 = Bool[0, 1, 1, 0, 1, 0]
    dom1  = (2:5,)
    cri1  = CartesianRunIndices(mask1, dom1)
    @test Set(collect(cri1)) == Set(CartesianIndex.(findall(view(mask1, 2:5)) .+ 1))
    @test collect(cri1) == CartesianIndex.(findall(@view mask1[2:5]) .+ 1)

    # 2D: restrict to a sub-rectangle
    mask2 = Bool[1 0 1; 0 1 0; 1 0 1]
    dom2  = (1:2, 1:2)
    cri2  = CartesianRunIndices(mask2, dom2)
    sub2  = view(mask2, 1:2, 1:2)
    @test collect(cri2) == findall(sub2)

    # Regression: domain starting at index > 1 must return original-space indices.
    # Previously (OffsetArray-based impl) pairs(IndexLinear(), view) returned
    # 1-based indices instead of the original mask coordinates.
    mask4 = Bool[0, 0, 1, 1, 0, 1]
    dom4  = (3:6,)
    cri4  = CartesianRunIndices(mask4, dom4)
    @test collect(cri4) == CartesianIndex.([3, 4, 6])

    # 2D case: neither dimension starts at 1
    mask5 = Bool[0 0 0; 0 1 1; 0 1 0]
    dom5  = (2:3, 2:3)
    cri5  = CartesianRunIndices(mask5, dom5)
    # All returned indices must lie within dom5 and be true in mask5
    @test all(ci -> mask5[ci] && ci[1] ∈ dom5[1] && ci[2] ∈ dom5[2], collect(cri5))
    @test length(cri5) == count(@view mask5[2:3, 2:3])


    mask3 = rand(Bool, 4, 5)
    @test collect(CartesianRunIndices(mask3, axes(mask3))) ==
          collect(CartesianRunIndices(mask3))

    # out-of-bounds domain throws
    @test_throws ArgumentError CartesianRunIndices(mask2, (0:2, 1:2))
    @test_throws ArgumentError CartesianRunIndices(mask2, (1:2, 1:4))
end

@testset "Base.in" begin
    for mask in (Bool[0, 1, 1, 0, 1, 0, 0, 1, 1, 1],
                 Bool[1 0 1; 0 1 0; 1 1 0; 0 0 1],
                 rand(Bool, 4, 5, 3))
        cri = CartesianRunIndices(mask)
        @test all(idx -> (idx ∈ cri) == mask[idx], CartesianIndices(mask))
    end
    @test CartesianIndex(1) ∉ CartesianRunIndices(Bool[])
    @test CartesianIndex(1, 1) ∉ CartesianRunIndices(falses(3, 3))
end

# ──────────────────────────────────────────────────────────────────────────────
# expand
# ──────────────────────────────────────────────────────────────────────────────

@testset "expand" begin
    for mask in (Bool[0, 1, 1, 0, 1, 0],
                 Bool[1 0 1; 0 1 0; 1 0 1],
                 rand(Bool, 4, 5, 3),
                 falses(3, 3), trues(3, 3))
        @test expand(CartesianRunIndices(mask), axes(mask)) == mask
    end

    # domain-restricted expand: output is length.(domain), not full mask size
    mask6 = Bool[0, 1, 1, 0, 1, 0]
    dom6  = (2:5,)
    cri6  = CartesianRunIndices(mask6, dom6)
    out6  = expand(cri6, dom6)
    @test size(out6) == (length(dom6[1]),)
    @test out6 == mask6[2:5]                    # sub-mask in 1-based output

    # 2D domain-restricted expand
    mask7 = Bool[1 0 1; 0 1 0; 1 0 1]
    dom7  = (1:2, 2:3)
    cri7  = CartesianRunIndices(mask7, dom7)
    @test expand(cri7, dom7) == mask7[1:2, 2:3]
end

@testset "getindex" begin
    mask = Bool[0, 1, 1, 0, 1]
    cri  = CartesianRunIndices(mask)
    trues_pos = CartesianIndex.(findall(mask))
    for k in 1:length(cri)
        @test cri[k] == trues_pos[k]
    end
    # 2D
    mask2 = Bool[1 0 1; 0 1 0; 1 0 1]
    cri2  = CartesianRunIndices(mask2)
    @test [cri2[k] for k in 1:length(cri2)] == findall(mask2)
end

# ──────────────────────────────────────────────────────────────────────────────
# Set operations — shared infrastructure
# ──────────────────────────────────────────────────────────────────────────────

# ──────────────────────────────────────────────────────────────────────────────
# complement
# ──────────────────────────────────────────────────────────────────────────────

@testset "complement: $label" for (label, mask) in (
    ("1D", Bool[0, 1, 1, 0, 1, 0]),
    ("2D", Bool[1 0 1; 0 1 0; 1 0 1]),
    ("3D", rand(Bool, 4, 5, 3)),
)
    cri = CartesianRunIndices(mask)
    c   = complement(cri, axes(mask))
    @test _to_set(c) == _to_set(CartesianRunIndices(.!mask))
    @test length(cri) + length(c) == prod(size(mask))
    @test _to_set(complement(c, axes(mask))) == _to_set(cri)   # involution
end

@testset "complement: trivial" begin
    @test _to_set(complement(CartesianRunIndices(falses(5)), (1:5,))) ==
          _to_set(CartesianRunIndices(trues(5)))
    @test isempty(complement(CartesianRunIndices(trues(5)), (1:5,)))
end

@testset "complement: restricted domain" begin
    # complement over a subdomain — only positions *within the domain* not in cri
    mask = Bool[0, 1, 1, 0, 1, 0]
    dom  = (2:5,)
    cri  = CartesianRunIndices(mask, dom)           # has {2,3,5}
    c    = complement(cri, dom)                     # {4} within 2:5
    @test _to_set(c) == Set([CartesianIndex(4)])
    @test length(cri) + length(c) == length(dom[1])
end

# ──────────────────────────────────────────────────────────────────────────────
# Binary set operations — intersect / setdiff / union
# ──────────────────────────────────────────────────────────────────────────────

@testset "binary ops: $label" for (label, masks) in (
    ("1D", (Bool[0, 1, 1, 0, 1, 1, 0, 1], Bool[1, 1, 0, 0, 1, 0, 1, 1])),
    ("2D", (Bool[1 0 1; 1 1 1; 0 1 0],    Bool[0 1 1; 1 1 0; 0 0 1])),
    ("3D", (rand(Bool, 4, 5, 3),           rand(Bool, 4, 5, 3))),
)
    A, B = CartesianRunIndices.(masks)
    SA, SB = _to_set(A), _to_set(B)
    for (op, ref) in ((intersect, intersect), (setdiff, setdiff), (union, union))
        C = op(A, B)
        @test _to_set(C) == ref(SA, SB)
    end
end

@testset "binary ops: 1D edge cases" begin
    A = CartesianRunIndices(Bool[0, 1, 1, 0, 1, 1, 0, 1])
    E = CartesianRunIndices(falses(8))
    SA = _to_set(A)

    @test _to_set(intersect(A, A)) == SA
    @test isempty(setdiff(A, A))
    @test _to_set(union(A, A)) == SA

    @test isempty(intersect(A, E))
    @test _to_set(setdiff(A, E)) == SA
    @test _to_set(union(A, E))   == SA
    @test isempty(setdiff(E, A))
    @test _to_set(union(E, A))   == SA

    # adjacent disjoint runs → merged into one
    X = CartesianRunIndices(Bool[1, 1, 0, 0])
    Y = CartesianRunIndices(Bool[0, 0, 1, 1])
    @test isempty(intersect(X, Y))
    @test _to_set(setdiff(X, Y)) == _to_set(X)
    @test length(union(X, Y)) == 4

    # subset
    Sub = CartesianRunIndices(Bool[0, 1, 0, 0, 1, 0, 0, 0])
    @test _to_set(intersect(A, Sub)) == _to_set(Sub)
    @test isempty(setdiff(Sub, A))
    @test _to_set(union(A, Sub)) == SA
end


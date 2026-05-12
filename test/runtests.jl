using CartesianRuns
using Test

_to_set(cri) = Set(collect(cri))

# ──────────────────────────────────────────────────────────────────────────────
# CartesianRunIndices — construction and indexing
# ──────────────────────────────────────────────────────────────────────────────

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
        @test expand(CartesianRunIndices(mask)) == mask
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# Set operations — shared infrastructure
# ──────────────────────────────────────────────────────────────────────────────

@testset "domain check" begin
    a = CartesianRunIndices(Bool[1, 0, 1])
    b = CartesianRunIndices(Bool[1, 1, 0, 0])
    @test_throws ArgumentError intersect(a, b)
    @test_throws ArgumentError setdiff(a, b)
    @test_throws ArgumentError union(a, b)
end

# ──────────────────────────────────────────────────────────────────────────────
# complement
# ──────────────────────────────────────────────────────────────────────────────

@testset "complement: $label" for (label, mask) in (
    ("1D", Bool[0, 1, 1, 0, 1, 0]),
    ("2D", Bool[1 0 1; 0 1 0; 1 0 1]),
    ("3D", rand(Bool, 4, 5, 3)),
)
    cri = CartesianRunIndices(mask)
    c   = complement(cri)
    @test domain(c) == domain(cri)
    @test _to_set(c) == _to_set(CartesianRunIndices(.!mask))
    @test length(cri) + length(c) == prod(length.(domain(cri)))
    @test _to_set(complement(c)) == _to_set(cri)          # involution
end

@testset "complement: trivial" begin
    @test _to_set(complement(CartesianRunIndices(falses(5)))) ==
          _to_set(CartesianRunIndices(trues(5)))
    @test isempty(complement(CartesianRunIndices(trues(5))))
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
        @test domain(C) == domain(A)
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


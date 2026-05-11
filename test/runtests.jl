using CartesianRuns
using Test

# ===== CartesianRunIndices =====

@testset "CartesianRunIndices: 1D" begin
    mask = Bool[0, 1, 1, 0, 1]
    cri  = CartesianRunIndices(mask)
    @test length(cri) == 3
    @test size(cri)   == (3,)
    @test length(cri.intervals[1]) == 2
    @test cri.intervals[1][1].mask == 2:3
    @test cri.intervals[1][2].mask == 5:5
    @test cri.offsets === ()

    trues = findall(mask)
    for k in 1:length(cri)
        @test cri[k] == CartesianIndex(trues[k])
    end
end

@testset "CartesianRunIndices: 2D" begin
    mask = Bool[
        1 0 0 1
        1 1 0 1
        0 1 0 1
    ]
    cri   = CartesianRunIndices(mask)
    trues = findall(mask)
    @test length(cri) == 7
    @test size(cri)   == (7,)
    @test length(cri.intervals[1]) == 3   # x-runs across active columns
    @test length(cri.intervals[2]) == 2   # y-runs in [T,T,F,T]
    @test length(cri.offsets[1])   == 4   # 3 active y-cells + 1 sentinel
    @test cri.offsets[1][begin]    == 1
    @test cri.offsets[1][end]      == 4   # each active col contributed 1 x-run
    for k in 1:length(cri)
        @test cri[k] == trues[k]
    end
end

@testset "CartesianRunIndices: 3D" begin
    mask = falses(2, 3, 2)
    mask[1, 1, 1] = mask[2, 1, 1] = mask[1, 2, 1] = mask[2, 3, 2] = true
    cri   = CartesianRunIndices(mask)
    trues = findall(mask)
    @test length(cri) == 4
    @test size(cri)   == (4,)
    for k in 1:length(cri)
        @test cri[k] == trues[k]
    end
end

@testset "CartesianRunIndices: edge cases" begin
    # 1D all-true
    cri_t = CartesianRunIndices(Bool[1, 1, 1])
    @test length(cri_t) == 3
    @test cri_t[1] == CartesianIndex(1)
    @test cri_t[3] == CartesianIndex(3)

    # 1D all-false
    cri_f = CartesianRunIndices(Bool[0, 0, 0])
    @test length(cri_f) == 0
    @test isempty(cri_f)

    # 2D all-false
    @test length(CartesianRunIndices(falses(3, 4))) == 0

    # 2D single true cell
    mask      = falses(4, 4)
    mask[2,3] = true
    cri       = CartesianRunIndices(mask)
    @test length(cri) == 1
    @test cri[1] == CartesianIndex(2, 3)
end

@testset "CartesianRunIndices: collect matches findall" begin
    mask = Bool[
        1 0 1 0 1 0
        0 1 0 1 0 1
        1 1 0 0 1 1
        0 0 1 1 0 0
        1 0 0 0 1 0
    ]
    @test collect(CartesianRunIndices(mask)) == findall(mask)
end

@testset "CartesianRunIndices: Base.in 1D" begin
    mask = Bool[0, 1, 1, 0, 1, 0, 0, 1, 1, 1]
    cri  = CartesianRunIndices(mask)
    for i in eachindex(mask)
        @test (CartesianIndex(i) ∈ cri) == mask[i]
    end
    @test CartesianIndex(1) ∉ CartesianRunIndices(Bool[])
end

@testset "CartesianRunIndices: Base.in 2D" begin
    mask = Bool[1 0 1; 0 1 0; 1 1 0; 0 0 1]
    cri  = CartesianRunIndices(mask)
    for idx in CartesianIndices(mask)
        @test (idx ∈ cri) == mask[idx]
    end
    @test CartesianIndex(1, 1) ∉ CartesianRunIndices(falses(3, 3))
end

@testset "CartesianRunIndices: Base.in 3D" begin
    mask = rand(Bool, 4, 5, 3)
    cri  = CartesianRunIndices(mask)
    for idx in CartesianIndices(mask)
        @test (idx ∈ cri) == mask[idx]
    end
end

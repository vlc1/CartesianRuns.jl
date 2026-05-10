using CartesianRuns
using Test
using OffsetArrays

@testset "empty and boundary cases" begin
    @test isempty(find_runs(Bool[]))
    @test isempty(find_runs(Bool[0, 0, 0, 0]))
    @test isempty(find_runs(Bool[0]))
    @test isempty(find_runs(Bool[], 100))
end

@testset "single elements" begin
    runs = find_runs(Bool[1])
    @test length(runs) == 1
    @test runs[1].range == 1:1
    @test runs[1].shift == 0

    runs = find_runs(Bool[0, 0, 1, 0, 0])
    @test length(runs) == 1
    @test runs[1].range == 3:3
    @test runs[1].shift == -2
end

@testset "all true" begin
    runs = find_runs(Bool[1, 1, 1])
    @test length(runs) == 1
    @test runs[1].range == 1:3
    @test runs[1].shift == 0

    runs = find_runs(Bool[1, 1, 1], 100)
    @test runs[1].shift == 100
end

@testset "two disjoint runs" begin
    runs = find_runs(Bool[0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0])
    @test length(runs) == 2

    @test runs[1].range == 5:7
    @test runs[1].shift == -4

    @test runs[2].range == 9:10
    @test runs[2].shift == -5
end

@testset "two disjoint runs with prior=100" begin
    # prior=100 simulates 100 trues seen on prior scanlines
    runs = find_runs(Bool[0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0], 100)
    @test length(runs) == 2

    @test runs[1].range == 5:7
    @test runs[1].shift == 96

    @test runs[2].range == 9:10
    @test runs[2].shift == 95
end

@testset "alternating pattern" begin
    runs = find_runs(Bool[1, 0, 1, 0, 1])
    @test length(runs) == 3

    @test runs[1].range == 1:1
    @test runs[1].shift == 0

    @test runs[2].range == 3:3
    @test runs[2].shift == -1

    @test runs[3].range == 5:5
    @test runs[3].shift == -2
end

@testset "runs at start, middle, end" begin
    runs = find_runs(Bool[1, 1, 0, 1, 0, 1, 1])
    @test length(runs) == 3

    @test runs[1].range == 1:2
    @test runs[1].shift == 0

    @test runs[2].range == 4:4
    @test runs[2].shift == -1

    @test runs[3].range == 6:7
    @test runs[3].shift == -2
end

@testset "long runs" begin
    runs = find_runs(Bool[1, 1, 1, 1, 1, 0, 0, 0])
    @test length(runs) == 1
    @test runs[1].range == 1:5
    @test runs[1].shift == 0

    runs = find_runs(Bool[0, 0, 0, 1, 1, 1, 1])
    @test length(runs) == 1
    @test runs[1].range == 4:7
    @test runs[1].shift == -3
end

@testset "runs with single-element gaps" begin
    runs = find_runs(Bool[1, 1, 0, 1, 1])
    @test length(runs) == 2
    @test runs[1].range == 1:2
    @test runs[1].shift == 0
    @test runs[2].range == 4:5
    @test runs[2].shift == -1

    runs = find_runs(Bool[1, 1, 1, 0, 1, 1])
    @test length(runs) == 2
    @test runs[1].range == 1:3
    @test runs[1].shift == 0
    @test runs[2].range == 5:6
    @test runs[2].shift == -1
end

@testset "shift semantics with prior=0 (explicit)" begin
    runs = find_runs(Bool[0, 0, 1, 1, 0, 1], 0)
    @test length(runs) == 2

    @test runs[1].range == 3:4
    @test runs[1].shift == -2

    @test runs[2].range == 6:6
    @test runs[2].shift == -3
end

@testset "shift semantics with negative prior" begin
    runs = find_runs(Bool[1, 1, 0, 1], -10)
    @test length(runs) == 2

    @test runs[1].range == 1:2
    @test runs[1].shift == -10

    @test runs[2].range == 4:4
    @test runs[2].shift == -11
end

@testset "shift semantics with positive prior" begin
    runs = find_runs(Bool[0, 1, 1], 5)
    @test length(runs) == 1
    @test runs[1].range == 2:3
    @test runs[1].shift == 4  # 0 - 2 + 5 + 1 = 4
end

@testset "large sparse vector" begin
    runs = find_runs(Bool[zeros(50)..., 1, 1, zeros(30)..., 1, zeros(20)...])
    @test length(runs) == 2

    @test runs[1].range == 51:52
    @test runs[1].shift == -50

    @test runs[2].range == 83:83
    @test runs[2].shift == -80
end

@testset "many short runs" begin
    runs = find_runs(Bool[1, 0, 1, 0, 1, 0, 1, 0])
    @test length(runs) == 4

    @test runs[1].range == 1:1
    @test runs[1].shift == 0

    @test runs[2].range == 3:3
    @test runs[2].shift == -1

    @test runs[3].range == 5:5
    @test runs[3].shift == -2

    @test runs[4].range == 7:7
    @test runs[4].shift == -3
end

@testset "shift formula: n - start + prior + 1" begin
    # Verify: for any i ∈ run.range, i + run.shift is the 1-based linear
    # position of mask[i] in the compact data array.
    runs = find_runs(Bool[1, 1, 0, 1])
    @test length(runs) == 2

    # First run: n=0, start=1, prior=0 => shift = 0 - 1 + 0 + 1 = 0
    @test runs[1].shift == 0

    # Second run: n=2, start=4, prior=0 => shift = 2 - 4 + 0 + 1 = -1
    @test runs[2].shift == -1
end

@testset "complex mixed-length pattern" begin
    runs = find_runs(Bool[1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1])
    @test length(runs) == 4

    @test runs[1].range == 1:1
    @test runs[1].shift == 0

    @test runs[2].range == 3:4
    @test runs[2].shift == -1

    @test runs[3].range == 7:9
    @test runs[3].shift == -3

    @test runs[4].range == 11:11
    @test runs[4].shift == -4
end

@testset "trailing run with positive prior" begin
    # Single open run that closes only at end-of-array, plus a non-zero prior:
    # exercises the trailing-push path under a non-default offset.
    runs = find_runs(Bool[0, 0, 1, 1, 1], 7)
    @test length(runs) == 1
    @test runs[1].range == 3:5
    @test runs[1].shift == 5  # 0 - 3 + 7 + 1 = 5
end

@testset "BitVector input" begin
    mask = BitVector([0, 0, 1, 1, 0, 1])
    runs = find_runs(mask)
    @test length(runs) == 2
    @test runs[1].range == 3:4
    @test runs[1].shift == -2
    @test runs[2].range == 6:6
    @test runs[2].shift == -3
end

@testset "OffsetArray input" begin
    # Indices -2..2; values at -2,-1,0,1,2 are 0,1,1,0,1
    mask = OffsetArray(Bool[0, 1, 1, 0, 1], -2:2)
    runs = find_runs(mask)
    @test length(runs) == 2

    # Run 1: range=-1:0, n=0, start=-1, shift = 0 - (-1) + 0 + 1 = 2
    @test runs[1].range == -1:0
    @test runs[1].shift == 2

    # Run 2: range=2:2, n=2, start=2, shift = 2 - 2 + 0 + 1 = 1
    @test runs[2].range == 2:2
    @test runs[2].shift == 1
end

@testset "1-based invariant: compact[i + shift] recovers cell" begin
    # Lock in the contract that i + run.shift indexes a 1-based compact array.
    # Build a hand-rolled compact buffer that stores the original index of each
    # true cell, in scan order; then verify compact[i + shift] == i.
    mask = Bool[0, 0, 1, 1, 0, 1, 1, 1]
    compact = [i for i in eachindex(mask) if mask[i]]
    @test compact == [3, 4, 6, 7, 8]

    runs = find_runs(mask)
    for run in runs, i in run.range
        @test compact[i + run.shift] == i
    end

    # Also exercise the invariant under a non-zero prior, where compact is
    # offset by `prior` cells worth of "earlier scanline" content.
    prior = 4
    compact_with_prior = vcat(fill(-1, prior), compact)
    runs = find_runs(mask, prior)
    for run in runs, i in run.range
        @test compact_with_prior[i + run.shift] == i
    end
end

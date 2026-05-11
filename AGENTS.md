# AGENTS.md

## Project goal

Build a Julia package providing a read-only `AbstractVector{CartesianIndex{N}}`
view of the positions where a boolean mask is `true`, using SAMURAI-style
interval compression for arbitrary `N`.

- Package name: **`CartesianRuns.jl`**.
- Exported type (planned): **`CartesianRunIndices{N} <: AbstractVector{CartesianIndex{N}}`**.
- Currently exported: `Interval`, `push_runs!`.

Reference: SAMURAI ([repo](https://github.com/hpc-maths/samurai),
[docs](https://hpc-math-samurai.readthedocs.io/)).

## Files

| File                          | Role                                          |
| ----------------------------- | --------------------------------------------- |
| `src/CartesianRuns.jl`        | Module entry; exports + `include`s            |
| `src/types.jl`                | `Interval` struct                             |
| `src/runs.jl`                 | 1D `push_runs!` + N-D stub                    |
| `test/runtests.jl`            | Test suite (run via `Pkg.test()`)             |
| `Project.toml`                | Deps: `OffsetArrays`; test extras: `Test`     |
| `.github/workflows/CI.yml`    | CI on Julia 1.10, 1.12, pre (Ubuntu)          |
| `.github/workflows/TagBot.yml`| Release tagging                               |
| `PLAN.md`                     | Short scratch notes for the next refactor     |

## Sticky decisions (do not re-litigate)

### `Interval`

- Plain struct, no type parameters: `range::UnitRange{Int}`, `shift::Int`.
- **Invariant**: for any `i ∈ run.range`, `i + run.shift` is the 1-based
  linear position of that cell in the compact data array.
- Field is `shift`, not `at`. A SAMURAI-style `@` operator/function can be
  added later as a separate API surface without renaming the field.

### `push_runs!(runs::AbstractVector{Interval}, mask::AbstractVector{Bool}, prior::Int)`

- Iterator-based via `pairs(IndexLinear(), mask)`. Not index-scan, not the
  dual-iterator pattern from the original.
- State: `prev` (was previous cell true?), `start` (open run start),
  `n` (cumulative trues in this scan), `m` (runs pushed so far).
- Two `push!` sites (transition close + guarded trailing close). Accepted as
  the honest cost: a closure helper would box `n`; an `Iterators.flatten`
  sentinel isn't type-stable.
- `prior::Int` is the cumulative count of trues from prior scanlines (no
  default; always passed explicitly). With `prior = 0`, the first true cell
  maps to compact position `1` (Julia's 1-based convention).
- Allocation is the caller's responsibility: `similar(mask, Interval, 0)`.
- Returns `(n, m)::Tuple{Int,Int}` — cumulative trues and run count for this
  scan. `n` is passed as `prior` to the next scanline in N-D contexts.
- No `@inbounds` — the iterator path already avoids bounds checks.
- Shift formula: `shift = n - start + prior + 1`.

### N-D direction

- Stub already in `src/runs.jl`:
  `push_runs!(runs, mask::AbstractArray{Bool,N}, line::NTuple{N-1,Int}, prior)`.
- Plan: scan along dim 1 per scanline; `line` carries the (N-1) indices into
  dims 2..N; `prior` propagates cumulative trues between scanlines;
  per-scanline offsets stored externally (SAMURAI layout).
- GPU support deferred.

## Pending work

1. Implement N-D `push_runs!`: scanline loop + offset table.
2. Build `CartesianRunIndices{N}` on top: `size`, `getindex`, iteration
   yielding `CartesianIndex{N}`.
3. Optional SAMURAI-flavored API: `at(::Interval, ::Int) = i + shift` plus
   a `@` operator/macro.

## Conventions

- Variable names in `push_runs!`: `prev`, `start`, `n`, `m`, `prior`.
- 1D mask: `AbstractVector{Bool}`. N-D mask: `AbstractArray{Bool,N}`.
- Run tests with `julia --project=. -e 'using Pkg; Pkg.test()'` from the
  package root (`~/.julia/dev/CartesianRuns`).

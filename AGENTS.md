# AGENTS.md

## Project goal

Build a Julia package providing a read-only `AbstractVector{CartesianIndex{N}}`
view of the positions where a boolean mask is `true`, using SAMURAI-style
interval compression for arbitrary `N`.

- Package name: **`CartesianRuns.jl`**.
- Exported: `Interval`, `CartesianRunIndices{N} <: AbstractVector{CartesianIndex{N}}`.

Reference: SAMURAI ([repo](https://github.com/hpc-maths/samurai),
[docs](https://hpc-math-samurai.readthedocs.io/)).

## Files

| File                          | Role                                                        |
| ----------------------------- | ----------------------------------------------------------- |
| `src/CartesianRuns.jl`        | Module entry; exports + `include`s                          |
| `src/types.jl`                | `Interval` struct; `CartesianRunIndices{N}` + `AbstractVector` interface |
| `src/construction.jl`         | `_build_runs!` (1D); `_build_fused!`; `_build_cartesian_runs` |
| `test/runtests.jl`            | Test suite (run via `Pkg.test()`)                           |
| `docs/make.jl`                | Documenter build + deploy script                            |
| `docs/src/index.md`           | Single-page API reference                                   |
| `Project.toml`                | Deps: `OffsetArrays`; test extras: `Test`                   |
| `.github/workflows/CI.yml`    | CI on Julia 1.10, 1.12, pre (Ubuntu)                        |
| `.github/workflows/Docs.yml`  | Docs build + deploy to GitHub Pages                         |
| `.github/workflows/TagBot.yml`| Release tagging                                             |
| `PLAN.md`                     | Short scratch notes for the next refactor                   |

## Sticky decisions (do not re-litigate)

### `Interval`

- Plain struct, no type parameters: `mask::UnitRange{Int}`, `compact::UnitRange{Int}`.
- Single constructor `Interval(mask, shift)` derives `compact = mask .+ shift`;
  invariant `length(mask) == length(compact)` is enforced by construction.
- `shift(iv::Interval)` (exported) returns `iv.compact.start - iv.mask.start`.

### `_build_runs!(runs, mask::AbstractVector{Bool}, prior::Int)` (internal)

- Iterator-based via `pairs(IndexLinear(), mask)`.
- State: `prev`, `start`, `n` (cumulative trues), `m` (runs pushed).
- `prior`: cumulative trues from prior scanlines; pass `0` for standalone 1D.
- Returns `(n, m)::Tuple{Int,Int}`. Shift formula: `shift = n - start + prior + 1`.
- No `@inbounds` — iterator avoids bounds checks.

### `CartesianRunIndices{N}`

- Stores `intervals::NTuple{N, AbstractVector{Interval}}` and
  `offsets::NTuple{N-1, AbstractVector{Int}}` (CSR, 1-based, pre-seeded with `[1]`).
- Construction: `_build_cartesian_runs` pre-allocates all buffers, then fills
  them in a single fused pass via `_build_fused!` (no intermediate allocations).
- `_build_fused!` is a two-method recursive function: D=1 base calls
  `_build_runs!`; D≥2 peels the outermost axis, keeping all run-detection state
  (`prev`, `start`, `n_d`, `m_d`, `inner_prior`) as local stack variables.
- `getindex(cri, k)`: binary-search `_search_compact` to locate the x-interval,
  then recurse via `_recover` through the offset tables to reconstruct
  `CartesianIndex{N}`.

### N-D construction — key invariants

- `offsets[d][p-1]+1 : offsets[d][p]` is the range of dim-`d` intervals that
  belong to the `p`-th active dim-`(d+1)` position.
- Dim-1 prior (`x_prior`) threads through ALL recursive calls; each level's own
  prior (`inner_prior`) is a local accumulator in that level's stack frame.
- GPU support deferred.

## Pending work

1. Optional SAMURAI-flavored API: `at(::Interval, ::Int) = i + shift` plus
   a `@` operator/macro.

## Conventions

- Variable names in `_build_runs!`: `prev`, `start`, `n`, `m`, `prior`.
- Variable names in `_build_fused!`: same, plus `n_d`, `m_d`, `inner_prior`, `d_prior`, `x_prior`.
- 1D mask: `AbstractVector{Bool}`. N-D mask: `AbstractArray{Bool,N}`.
- Run tests: `julia --project=. -e 'using Pkg; Pkg.test()'` from the package root.
- Build docs: `julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=".")); Pkg.instantiate()'`
  then `julia --project=docs/ docs/make.jl`.

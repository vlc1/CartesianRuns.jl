# AGENTS.md

## Project goal

Build a Julia package providing a read-only `AbstractVector{CartesianIndex{N}}`
view of the positions where a boolean mask is `true`, using SAMURAI-style interval
compression for arbitrary `N`.

- Package name: **`CartesianRuns.jl`**.
- Exported: `Interval`, `CartesianRunIndices{N} <: AbstractVector{CartesianIndex{N}}`,
  `shift`, `expand`, `complement`, `intersect`, `setdiff`, `union`.

Reference: [SAMURAI](https://github.com/hpc-maths/samurai)
([docs](https://hpc-math-samurai.readthedocs.io/)).

## Files

| File                          | Role                                                        |
| ----------------------------- | ----------------------------------------------------------- |
| `src/CartesianRuns.jl`        | Module entry; exports + `include`s                          |
| `src/types.jl`                | `Interval` struct; `CartesianRunIndices{N}` + `AbstractVector` interface |
| `src/construction.jl`         | `_build_runs!` (1D); `_build_fused!`; `_build_cartesian_runs` |
| `src/common.jl`               | Shared helpers: `_inner_slice`                                  |
| `src/expand.jl`          | `expand`: `_expand_runs!`, `_expand_fused!`      |
| `src/complement.jl`           | `complement`: `_complement_runs!`, `_complement_fused!`         |
| `src/intersect.jl`            | `intersect`: `_intersect_runs!`, `_intersect_fused!`            |
| `src/setdiff.jl`              | `setdiff`: `_setdiff_runs!`, `_setdiff_fused!`                  |
| `src/union.jl`                | `union`: `_union_runs!`, `_union_fused!`                        |
| `ext/OffsetArraysExt.jl` | Package extension: `expand` for non-1-based domains (requires `using OffsetArrays`) |
| `docs/make.jl`                | Documenter build + deploy script                            |
| `docs/src/index.md`           | Single-page API reference                                   |
| `Project.toml`                | Weak dep: `OffsetArrays`; test extras: `Test`, `OffsetArrays`   |
| `test/runtests.jl`            | Test suite (run via `Pkg.test()`)                           |
| `.github/workflows/CI.yml`    | CI on Julia 1.10, 1.12, pre (Ubuntu); compat minimum is 1.9  |
| `.github/workflows/Docs.yml`  | Docs build + deploy to GitHub Pages                         |
| `.github/workflows/TagBot.yml`| Release tagging                                             |

## Sticky decisions (do not re-litigate)

### `Interval`

- Plain struct, no type parameters: `mask::UnitRange{Int}`, `compact::UnitRange{Int}`.
- Single constructor `Interval(mask, shift)` derives `compact = mask .+ shift`;
  invariant `length(mask) == length(compact)` is enforced by construction.
- `shift(iv::Interval)` (exported) returns `iv.compact.start - iv.mask.start`.

### `_build_runs!(runs, mask::AbstractVector{Bool}, range::AbstractUnitRange{Int}, prior::Int)` (internal)

- Iterates `for i in range; val = mask[i]` — `i` is always in original mask-space coordinates.
- State: `prev`, `start`, `n` (cumulative trues), `m` (runs pushed).
- `range`: subset of `axes(mask, 1)` to scan; passing `axes(mask, 1)` gives the full scan.
- `prior`: cumulative trues from prior scanlines; pass `0` for standalone 1D.
- Returns `(n, m)::Tuple{Int,Int}`. Shift formula: `shift = n - start + prior + 1`.
- Backward-compat 3-arg wrapper: `_build_runs!(runs, mask, prior)` passes `axes(mask, 1)`.

### `CartesianRunIndices{N}`

- Stores `intervals::NTuple{N, AbstractVector{Interval}}` and
  `offsets::NTuple{N-1, AbstractVector{Int}}` (CSR, 1-based, pre-seeded with `[1]`).
- Construction: `_build_cartesian_runs(mask[, domain])` pre-allocates all buffers,
  then fills them in a single fused pass via `_build_fused!(ivs, offs, 0, 0, mask, domain)`
  (no intermediate allocations).
- `_build_fused!` is a two-method recursive function dispatching on NTuple length
  of `domain`: D=1 base calls `_build_runs!(ivs[1], mask, domain[1], x_prior)`;
  D≥2 iterates `domain[D]`, keeps all run-detection state (`prev`, `start`, `n_d`,
  `m_d`, `inner_prior`) as local stack variables, and recurses with `Base.front(domain)`.
- `getindex(cri, k)`: binary-search `_search_compact` to locate the x-interval,
  then recurse via `_recover` through the offset tables to reconstruct
  `CartesianIndex{N}`.

#### Domain-restricted constructor

```julia
CartesianRunIndices(mask, domain::NTuple{N,AbstractUnitRange{Int}})
```

Restricts the scan to `true` cells within `domain`. Implementation: validates
that each `domain[d]` is a subset of `axes(mask, d)`, converts to plain
`UnitRange{Int}`, then calls `_build_cartesian_runs(mask, dom)` which passes
`dom` as the `domain` NTuple through `_build_fused!`, so the kernels iterate
only over `domain[d]` without any array wrapping or post-processing.
The stored `domain` field is gone; the restricted constructor calls
`_build_cartesian_runs(mask, dom)` directly and does not retain `dom`.
Validation: `ArgumentError` if any `domain[d]` is out of `axes(mask, d)`.

### N-D construction — key invariants

- `offsets[d][p-1]+1 : offsets[d][p]` is the range of dim-`d` intervals that
  belong to the `p`-th active dim-`(d+1)` position.
- Dim-1 prior (`x_prior`) threads through ALL recursive calls; each level's own
  prior (`inner_prior`) is a local accumulator in that level's stack frame.
- GPU support deferred.

## Pending work

1. Optional SAMURAI-flavored API: `at(::Interval, ::Int) = i + shift` plus
   a `@` operator/macro.

## Set operations, complement and expand

`intersect`, `setdiff`, `union`, `complement`, and `expand` each live in
their own source file (`src/intersect.jl`, `src/setdiff.jl`, `src/union.jl`,
`src/complement.jl`, `src/expand.jl`). The shared helper `_inner_slice` is
in `src/common.jl`, which is included first.

`complement` and `expand` take an explicit `domain` argument. `expand` requires `domain` to be
`NTuple{N,Base.OneTo{Int}}` (i.e., 1-based); for non-1-based domains load
`OffsetArrays.jl`, which activates the package extension `OffsetArraysExt`
and returns an `OffsetArray{Bool,N}`. `complement` accepts any
`NTuple{N,<:AbstractUnitRange{Int}}` without `OffsetArrays` because it
operates entirely in mask-space coordinates and returns a `CartesianRunIndices`.
Binary operations do not check domain compatibility — that is the caller's
responsibility.

### Three-layer dispatch pattern (all operations)

| Layer | Name pattern | Role |
|-------|-------------|------|
| Public API | `op(cri, ...)` | Buffer allocation, single top-level call |
| 1D kernel | `_op_runs!(out, a_ivs, a_lo, a_hi, ..., prior)` | Flat sweep on one slice of `AbstractVector{Interval}` |
| N-D kernel | `_op_fused!(out_ivs, out_offs, a_ivs, a_offs, ...)` | Recursive dimensional peeling; dispatches to 1D base at bottom |

Dimension peeling passes cost-free `SubArray`s (via `view`) rather than copies.

### `_op_fused!` dispatch

Two methods, selected by Julia on tuple length:

- **1D base**: `a_ivs::NTuple{1,...}`, `a_offs::Tuple{}` → calls `_op_runs!` on `a_ivs[1]`, returns `(x_prior + n, n, m)`.
- **N-D case**: `a_ivs::NTuple{N,...}`, `a_offs::NTuple{M,...}` where N,M≥2 → peels last dimension:
  - `a_outer = last(a_ivs)`, `a_off = last(a_offs)` — outermost interval vector + its CSR offsets.
  - `front_a = Base.front(a_ivs)`, `front_ao = Base.front(a_offs)` — passed to the recursive call.
  - `_inner_slice(iv, offs, r)` maps row `r` within outer interval `iv` to a `(lo, hi)` slice of the inner interval vector.
  - Outer run detection: `prev`, `start`, `n_d`, `m_d`; shift formula `n_d - start + d_prior + 1`.

### Loop structure by operation

| Operation | Category | Outer loop | Recurse when |
|-----------|----------|-----------|--------------|
| `_build_fused!` | unary | `for r in domain[D]` — all rows | always |
| `_expand_fused!` | unary | same all-rows loop | always |
| `_complement_fused!` | unary | `for r in last(domain)` — all rows; single pointer `i` advances through A | always |
| `_intersect_fused!` | binary | `while ar < typemax \|\| br < typemax` | both A and B cover `r` |
| `_setdiff_fused!` | binary | same ar/br | A covers `r` |
| `_union_fused!` | binary | same ar/br | always (at least one covers `r`) |

### Key invariants

- **`x_prior`**: cumulative dim-1 (x) compact count, threaded through every recursive call as the first return element.  Each level's `inner_prior` is a local accumulator for the level below.
- **Gap check** (binary ops only): at the TOP of the `while` loop body, before computing inner slices: `if prev && r > last_r + 1` → close the open outer run at `last_r`. Prevents merging runs across rows skipped by the ar/br sweep.
- **`last_r`** is updated unconditionally on every visited row; `prev` guards whether the gap check fires.

## Conventions

- Variable names in `_build_runs!`: `prev`, `start`, `n`, `m`, `prior`.
- Variable names in `_build_fused!`: same, plus `n_d`, `m_d`, `inner_prior`, `d_prior`, `x_prior`.
- 1D mask: `AbstractVector{Bool}`. N-D mask: `AbstractArray{Bool,N}`.
- Run tests: `julia --project=. -e 'using Pkg; Pkg.test()'` from the package root.
- Build docs: `julia --project=docs/ -e 'using Pkg; Pkg.develop(PackageSpec(path=".")); Pkg.instantiate()'`
  then `julia --project=docs/ docs/make.jl`.

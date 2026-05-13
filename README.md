# CartesianRuns.jl

[![docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://vlc1.github.io/CartesianRuns.jl/dev/)

A Julia package providing a read-only `AbstractVector{CartesianIndex{N}}`
view of the positions where a boolean mask is `true`, using
[SAMURAI](https://github.com/hpc-maths/samurai)-style interval compression for
arbitrary `N`.

Two coordinate spaces are distinguished throughout:

- **Mask space**: the Cartesian indices of the original array (e.g.
  `CartesianIndex(i, j)` for a 2D mask). Contains gaps — the `false` cells.
- **Compact space**: the 1-based linear positions 1…`count(mask)`,
  enumerating only the `true` cells in column-major order. No gaps.

`CartesianRunIndices` exposes three fundamental entry points:

| Operation                    | Direction                          | Julia API                   |
|------------------------------|------------------------------------|-----------------------------|
| `CartesianRunIndices(mask)`  | boolean mask → compressed index set | constructor                |
| `expand(cri, axes(mask))`    | compressed index set → boolean mask | `CartesianRuns.expand`     |
| `cri[k]`                     | compact → mask space               | `Base.getindex`             |
| `ci ∈ cri`                   | mask → compact (membership test)   | `Base.in`                   |

`CartesianRunIndices(mask)` and `expand(cri, axes(mask))` are exact inverses:
`expand(CartesianRunIndices(mask), axes(mask)) == mask`. Both are backed by
`Interval` objects — each storing a contiguous mask-space run
(`mask::UnitRange{Int}`) alongside the corresponding compact-space run
(`compact::UnitRange{Int}`), always of equal length.

Set-algebraic operations return a new `CartesianRunIndices` without ever
materialising a boolean mask:

| Operation               | Returns                                   | Julia API                  |
|-------------------------|-------------------------------------------|----------------------------|
| `complement(A, domain)` | positions in `domain` not in `A`          | `CartesianRuns.complement` |
| `intersect(A, B)`       | positions in both `A` and `B`             | `Base.intersect`           |
| `union(A, B)`           | positions in `A`, `B`, or both            | `Base.union`               |
| `setdiff(A, B)`         | positions in `A` but not in `B`           | `Base.setdiff`             |

## Example

```julia
using CartesianRuns

mask = Bool[0 1 0; 1 1 0; 0 1 1]
cri  = CartesianRunIndices(mask)
collect(cri) == findall(mask)      # true — compact → mask
CartesianIndex(2, 1) ∈ cri        # true  — mask → compact membership
```

## How operations work

Every operation on `CartesianRunIndices` is implemented in three layers:

```
Public API    allocates output buffers, calls ↓
_op_fused!    recursive N-D kernel, calls ↓ at the bottom
_op_runs!     1D sweep on a flat slice of intervals
```

The key idea is **dimensional peeling**. An N-D boolean mask can be thought
of as a stack of (N-1)-D slices. Each slice either contributes to an active
run in the outermost dimension or it doesn't. `_op_fused!` exploits this
recursively:

1. Look at the outermost dimension. Walk its intervals (one per active
   contiguous group of slices).
2. For each row in that dimension, call `_op_fused!` on the next dimension
   inward — passing a narrower slice of the interval tables as a cost-free
   `SubArray` via `view`.
3. When only one dimension remains, call `_op_runs!`, which does a plain
   1D sweep.

The recursion bottoms out via Julia dispatch: `_op_fused!` has two methods
distinguished by the length of the interval tuple — `NTuple{1,...}` +
`Tuple{}` offsets dispatches to the 1D base; `NTuple{N,...}` +
`NTuple{M,...}` dispatches to the recursive case.

### Unary operations

**Construction** (`CartesianRunIndices(mask)`), **expand**
(`expand(cri, domain)`), and **complement** (`complement(cri, domain)`) all
iterate *every* row of the outermost axis with a single advancing pointer
through the input's interval table.

- `CartesianRunIndices(mask)` builds the interval representation from a
  boolean mask (`_build_fused!` / `_build_runs!`).
- `expand(cri, domain)` is the exact inverse: it allocates a `Bool` array via
  `similar` (with shape `length.(domain)`), fills it `false`, then writes
  `true` for each interval's mask range (`_expand_fused!` / `_expand_runs!`).
- `complement(cri, domain)` produces a new `CartesianRunIndices` covering every
  row in `domain` that is *not* covered by the input
  (`_complement_fused!` / `_complement_runs!`).

### Binary operations

**`intersect`**, **`setdiff`**, and **`union`** use an `ar/br` two-pointer
sweep that skips rows covered by neither operand — more efficient, but
requires an extra invariant: when a gap is detected (`r > last_r + 1` while
`prev` is `true`), the currently open outer run must be closed before
advancing. Without this check, runs across non-adjacent rows would be
silently merged.

## Documentation

Build locally with:

```
julia --project=docs/ -e \
  'using Pkg; Pkg.develop(PackageSpec(path=".")); Pkg.instantiate()'
julia --project=docs/ docs/make.jl
```

Inspired by [SAMURAI](https://github.com/hpc-maths/samurai).

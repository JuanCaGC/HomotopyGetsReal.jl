# HomotopyGetsReal

[![CI](https://github.com/JuanCaGC/HomotopyGetsReal.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuanCaGC/HomotopyGetsReal.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/JuanCaGC/HomotopyGetsReal.jl/actions/workflows/Documentation.yml/badge.svg)](https://JuanCaGC.github.io/HomotopyGetsReal.jl)

Julia reimplementation of **Homotopy gets real** — numerical algebraic geometry for decomposing real algebraic curves and surfaces using [HomotopyContinuation.jl](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl).

Full documentation: [https://JuanCaGC.github.io/HomotopyGetsReal.jl](https://JuanCaGC.github.io/HomotopyGetsReal.jl)

## Install

Registered in Julia's General registry as v0.2.1, including the full feature set described below (isosingular deflation, generic projection support, incidence-based mesh welding).

```julia
using Pkg
Pkg.add("HomotopyGetsReal")
```

Or develop a local clone:

```julia
using Pkg
Pkg.activate("/path/to/HomotopyGetsReal")  # or `] dev /path/to/HomotopyGetsReal`
Pkg.instantiate()
```

Requires **Julia 1.12+** (enforced by `Project.toml` compat on `julia` / `LinearAlgebra`). Licensed under the MIT License (see `LICENSE`).

## Quick start — unit sphere end-to-end

```julia
using HomotopyContinuation
using HomotopyGetsReal
using GLMakie

@var x y z
F = System([x^2 + y^2 + z^2 - 1], variables = [x, y, z])
cfg = HomotopyConfig{Float64}()  # defaults: edge_sample_density=50, midslice_sample_density=100

vertices, edges, faces, mesh = decompose_3d_surface(F, cfg)

fig = plot_surface_decomposition(mesh; color_by = :z, cfg = cfg)
GLMakie.save("sphere.png", fig)
```

For a 2D slice only (standalone — does not depend on the block above):

```julia
using HomotopyContinuation
using HomotopyGetsReal
using GLMakie

@var x y z
F = System([x^2 + y^2 + z^2 - 1], variables = [x, y, z])
cfg = HomotopyConfig{Float64}()

v2d, e2d = slice_at_z(F, 0.0, cfg)
fig2d = plot_curve_decomposition(v2d, e2d; cfg = cfg)
GLMakie.save("sphere_equator.png", fig2d)
```

## Testing

**Fast suite** (default, skips slow Taubin integration):

```julia
using Pkg
Pkg.test("HomotopyGetsReal")
```

**Full suite** (adds `test_taubin.jl` and other slow-gated tests; 537 tests total, may vary by ±1 run-to-run (see `docs/DESIGN_NOTES.md`), ~30-34 min in this environment as of 2026-08):

```bash
HOMOTOPYGETSREAL_RUN_SLOW_TESTS=1 julia --project -e 'using Pkg; Pkg.test()'
```

Historical validation scripts from the Phases 1–6 rebuild live under `dev/scratch/` (superseded by `test/`).

## Layout

| Path | Purpose |
|------|---------|
| `src/` | Package source (flat module) |
| `test/` | Formal `Test.jl` suite |
| `dev/scratch/` | Archived scratch-phase validation scripts |
| `docs/` | Documenter.jl site (`docs/make.jl`) |

## AI-assisted development

Portions of this package's implementation, testing, and documentation were developed with the assistance of AI coding tools (Claude and Cursor), used under direct human supervision for code generation, refactoring, and independent review and auditing of the author's own design decisions. Commits with substantial AI-assisted content carry a `Co-Authored-By` trailer in the git history. All core design decisions were made by the author, and every correctness claim in this repository is backed by the verification methods described in `docs/DESIGN_NOTES.md` (residual checks, cross-validation against published ground-truth examples, and the package's own test suite) -- not by the AI tools' own output alone.

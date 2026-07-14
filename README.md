# HomotopyGetsReal

Julia reimplementation of **Homotopy gets real** — numerical algebraic geometry for decomposing real algebraic curves and surfaces using [HomotopyContinuation.jl](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl).

## Install

Not yet registered in Julia’s General registry. Install from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/JuanCaGC/HomotopyGetsReal.jl")
```

Or develop a local clone:

```julia
using Pkg
Pkg.activate("/path/to/HomotopyGetsReal")  # or `] dev /path/to/HomotopyGetsReal`
Pkg.instantiate()
```

Requires **Julia 1.12+** (enforced by `Project.toml` compat on `julia` / `LinearAlgebra` / `Statistics`). Licensed under the MIT License (see `LICENSE`).

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

**Full suite** (adds `test_taubin.jl`, ~77s for Taubin alone in this environment):

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
| `docs/` | Minimal docs pointer (API lives in docstrings) |

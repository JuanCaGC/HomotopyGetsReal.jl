# HomotopyGetsReal documentation

Built with [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) from package docstrings.

## Build locally

From the repository root (needs a normal desktop session so `GLMakie` can load):

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Then open `docs/build/index.html`.

## CI / headless Linux

`HomotopyGetsReal` loads `GLMakie` at package load time. On GitHub Actions (Ubuntu), follow Makie's [headless](https://docs.makie.org/stable/explanations/headless) / docs CI pattern: install `xvfb` + OpenGL packages and wrap **both** `Pkg.instantiate` and `docs/make.jl` with `xvfb-run` (see `.github/workflows/Documentation.yml`).

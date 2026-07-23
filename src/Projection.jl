# src/Projection.jl
#
# Phase 8: generic-projection support for decompose_3d_surface. The pipeline
# itself is coordinate-free (every frame-dependent decision flows through
# cfg.bbox_* or positional variable indices), so generic projections are a
# change-of-coordinates WRAPPER at a single boundary: rotate the system into a
# working chart, run the existing pipeline unchanged, map every output back to
# world coordinates. Nothing in Phase 1-7 reads world coordinates directly.
#
# Why this exists (the crash class): decompose_3d_surface's auto-augmented
# critical system {f, ∂f/∂x, ∂f/∂y} contains an identically-zero polynomial
# whenever the surface is independent of an augmenting axis (e.g. z - x^2, or
# any y-independent surface), and HomotopyContinuation's polyhedral start
# system then throws `OverflowError: Cannot compute a start system` from deep
# inside solve. A generic rotation makes that configuration measure-zero;
# `_verify_projection_ok` makes the residual non-generic cases fail loudly
# with an ArgumentError instead (BertiniReal's verify_projection_ok analogue).
#
# Orientation: projections are restricted to SO(3) (det = +1), not general
# O(3) -- weld_mesh aligns triangle windings with +∇f in the CHART, and for an
# orthogonal map cross(Qa, Qb) = det(Q)·Q·cross(a, b) while gradients map as
# ∇f_world = Q·∇f_chart, so a reflection (det = -1) would silently flip every
# mapped-back normal against the winding convention. `_resolve_projection`
# therefore rejects det < 0 matrices rather than auto-fixing them.

using LinearAlgebra
using Random

"""
    random_orthogonal_matrix(::Type{T}, n::Int; rng = Random.default_rng()) where {T<:AbstractFloat}
        -> Matrix{T}

Draw a Haar-uniform random rotation matrix in SO(n).

Constructed by QR decomposition of a Gaussian random matrix with the sign
correction that makes the distribution Haar-uniform on O(n), then mapped onto
SO(n). Distributionally equivalent to the Stewart Householder-reflection
construction BertiniReal/Bertini1 use (`make_matrix_random_real_mp`). Pass a
seeded `rng` for reproducibility.
"""
function random_orthogonal_matrix(::Type{T}, n::Int; rng::Random.AbstractRNG = Random.default_rng()) where {T<:AbstractFloat}
    # Plain LAPACK qr() of a Gaussian matrix is NOT Haar-distributed (its sign
    # convention biases the distribution); multiplying columns by the signs of
    # diag(R) restores Haar measure on O(n) (Mezzadri, Notices AMS 54 (2007)).
    # Negating one column when det < 0 is a measure-preserving right-translation
    # of the det=-1 coset onto SO(n), so uniformity survives. Verified
    # empirically (20k samples): det always +1, ||mean(column)|| consistent
    # with the 1/sqrt(N) null, E[(col)_z^2] = 0.3323 vs 1/3 for uniform.
    # Always generated in Float64 and converted to T, so chart systems built
    # from these entries keep exactly-Float64-representable coefficients (the
    # same precision boundary as path tracking; see Solver.jl's module header).
    A = randn(rng, n, n)
    Fq = qr(A)
    Q = Matrix(Fq.Q)
    signs = [r < 0 ? -1.0 : 1.0 for r in diag(Fq.R)]
    Q = Q .* signs'
    det(Q) < 0 && (Q[:, 1] .= -Q[:, 1])
    return Matrix{T}(Q)
end

"""
    _resolve_projection(projection, rng, cfg::HomotopyConfig{T}) where {T<:AbstractFloat} -> Matrix{Float64}

Build or validate the projection matrix for [`decompose_3d_surface`](@ref):
`:random` draws a fresh Haar-uniform SO(3) rotation from `rng`; a 3x3
orthonormal, det = +1 matrix is accepted as-is (converted to Float64 -- the
substitution precision path tracking uses anyway); anything else throws an
`ArgumentError`. Reflections (det < 0) are rejected rather than auto-fixed;
see the module header for the winding-convention rationale.

`cfg` (2026-07-23, Audit 1 Item 4/2a fix): the orthonormality acceptance
threshold now comes from `cfg.projection_orthonormality_tol` instead of a
bare `1e-8` literal with no way for a caller to source it from `cfg` at
all -- same default value as before, genuinely configurable now. Compared
against the Float64-computed `ortho_defect` via `Float64(cfg....)`, matching
this codebase's established `btol64`-style pattern (see
`intersect_bounding_object`) for a T-generic tolerance gating a strictly
Float64-only check.
"""
function _resolve_projection(projection, rng::Random.AbstractRNG, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    if projection === :random
        return random_orthogonal_matrix(Float64, 3; rng = rng)
    elseif projection isa AbstractMatrix
        size(projection) == (3, 3) || throw(ArgumentError(
            "decompose_3d_surface: projection matrix must be 3x3; got size $(size(projection)).",
        ))
        Q = Matrix{Float64}(projection)
        ortho_defect = norm(Q' * Q - I)
        ortho_tol64 = Float64(cfg.projection_orthonormality_tol)
        ortho_defect <= ortho_tol64 || throw(ArgumentError(
            "decompose_3d_surface: projection matrix is not orthonormal " *
            "(norm(Q'Q - I) = $(ortho_defect) > $(ortho_tol64)).",
        ))
        det(Q) > 0 || throw(ArgumentError(
            "decompose_3d_surface: projection matrix has det = $(det(Q)) < 0 (a reflection), which " *
            "would flip weld_mesh's normal-vs-gradient winding convention on the mapped-back mesh. " *
            "Negate one column to obtain a rotation.",
        ))
        return Q
    else
        throw(ArgumentError(
            "decompose_3d_surface: projection must be nothing, :random, or a 3x3 orthonormal " *
            "matrix; got $(repr(projection)).",
        ))
    end
end

"""
    _rotate_system(F::System, Q::Matrix{Float64}) -> System

Build the working-chart system `F'(x') := F(Q x')` by simultaneous
substitution, keeping the same variable symbols and order. A chart point `x'`
corresponds to the world point `x = Q x'`.
"""
function _rotate_system(F::System, Q::Matrix{Float64})
    vars = F.variables
    rotated = Q * vars
    exprs = [subs(e, vars => rotated) for e in F.expressions]
    return System(exprs, variables = vars)
end

"""
    _chart_config(cfg::HomotopyConfig{T}, Q::Matrix{Float64}) where {T<:AbstractFloat}
        -> HomotopyConfig{T}

Copy `cfg` with `bbox_*` replaced by the enclosing axis-aligned box of the
rotated world bbox (`Q' * corners`). The chart working domain is therefore
slightly LARGER than the world box at the same cfg -- documented on
[`decompose_3d_surface`](@ref); all other fields are carried over unchanged.
"""
function _chart_config(cfg::HomotopyConfig{T}, Q::Matrix{Float64}) where {T<:AbstractFloat}
    corners = Vector{T}[]
    for bx in cfg.bbox_x, by in cfg.bbox_y, bz in cfg.bbox_z
        push!(corners, T[bx, by, bz])
    end
    QT = Matrix{T}(Q')
    transformed = [QT * c for c in corners]
    lims(k) = (minimum(p[k] for p in transformed), maximum(p[k] for p in transformed))
    return reconstruct(cfg, bbox_x = lims(1), bbox_y = lims(2), bbox_z = lims(3))
end

# Fixed, deterministic, generic-looking complex probe points for the
# polynomial-identity test in _verify_projection_ok: a nonzero polynomial is
# nonzero at a generic point, and an identically-zero one evaluates to exactly
# 0 everywhere, so 3 probes give a reliable, reproducible discriminator with
# no rng involvement.
const _PROJECTION_PROBES = [
    ComplexF64[0.754877666+0.31im, -0.410805062+0.27im, 0.228756604-0.53im],
    ComplexF64[-0.318309886+0.15im, 0.577215665-0.42im, -0.693147181+0.36im],
    ComplexF64[0.267949192-0.24im, -0.866025404+0.19im, 0.414213562+0.41im],
]

"""
    _verify_projection_ok(F_chart::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}

BertiniReal's `verify_projection_ok` analogue: reject projections that are
degenerate for THIS surface before the pipeline runs, so the failure is a
loud `ArgumentError` naming the problem instead of an `OverflowError` from
deep inside HomotopyContinuation's start-system construction.

Checks exactly the crash class: whether either AUGMENTING partial
(`∂f'/∂x`, `∂f'/∂y` -- the equations `compute_critical_points` appends)
vanishes identically, via evaluation at 3 fixed generic complex probe points
with a relative threshold of `cfg.jacobian_rank_tol` against the largest
partial's probe magnitude. Deliberately scoped: `∂f'/∂z ≡ 0` is NOT checked
(a chart-z-independent surface, e.g. a cylinder, has an empty critical system
and decomposes fine), and positive-dimensional critical loci are left to the
downstream machinery that already detects them. Measured discrimination:
`z - x^2` with `Q = I` flagged (probe max exactly 0), a 1e-9 near-degenerate
rotation flagged (ratio 6e-10), a 1e-7 rotation passes, and the sphere /
Taubin heart / cylinder all pass.
"""
function _verify_projection_ok(F_chart::System, cfg::HomotopyConfig{T}) where {T<:AbstractFloat}
    vars = F_chart.variables
    f = F_chart.expressions[1]
    partials = [differentiate(f, v) for v in vars]
    probe_max = [
        maximum(abs(ComplexF64(evaluate(g, vars => p))) for p in _PROJECTION_PROBES)
        for g in partials
    ]
    scale = maximum(probe_max)
    scale > 0 || throw(ArgumentError(
        "decompose_3d_surface: the transformed surface equation has an identically-zero gradient " *
        "(constant f?); nothing to decompose.",
    ))
    tol = Float64(cfg.jacobian_rank_tol)
    for i in 1:2
        if probe_max[i] <= tol * scale
            throw(ArgumentError(
                "decompose_3d_surface: the chosen projection is degenerate for this surface -- " *
                "∂f/∂$(vars[i]) of the transformed equation vanishes identically (probe max " *
                "$(probe_max[i]) vs gradient scale $(scale), relative threshold " *
                "cfg.jacobian_rank_tol = $(tol)). compute_critical_points would receive an " *
                "identically-zero augmenting equation, which HomotopyContinuation cannot build a " *
                "start system for. Use projection = :random or a different matrix.",
            ))
        end
    end
    return nothing
end

"""
    _map_to_world(vertices, edges, faces, mesh, Q::Matrix{Float64})

Map a chart-frame decomposition back to world coordinates via `p = Q p'`:
vertex coordinates (complex), edge sample points, face `mesh_vertices` rows,
and the welded mesh's points (triangles unchanged; det(Q) = +1 preserves the
winding convention). `Face.mid_slice_z`, `Face.boundary_edges`, ids,
singularity flags, and vertex metadata (including the chart-frame
`:fixed_variable`/`:fixed_value` boundary annotations) are carried through
unchanged.
"""
function _map_to_world(
    vertices::Vector{NativeVertex{T}},
    edges::Vector{Edge{T}},
    faces::Vector{Face{T}},
    mesh::GeometryBasics.Mesh,
    Q::Matrix{Float64},
) where {T<:AbstractFloat}
    QT = Matrix{T}(Q)

    world_vertices = NativeVertex{T}[
        NativeVertex{T}(
            id = v.id,
            coordinates = QT * v.coordinates,
            v_type = v.v_type,
            metadata = v.metadata,
        )
        for v in vertices
    ]

    world_edges = Edge{T}[
        Edge{T}(
            id = e.id,
            left_vertex_id = e.left_vertex_id,
            right_vertex_id = e.right_vertex_id,
            sampled_points = [QT * p for p in e.sampled_points],
            is_singular = e.is_singular,
        )
        for e in edges
    ]

    # mesh_vertices rows are points: (P * Q^T) row i == (Q * p_i)^T.
    world_faces = Face{T}[
        Face{T}(
            id = fc.id,
            mid_slice_z = fc.mid_slice_z,
            boundary_edges = fc.boundary_edges,
            mesh_vertices = fc.mesh_vertices * transpose(QT),
            mesh_topology = fc.mesh_topology,
        )
        for fc in faces
    ]

    Q32 = Matrix{Float32}(Q)
    world_points = [
        begin
            w = Q32 * Float32[p[1], p[2], p[3]]
            GeometryBasics.Point3f(w[1], w[2], w[3])
        end
        for p in GeometryBasics.coordinates(mesh)
    ]
    world_mesh = GeometryBasics.Mesh(world_points, collect(GeometryBasics.faces(mesh)))

    return world_vertices, world_edges, world_faces, world_mesh
end

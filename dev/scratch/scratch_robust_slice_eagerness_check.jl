# scratch_robust_slice_eagerness_check.jl
#
# Evidence record for the Phase 7.5 hardening pass (2026-07): WHY
# `SurfaceDecomposition._robust_slice_at_z`'s gradient gate runs EAGERLY
# on the very first (naive-midpoint) attempt, and why the retry-armed
# lazy variant (gradient gate active only from the first retry onward,
# which the docstring originally -- wrongly -- claimed was the
# implemented behavior) was evaluated and rejected.
#
# Construction: Taubin heart with bbox_z = (-0.96, 1.3). The critical z
# at -1.0 falls outside bbox_z, so the bottom slab becomes [-0.96, 1.0]
# and its naive midpoint is z_mid = 0.02 -- exactly the candidate
# test/test_taubin.jl section 6 already proves is topologically CLEAN
# but catastrophically gradient-degenerate. The degenerate plane (z=0)
# is NOT a critical z-value, so nothing pins it to a slab's exact
# center: any asymmetric bbox_z crop is enough to move the naive
# midpoint right next to it.
#
# Measured against the eager implementation (2026-07-21):
#   (a) topology gate on the naive z=0.02 slice: does NOT fire
#       (2 ordinary Critical vertices -- a lazy gate would accept it
#       silently, with no error, warning, or failing residual check
#       anywhere downstream);
#   (b) eager gradient gate: rejects z=0.02, retries land on z ≈ 0.0592;
#   (c) downstream sweep max |f|:  1.60 if z=0.02 were accepted,
#       2.4e-7 at the eager gate's z ≈ 0.0592.
#
# Run from the repo root:  julia --project=. dev/scratch/scratch_robust_slice_eagerness_check.jl

using HomotopyContinuation
using HomotopyGetsReal
const HGR = HomotopyGetsReal

@var x y z
f = (x^2 + (1.2 * y)^2 + z^2 - 1)^3 - x^2 * z^3 - 0.1 * (1.2 * y)^2 * z^3
Fh = System([f], variables = [x, y, z])
cfg = HomotopyConfig{Float64}(
    bbox_x = (-1.5, 1.5), bbox_y = (-1.5, 1.5), bbox_z = (-0.96, 1.3),
    edge_sample_density = 8, midslice_sample_density = 8,
)
patch = HGR.build_patch_system(Fh)

z_crits = compute_critical_z_slices(Fh, cfg)
z_in = filter(zc -> cfg.bbox_z[1] <= zc <= cfg.bbox_z[2], z_crits)
z_bounds = sort(unique(vcat([cfg.bbox_z[1]], z_in, [cfg.bbox_z[2]])))
println("z_crits (all)   = ", round.(sort(z_crits); digits = 4))
println("z_bounds (bbox) = ", round.(z_bounds; digits = 4))

zb, zt = z_bounds[1], z_bounds[2]
z_naive = (zb + zt) / 2
println("bottom slab = [$(zb), $(zt)],  naive z_mid = $(z_naive)")

# (a) topology gate alone on the naive midpoint
v_naive, e_naive = slice_at_z(Fh, z_naive, cfg)
topo_fires = any(
    v -> v.v_type == HGR.Artificial && get(v.metadata, :origin, nothing) == :endpoint_fallback,
    v_naive,
) && any(v -> v.v_type == HGR.Singular, v_naive)
println()
println("(a) naive midpoint slice: $(length(v_naive)) vertices ($(join(unique(string(v.v_type) for v in v_naive), ", "))), $(length(e_naive)) edges")
println("    topology gate fires on naive midpoint: $(topo_fires)")
println("    => a retry-armed lazy gate would $(topo_fires ? "still retry" : "ACCEPT z_mid=$(z_naive) with no gradient check")")

# (b) current eager behavior
t_rob = @elapsed (v_r, e_r, z_used) = HGR._robust_slice_at_z(Fh, patch, zb, zt, cfg)
println()
println("(b) eager _robust_slice_at_z: z_mid_used = $(z_used)  (naive was $(z_naive))  [$(round(t_rob; digits = 2)) s]")
println("    gradient gate rejected the naive midpoint: $(!isapprox(z_used, z_naive; atol = 1e-12))")

# (c) downstream damage if z=0.02 were accepted (single edge sweep suffices)
if !isempty(e_naive)
    face_naive = HGR.track_face(Fh, patch, e_naive[1], z_naive, zb, zt, 1, cfg)
    resid_naive = maximum(
        abs(HGR._residual_at(patch, face_naive.mesh_vertices[i, 1], face_naive.mesh_vertices[i, 2], face_naive.mesh_vertices[i, 3], cfg))
        for i in 1:size(face_naive.mesh_vertices, 1)
    )
    face_used = HGR.track_face(Fh, patch, e_r[1], z_used, zb, zt, 2, cfg)
    resid_used = maximum(
        abs(HGR._residual_at(patch, face_used.mesh_vertices[i, 1], face_used.mesh_vertices[i, 2], face_used.mesh_vertices[i, 3], cfg))
        for i in 1:size(face_used.mesh_vertices, 1)
    )
    println()
    println("(c) sweep max|f| if naive z=$(z_naive) accepted : $(resid_naive)")
    println("    sweep max|f| at eager-gate's  z=$(z_used) : $(resid_used)")
end

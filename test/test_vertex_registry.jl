@testset "VertexRegistry (Phase 10, Stage 1)" begin

println("=" ^ 70)
println("1. Empty registry, first registration, and non-matching second registration")
println("=" ^ 70)

reg = VertexRegistry{Float64}(1e-3)
@test length(reg) == 0

v1 = NativeVertex{Float64}(id = 7, coordinates = [1.0 + 0im, 2.0 + 0im], v_type = Critical,
    metadata = Dict{Symbol,Any}(:tag => :first))
id1 = register!(reg, v1)
@test id1 == 1                                   # registry-assigned id, NOT candidate.id (7)
@test length(reg) == 1
@test reg.vertices[1].coordinates == v1.coordinates
@test reg.vertices[1].v_type == Critical
@test reg.vertices[1].metadata[:tag] == :first    # sole member: :tag IS the (trivially) shared key
@test reg.vertices[1].metadata[:cluster_size] == 1   # uniform bookkeeping even pre-merge
@test reg.vertices[1].metadata[:cluster_member_ids] == [7]
println("first registration: id=$id1 (candidate carried id=7 -- registry ids are independent)")

v2_far = NativeVertex{Float64}(id = 1, coordinates = [50.0 + 0im, 50.0 + 0im], v_type = Boundary)
id2 = register!(reg, v2_far)
@test id2 == 2
@test length(reg) == 2
@test reg.vertices[1].coordinates == v1.coordinates   # first entry untouched by the second registration
println("far registration (beyond tol): id=$id2, new entry, first entry unchanged.")

println()
println("=" ^ 70)
println("2. Matching registration: same id returned, no new entry, coordinates averaged")
println("=" ^ 70)

v3_near = NativeVertex{Float64}(id = 99, coordinates = [1.0002 + 0im, 2.0001 + 0im], v_type = Critical)
id3 = register!(reg, v3_near)
@test id3 == id1                                 # merged into the FIRST entry, same id returned
@test length(reg) == 2                           # no new canonical entry created
expected_centroid = (v1.coordinates .+ v3_near.coordinates) ./ 2
@test isapprox(real.(reg.vertices[1].coordinates), real.(expected_centroid); atol = 1e-12)
@test reg.vertices[1].metadata[:cluster_size] == 2
@test sort(reg.vertices[1].metadata[:cluster_member_ids]) == [7, 99]
println("matching registration: id=$id3 == $id1 (merged, not new), centroid updated, cluster_size=2.")

println()
println("=" ^ 70)
println("3. THIRD incremental merge into the same entry: running centroid + cluster_size,")
println("   confirming the full-recompute design (not a lossy two-at-a-time merge)")
println("=" ^ 70)

v4_near = NativeVertex{Float64}(id = 42, coordinates = [0.9998 + 0im, 1.9999 + 0im], v_type = Critical)
id4 = register!(reg, v4_near)
@test id4 == id1
@test length(reg) == 2
expected_centroid3 = (v1.coordinates .+ v3_near.coordinates .+ v4_near.coordinates) ./ 3
@test isapprox(real.(reg.vertices[1].coordinates), real.(expected_centroid3); atol = 1e-12)
@test reg.vertices[1].metadata[:cluster_size] == 3
@test sort(reg.vertices[1].metadata[:cluster_member_ids]) == [7, 42, 99]
println("third merge: cluster_size=3, centroid = true 3-way average (not a stale 2-item merge).")

println()
println("=" ^ 70)
println("4. Type-resolution rule (mirrors cluster_vertices' Singular > common > Artificial)")
println("=" ^ 70)

reg_t = VertexRegistry{Float64}(1e-3)
a = NativeVertex{Float64}(id = 1, coordinates = [0.0 + 0im], v_type = Critical)
b = NativeVertex{Float64}(id = 2, coordinates = [0.0001 + 0im], v_type = Singular)
ida = register!(reg_t, a)
idb = register!(reg_t, b)
@test idb == ida
@test reg_t.vertices[ida].v_type == Singular
println("Critical + Singular -> Singular (matches cluster_vertices' Singular-wins rule).")

c = NativeVertex{Float64}(id = 3, coordinates = [0.0002 + 0im], v_type = Critical)
idc = register!(reg_t, c)
@test idc == ida
@test reg_t.vertices[ida].v_type == Singular      # Singular-wins PERSISTS across a further merge
println("further Critical merge does not revert the entry from Singular -- persists correctly.")

reg_t2 = VertexRegistry{Float64}(1e-3)
d1 = NativeVertex{Float64}(id = 1, coordinates = [0.0 + 0im], v_type = Boundary)
d2 = NativeVertex{Float64}(id = 2, coordinates = [0.0001 + 0im], v_type = Boundary)
register!(reg_t2, d1)
idd = register!(reg_t2, d2)
@test reg_t2.vertices[idd].v_type == Boundary     # common type preserved, not forced to Artificial
println("common non-Singular type (Boundary + Boundary) preserved, not forced to Artificial.")

reg_t3 = VertexRegistry{Float64}(1e-3)
e1 = NativeVertex{Float64}(id = 1, coordinates = [0.0 + 0im], v_type = Critical)
e2 = NativeVertex{Float64}(id = 2, coordinates = [0.0001 + 0im], v_type = Boundary)
register!(reg_t3, e1)
ide = register!(reg_t3, e2)
@test reg_t3.vertices[ide].v_type == Artificial   # disagreeing, neither Singular -> Artificial
println("disagreeing non-Singular types (Critical + Boundary) -> Artificial, the third branch.")

println()
println("=" ^ 70)
println("5. Metadata combination delegates to _combine_values (mean for floats, min for jacobian_rank)")
println("=" ^ 70)

reg_m = VertexRegistry{Float64}(1e-3)
f1 = NativeVertex{Float64}(id = 1, coordinates = [0.0 + 0im], v_type = Critical,
    metadata = Dict{Symbol,Any}(:tolerance_used => 1e-6))
f2 = NativeVertex{Float64}(id = 2, coordinates = [0.0001 + 0im], v_type = Critical,
    metadata = Dict{Symbol,Any}(:tolerance_used => 2e-6))
register!(reg_m, f1)
idf = register!(reg_m, f2)
@test isapprox(reg_m.vertices[idf].metadata[:tolerance_used], 1.5e-6; atol = 1e-18)
println("shared Float64 key (:tolerance_used) combined via arithmetic mean (delegated to _combine_values).")

println()
println("=" ^ 70)
println("6. Cross-check against test_solver.jl's own cluster_vertices fixture:")
println("   incremental register! must reach the SAME domain answer as batch cluster_vertices")
println("   on the identical near-dupe/jacobian_rank vertices (registered ONE AT A TIME here).")
println("=" ^ 70)

reg_x = VertexRegistry{Float64}(1e-3)
nd1 = NativeVertex{Float64}(id = 1, coordinates = [1.0 + 0im, 1.0 + 0im], v_type = Critical,
    metadata = Dict{Symbol,Any}(:jacobian_rank => 2))
nd2 = NativeVertex{Float64}(id = 2, coordinates = [1.00001 + 0im, 0.99999 + 0im], v_type = Critical,
    metadata = Dict{Symbol,Any}(:jacobian_rank => 2))
nd3 = NativeVertex{Float64}(id = 3, coordinates = [1.00002 + 0im, 1.00001 + 0im], v_type = Singular,
    metadata = Dict{Symbol,Any}(:jacobian_rank => 1))
nd4 = NativeVertex{Float64}(id = 4, coordinates = [5.0 + 0im, 5.0 + 0im], v_type = Critical,
    metadata = Dict{Symbol,Any}(:jacobian_rank => 2))
idx1 = register!(reg_x, nd1)
idx2 = register!(reg_x, nd2)
idx3 = register!(reg_x, nd3)
idx4 = register!(reg_x, nd4)
@test idx1 == idx2 == idx3               # 1,2,3 merge (same as cluster_vertices' clustered_loose)
@test idx4 != idx1                       # 4 stays separate
@test length(reg_x) == 2                 # matches cluster_vertices([nd1..nd4], 1e-3) -> length 2
merged_entry = reg_x.vertices[idx1]
@test merged_entry.v_type == Singular
@test merged_entry.metadata[:jacobian_rank] == 1     # MIN(2,2,1), not averaged to 1.667
@test merged_entry.metadata[:jacobian_rank] isa Integer
separate_entry = reg_x.vertices[idx4]
@test separate_entry.v_type == Critical
@test separate_entry.metadata[:jacobian_rank] == 2
println("register!, called one-at-a-time, reproduces cluster_vertices' batch answer exactly:")
println("  2 canonical entries; merged entry Singular with jacobian_rank=min(2,2,1)=1 (not averaged).")

println()
println("=" ^ 70)
println("7. Nearest-match determinism: candidate merges into the CLOSEST entry, not first/last")
println("=" ^ 70)

reg_n = VertexRegistry{Float64}(0.05)
n1 = NativeVertex{Float64}(id = 1, coordinates = [0.0 + 0im], v_type = Critical)
n2 = NativeVertex{Float64}(id = 2, coordinates = [1.0 + 0im], v_type = Critical)
n3 = NativeVertex{Float64}(id = 3, coordinates = [2.0 + 0im], v_type = Critical)
i1 = register!(reg_n, n1)
i2 = register!(reg_n, n2)
i3 = register!(reg_n, n3)
@test length(reg_n) == 3
cand = NativeVertex{Float64}(id = 4, coordinates = [1.02 + 0im], v_type = Critical)  # closest to n2
idn = register!(reg_n, cand)
@test idn == i2
@test length(reg_n) == 3
println("candidate at 1.02 merges into the entry at 1.0 (nearest), not 0.0 or 2.0.")

println()
println("=" ^ 70)
println("8. T-genericity: BigFloat")
println("=" ^ 70)

reg_bf = VertexRegistry{BigFloat}(BigFloat(1e-4))
@test reg_bf isa VertexRegistry{BigFloat}
bf1 = NativeVertex{BigFloat}(id = 1, coordinates = [BigFloat(1.0) + 0im], v_type = Critical)
bf2 = NativeVertex{BigFloat}(id = 2, coordinates = [BigFloat(1.00001) + 0im], v_type = Critical)
idbf1 = register!(reg_bf, bf1)
idbf2 = register!(reg_bf, bf2)
@test idbf1 == idbf2
@test length(reg_bf) == 1
@test eltype(reg_bf.vertices[idbf1].coordinates) == Complex{BigFloat}
println("BigFloat registry: merges correctly, coordinates stay Complex{BigFloat}.")

println()
println("=" ^ 70)
println("9. Degenerate tol = 0: only bit-identical coordinates merge")
println("=" ^ 70)

reg_z = VertexRegistry{Float64}(0.0)
z1 = NativeVertex{Float64}(id = 1, coordinates = [3.0 + 0im], v_type = Critical)
z2 = NativeVertex{Float64}(id = 2, coordinates = [3.0 + 0im], v_type = Critical)       # exact duplicate
z3 = NativeVertex{Float64}(id = 3, coordinates = [3.0 + 1e-15im], v_type = Critical)   # NOT bit-identical
idz1 = register!(reg_z, z1)
idz2 = register!(reg_z, z2)
idz3 = register!(reg_z, z3)
@test idz1 == idz2
@test idz3 != idz1
@test length(reg_z) == 2
println("tol=0: exact duplicate merges, epsilon-different coordinate does not.")

end

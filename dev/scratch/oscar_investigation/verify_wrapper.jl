# Step 2: verify examples/oscar_integration.jl's oscar_ideal_to_system
# against the exact two repro cases confirmed broken in Task 3.
# (End-to-end decompose_3d_surface/decompose_1d_curve smoke test skipped
# this round -- blocked by the pre-existing display-sleep GLMakie
# segfault; Juan explicitly chose to accept the two repro-case fixes as
# sufficient verification for this round rather than wait.)

using Oscar
using OscarHomotopyContinuation
using HomotopyContinuation

include("/Users/juancagc/HomotopyGetsReal/examples/oscar_integration.jl")

println("="^70)
println("Repro case 1: astroid, ring order [z_var, x_var, y_var], z_var unused")
println("(Task 3 original finding: bridge's bare System(I) DROPPED z_var")
println("entirely -- sys.variables came back as only [x_var, y_var])")
println("="^70)

R, (z_var, x_var, y_var) = Oscar.polynomial_ring(Oscar.QQ, ["z_var", "x_var", "y_var"])
f_astroid = (x_var^2 + y_var^2 - 1)^3 + 27 * x_var^2 * y_var^2
I = Oscar.ideal(R, [f_astroid])

sys_bare = OscarHomotopyContinuation.System(I)
println("Bare bridge System(I).variables (for comparison, still buggy): ", sys_bare.variables)

sys_wrapped = oscar_ideal_to_system(I)
println("oscar_ideal_to_system(I).variables:                             ", sys_wrapped.variables)
println("gens(base_ring(I)):                                             ", Oscar.gens(Oscar.base_ring(I)))

test1_fixed = length(sys_wrapped.variables) == 3 && string.(sys_wrapped.variables) == string.(Oscar.gens(Oscar.base_ring(I)))
println("\nFIXED (3 variables, exact gens() order)? ", test1_fixed)

println("\n" * "="^70)
println("Repro case 2: ellipsoid, ring order [c_var, a_var, b_var]")
println("(Task 3 original finding: bridge's bare System(I2) silently")
println("REORDERED to alphabetical [a_var, b_var, c_var])")
println("="^70)

R2, (c_var, a_var, b_var) = Oscar.polynomial_ring(Oscar.QQ, ["c_var", "a_var", "b_var"])
f_ellipsoid = a_var^2 + 4 * b_var^2 + 9 * c_var^2 - 1
I2 = Oscar.ideal(R2, [f_ellipsoid])

sys2_bare = OscarHomotopyContinuation.System(I2)
println("Bare bridge System(I2).variables (for comparison, still buggy): ", sys2_bare.variables)

sys2_wrapped = oscar_ideal_to_system(I2)
println("oscar_ideal_to_system(I2).variables:                             ", sys2_wrapped.variables)
println("gens(base_ring(I2)):                                             ", Oscar.gens(Oscar.base_ring(I2)))

test2_fixed = string.(sys2_wrapped.variables) == string.(Oscar.gens(Oscar.base_ring(I2)))
println("\nFIXED (exact gens() order [c_var, a_var, b_var])? ", test2_fixed)

println("\n" * "="^70)
println("HEADLINE")
println("="^70)
println("Test 1 (astroid, dropping bug):    FIXED = ", test1_fixed)
println("Test 2 (ellipsoid, reordering bug): FIXED = ", test2_fixed)
if test1_fixed && test2_fixed
    println("\nBoth confirmed bug cases from Task 3 are fixed by oscar_ideal_to_system.")
else
    println("\nAt least one case is STILL BROKEN -- do not treat the wrapper as verified.")
end
println("\nNOTE: end-to-end decompose_3d_surface/decompose_1d_curve smoke test")
println("NOT run this round (GLMakie precompile blocked by display sleep;")
println("skipped per explicit decision) -- outstanding for a future check.")
println("\nVERIFY DONE")

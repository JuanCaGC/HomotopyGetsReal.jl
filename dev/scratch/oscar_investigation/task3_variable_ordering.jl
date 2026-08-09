# Task 3: empirical variable-ordering test. Does NOT load HomotopyGetsReal
# (sidesteps the environmental GLMakie precompile issue) -- only needs
# Oscar, the bridge, and raw HomotopyContinuation to inspect the
# resulting System's .variables field order.

using Oscar
using OscarHomotopyContinuation
using HomotopyContinuation

println("="^70)
println("Test 1: astroid, ring order [z_var, x_var, y_var] (z FIRST -- stress")
println("test against decompose_3d_surface's own \"z LAST\" convention)")
println("="^70)

R, (z_var, x_var, y_var) = Oscar.polynomial_ring(Oscar.QQ, ["z_var", "x_var", "y_var"])
println("gens(R) order: ", Oscar.gens(R))

# Astroid, only 2 of the 3 ring variables actually appear in the equation
# (z_var is a genuine extra unused ring variable -- also tests whether the
# bridge handles a variable with zero appearances in the polynomial).
f_astroid = (x_var^2 + y_var^2 - 1)^3 + 27 * x_var^2 * y_var^2
I = Oscar.ideal(R, [f_astroid])

sys = OscarHomotopyContinuation.System(I)
println("\nBridge-built HomotopyContinuation.System(I).variables: ", sys.variables)
println("gens(parent(f_astroid)) order:                          ", Oscar.gens(Oscar.parent(f_astroid)))

matches_1 = string.(sys.variables) == string.(Oscar.gens(R))
println("\nDoes sys.variables match gens(R) order EXACTLY? ", matches_1)

println("\n" * "="^70)
println("Test 2: ellipsoid, ring order [c_var, a_var, b_var] (arbitrary,")
println("non-alphabetical, 3 variables all actually used -- second")
println("independent check with a genuinely 3-variable surface)")
println("="^70)

R2, (c_var, a_var, b_var) = Oscar.polynomial_ring(Oscar.QQ, ["c_var", "a_var", "b_var"])
println("gens(R2) order: ", Oscar.gens(R2))

# Ellipsoid x^2+4y^2+9z^2-1, relabeled a_var<->x, b_var<->y, c_var<->z to
# keep the SAME already-validated HGR fixture/equation shape while using
# non-x/y/z names (so this is a genuine ring-order test, not just an
# accidental x/y/z alphabetical coincidence).
f_ellipsoid = a_var^2 + 4 * b_var^2 + 9 * c_var^2 - 1
I2 = Oscar.ideal(R2, [f_ellipsoid])

sys2 = OscarHomotopyContinuation.System(I2)
println("\nBridge-built HomotopyContinuation.System(I2).variables: ", sys2.variables)
println("gens(parent(f_ellipsoid)) order:                          ", Oscar.gens(Oscar.parent(f_ellipsoid)))

matches_2 = string.(sys2.variables) == string.(Oscar.gens(R2))
println("\nDoes sys2.variables match gens(R2) order EXACTLY? ", matches_2)

println("\n" * "="^70)
println("HEADLINE RESULT")
println("="^70)
println("Test 1 (astroid, z_var first, 1 unused variable):  matches gens() order = ", matches_1)
println("Test 2 (ellipsoid, arbitrary 3-var order):          matches gens() order = ", matches_2)
if matches_1 && matches_2
    println("\nBridge's System(I) DOES preserve gens(parent(f)) order in both tests.")
else
    println("\nBridge's System(I) does NOT reliably preserve gens(parent(f)) order --")
    println("this is the silent-reordering bug, confirmed live.")
end
println("\nTASK 3 DONE")

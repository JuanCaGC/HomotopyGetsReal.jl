# Task 1: satellite environment setup for the OSCAR compatibility investigation.
# Times Pkg.add("Oscar"), Pkg.add(url=bridge), and first `using Oscar` +
# bridge precompile, since Oscar.jl's precompile time is a real, not
# assumed, input to the Task 4 architecture recommendation.

using Pkg
Pkg.activate("/Users/juancagc/HomotopyGetsReal/dev/scratch/oscar_investigation")

println("=== Pkg.develop(HGR) ===")
t0 = time()
Pkg.develop(path = "/Users/juancagc/HomotopyGetsReal")
println("Pkg.develop(HGR): $(round(time()-t0, digits=1))s")

println("\n=== Pkg.add(\"Oscar\") ===")
t0 = time()
Pkg.add("Oscar")
t_add_oscar = time() - t0
println("Pkg.add(Oscar): $(round(t_add_oscar, digits=1))s")

println("\n=== Pkg.add(url=OscarHomotopyContinuation bridge) ===")
t0 = time()
try
    Pkg.add(url = "https://github.com/taboege/OscarHomotopyContinuation")
    global t_add_bridge = time() - t0
    println("Pkg.add(url=bridge): $(round(t_add_bridge, digits=1))s -- SUCCEEDED")
catch e
    global t_add_bridge = time() - t0
    println("Pkg.add(url=bridge): FAILED after $(round(t_add_bridge, digits=1))s")
    println(sprint(showerror, e))
    rethrow(e)
end

println("\n=== First `using Oscar` (precompile timing) ===")
t0 = time()
@eval using Oscar
t_using_oscar = time() - t0
println("using Oscar (first, incl. precompile): $(round(t_using_oscar, digits=1))s")

println("\n=== `using OscarHomotopyContinuation` ===")
t0 = time()
@eval using OscarHomotopyContinuation
t_using_bridge = time() - t0
println("using OscarHomotopyContinuation: $(round(t_using_bridge, digits=1))s")

println("\n=== `using HomotopyGetsReal` (alongside Oscar) ===")
t0 = time()
@eval using HomotopyGetsReal
t_using_hgr = time() - t0
println("using HomotopyGetsReal: $(round(t_using_hgr, digits=1))s")

println("\n=== TIMING SUMMARY ===")
println("Pkg.add(Oscar):              $(round(t_add_oscar, digits=1))s")
println("Pkg.add(url=bridge):         $(round(t_add_bridge, digits=1))s")
println("using Oscar (precompile):    $(round(t_using_oscar, digits=1))s")
println("using OscarHomotopyContinuation: $(round(t_using_bridge, digits=1))s")
println("using HomotopyGetsReal:      $(round(t_using_hgr, digits=1))s")
println("TOTAL setup wall time:       $(round(t_add_oscar+t_add_bridge+t_using_oscar+t_using_bridge+t_using_hgr, digits=1))s")
println("\nSETUP DONE")

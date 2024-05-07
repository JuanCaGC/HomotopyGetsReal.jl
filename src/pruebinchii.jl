include("plotbr.jl")
using .plotbrModule
include("curve.jl")
using .curveModule
include("parser.jl")
using .parserModule


path = "/home/juancagc/.julia/dev/BertiniReal/pruebas/astroid/Dir_Name"
directory, MPtype, dimension = parse_directory_name(path)

cur = Curve(directory)

plotbr(cur)
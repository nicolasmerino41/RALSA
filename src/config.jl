const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW = joinpath(ROOT, "data", "raw")
const INPUT = joinpath(ROOT, "data", "processed")
const OUT = joinpath(ROOT, "outputs")
const DATA = ["Quercus", "Nahuel", "Salix_Galpar", "Gottin_HP", "Gottin_PP",
              "Garraf_HP", "Garraf_PP", "Garraf_PP2", "Olot", "Montseny"]
const REMOVAL = collect(0.0:0.1:0.8)
const LEVELS = REMOVAL
const NREPS = 500
const SEED = 4300
const EXAMPLE = "Salix_Galpar"
const GROUPS = ["Specialists", "Generalists"]
const DISPLAY = Dict(s => replace(s, "_" => " ") for s in DATA)
DISPLAY["Salix_Galpar"] = "Salix–Galpar"

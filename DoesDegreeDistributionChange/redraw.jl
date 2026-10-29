using CSV,DataFrames,Statistics,CairoMakie,Printf
const OUT=joinpath(@__DIR__,"outputs")
const B=2000
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
include("plots_and_report.jl")

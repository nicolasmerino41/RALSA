using CSV,DataFrames,Statistics,CairoMakie,Printf
dest=abspath(ARGS[1])
datasets=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
species=CSV.read(joinpath(dest,"species.csv"),DataFrame);links=CSV.read(joinpath(dest,"links.csv"),DataFrame)
stats=CSV.read(joinpath(dest,"species_associations.csv"),DataFrame);groups=CSV.read(joinpath(dest,"degree_groups.csv"),DataFrame)
origins=CSV.read(joinpath(dest,"rare_link_origins.csv"),DataFrame);sensitivity=CSV.read(joinpath(dest,"specialist_thresholds.csv"),DataFrame)
td=CSV.read(joinpath(dest,"conditional_tests.csv"),DataFrame)
function ranks(x)
 p=sortperm(x);r=zeros(length(x));i=1
 while i<=length(x)
  j=i;while j<length(x) && x[p[j+1]]==x[p[i]];j+=1;end
  r[p[i:j]].=(i+j)/2;i=j+1
 end
 r
end
rho(x,y)=length(x)>2 && std(x)>0 && std(y)>0 ? cor(ranks(x),ranks(y)) : NaN
include("plots.jl")
include("report.jl")

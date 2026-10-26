using CSV,DataFrames,CairoMakie,Statistics,Printf,SHA
root=ARGS[1];dest=ARGS[2];mkpath(dest)
source=joinpath(root,"outputs/statistical description of local support/opportunity_theory/conditioned_M/curves.csv")
d=CSV.read(source,DataFrame);d=d[d.model.=="Common",:]
datasets=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
fig=Figure(size=(1800,900));summary=NamedTuple[];curves=NamedTuple[]
for (i,ds) in enumerate(datasets)
 v=d[d.dataset.==ds,:];sort!(v,:m);x=100v.removal_fraction;local_retention=100 .- x;regional_retention=local_retention.+v.empirical_gap_pp
 @assert isapprox(first(regional_retention),100;atol=1e-8) && isapprox(last(regional_retention),0;atol=1e-8)
 @assert all(regional_retention.>=-1e-8) && all(regional_retention.>=local_retention.-1e-8)
 gap=mean(v.empirical_gap_pp[v.removal_fraction.<=.8]);push!(summary,(dataset=ds,mean_gap_0_80_pp=gap))
 ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" ")*"\n"*@sprintf("Mean gap (0–80%%): %.1f pp",gap),titlesize=18,xlabel="Sites removed (%)",ylabel=mod1(i,5)==1 ? "Retained (% of original total)" : "",xticks=0:20:100,yticks=0:20:100)
 band!(ax,x,local_retention,regional_retention;color=("#CC79A7",.20))
 lines!(ax,x,regional_retention;color="#0072B2",linewidth=3)
 lines!(ax,x,local_retention;color="#D55E00",linewidth=3)
 xlims!(ax,0,100);ylims!(ax,0,100)
 for j in eachindex(x);push!(curves,(dataset=ds,m=v.m[j],removed_percent=x[j],regional_retention_percent=regional_retention[j],local_retention_percent=local_retention[j],gap_pp=v.empirical_gap_pp[j]));end
end
Legend(fig[3,1:5],[LineElement(color="#0072B2",linewidth=3),LineElement(color="#D55E00",linewidth=3),PolyElement(color=("#CC79A7",.20))],["Regional links retained","Local interaction occurrences retained","Regional–local gap"],orientation=:horizontal,framevisible=false)
save(joinpath(dest,"01_regional_local_retention.png"),fig;px_per_unit=2)
CSV.write(joinpath(dest,"mean_gaps.csv"),DataFrame(summary));CSV.write(joinpath(dest,"curves.csv"),DataFrame(curves))


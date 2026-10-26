using CSV,DataFrames,CairoMakie,Statistics,Printf,SHA
root=ARGS[1];dest=ARGS[2];mkpath(dest)
source=joinpath(root,"outputs/statistical description of local support/opportunity_theory/conditioned_M/curves.csv")
d=CSV.read(source,DataFrame);d=d[d.model.=="Common",:]
datasets=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
fig=Figure(size=(1800,900));summary=NamedTuple[];curves=NamedTuple[]
for (i,ds) in enumerate(datasets)
 v=d[d.dataset.==ds,:];sort!(v,:m);x=100v.removal_fraction;local_loss=x;regional_loss=x.-v.empirical_gap_pp
 @assert isapprox(first(regional_loss),0;atol=1e-8) && isapprox(last(regional_loss),100;atol=1e-8)
 @assert all(regional_loss.>=-1e-8) && all(regional_loss.<=local_loss.+1e-8)
 gap=mean(v.empirical_gap_pp[v.removal_fraction.<=.8]);push!(summary,(dataset=ds,mean_gap_0_80_pp=gap))
 ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" ")*"\n"*@sprintf("Mean gap (0–80%%): %.1f pp",gap),titlesize=18,xlabel="Sites removed (%)",ylabel=mod1(i,5)==1 ? "Loss (% of original total)" : "",xticks=0:20:100,yticks=0:20:100)
 band!(ax,x,regional_loss,local_loss;color=("#CC79A7",.20))
 lines!(ax,x,regional_loss;color="#0072B2",linewidth=3)
 lines!(ax,x,local_loss;color="#D55E00",linewidth=3)
 xlims!(ax,0,100);ylims!(ax,0,100)
 for j in eachindex(x);push!(curves,(dataset=ds,m=v.m[j],removed_percent=x[j],regional_loss_percent=regional_loss[j],local_loss_percent=local_loss[j],gap_pp=v.empirical_gap_pp[j]));end
end
Legend(fig[3,1:5],[LineElement(color="#0072B2",linewidth=3),LineElement(color="#D55E00",linewidth=3),PolyElement(color=("#CC79A7",.20))],["Regional link loss","Local interaction-occurrence loss","Regional–local gap"],orientation=:horizontal,framevisible=false)
save(joinpath(dest,"01_regional_local_loss.png"),fig;px_per_unit=2)
CSV.write(joinpath(dest,"mean_gaps.csv"),DataFrame(summary));CSV.write(joinpath(dest,"curves.csv"),DataFrame(curves))
write(joinpath(dest,"README.md"),"# Regional and local loss\n\nExact expected loss under uniform removal of a fixed number of sites, using the saved empirical gap curves (the empirical column is repeated across models; Common rows are used once). Regional loss is relative to initial unique links L; local occurrence loss is relative to initial total support M. Both increase from 0 to 100%. Shading spans their difference. Panel annotations use the established arithmetic mean of gap over integer removal levels m/N≤0.8, including zero, in percentage points. This is the established mean-gap summary, not a continuous integrated area. Curves cover 0–100%. No simulations or model refits.\n\nSource: "*relpath(source,root)*"\nSHA256: "*bytes2hex(sha256(read(source)))*"\n")

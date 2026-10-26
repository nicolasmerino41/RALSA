using CSV,DataFrames,CairoMakie,SHA
root=ARGS[1];dest=ARGS[2];mkpath(dest)
source=joinpath(root,"outputs/statistical description of local support/opportunity_theory/opportunity_support_organisation/distributions.csv")
d=CSV.read(source,DataFrame);d=d[d.quantity.=="Opportunity",:]
datasets=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
set_theme!(Theme(font="Arial",fontsize=20,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=17,yticklabelsize=17)))
fig=Figure(size=(1800,850))
for (i,ds) in enumerate(datasets)
 v=d[d.dataset.==ds,:];sort!(v,:K);@assert isapprox(sum(v.frequency),1;atol=1e-8)
 maxn=maximum(v.K[v.frequency.>0]);v=v[v.K.<=maxn,:]
 ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Co-occurring sites, n",ylabel=mod1(i,5)==1 ? "Proportion of pairs" : "")
 barplot!(ax,v.K,v.frequency;color="#D55E00",width=.85)
 xlims!(ax,.4,maxn+.6);ylims!(ax,0,.6)
end
save(joinpath(dest,"01_CO_frequency.png"),fig;px_per_unit=2)
CSV.write(joinpath(dest,"frequencies.csv"),select(d,:dataset,:K=>:n,:frequency=>:proportion))
write(joinpath(dest,"README.md"),"# Co-occurrence frequency\n\nFor the same known regional interacting pairs (K>0) as the previous comparison, n counts sites with both species present. Each bar is the fraction of these pairs with exactly n co-occurring sites. Pairs never observed interacting are excluded. All panels use a common y scale. No model fitting. Confirmed-occurrence inputs and existing quality flags are inherited from opportunity_support_organisation.\nSource: "*relpath(source,root)*"\nSHA256: "*bytes2hex(sha256(read(source)))*"\n")

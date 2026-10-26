using CSV,DataFrames,CairoMakie,SHA
root=ARGS[1];dest=ARGS[2];mkpath(dest)
source=joinpath(root,"outputs/support_prediction_zero/confirmed/pairs.csv")
d=CSV.read(source,DataFrame);d=d[d.n.>=1,:]
datasets=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
set_theme!(Theme(font="Arial",fontsize=20,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=17,yticklabelsize=17)))
fig=Figure(size=(1800,850));rows=NamedTuple[]
for (i,ds) in enumerate(datasets)
 v=d[d.dataset.==ds,:];@assert all(0 .<=v.K .<=v.n .<=v.N)
 x=collect(0:maximum(v.K));p=[count(==(k),v.K)/nrow(v) for k in x];@assert isapprox(sum(p),1)
 ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Local interaction support, K",ylabel=mod1(i,5)==1 ? "Proportion of co-occurring pairs" : "")
 barplot!(ax,x,p;color=[k==0 ? "#999999" : "#0072B2" for k in x],width=.85)
 xlims!(ax,-.6,maximum(x)+.6);ylims!(ax,0,1)
 step=maximum(x)<=10 ? 2 : 5;ax.xticks=0:step:maximum(x)
 for (k,f) in zip(x,p);push!(rows,(dataset=ds,K=k,count=count(==(k),v.K),proportion=f,total_cooccurring_pairs=nrow(v)));end
end
Legend(fig[3,1:5],[PolyElement(color="#999999"),PolyElement(color="#0072B2")],["Never observed interacting (K = 0)","Observed interactions (K ≥ 1)"],orientation=:horizontal,framevisible=false)
save(joinpath(dest,"01_K_frequency_including_zero.png"),fig;px_per_unit=2)
CSV.write(joinpath(dest,"frequencies.csv"),DataFrame(rows))
write(joinpath(dest,"README.md"),"# Local support including zeros\n\nAll pairs with at least one co-occurrence (n≥1) in the existing confirmed-occurrence pairs table. K=0 means no observed interaction anywhere in the dataset, not proven biological impossibility. Each dataset denominator is all co-occurring pairs, including zeros. Bars sum to one. Previous confirmed-occurrence corrections and unresolved occupancy flags remain unchanged.\nSource: "*relpath(source,root)*"\nSHA256: "*bytes2hex(sha256(read(source)))*"\n")

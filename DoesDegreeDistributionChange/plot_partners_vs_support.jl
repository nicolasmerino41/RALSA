using CSV, DataFrames, CairoMakie, SHA
root=ARGS[1]
out=length(ARGS)>1 ? ARGS[2] : joinpath(@__DIR__,"partners_vs_support")
mkpath(out)
input=joinpath(root,"outputs/support_prediction_zero/confirmed/pairs.csv")
p=CSV.read(input,DataFrame);p=p[p.K.>0,:]
datasets=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
colours=["#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#332288","#882255","#44AA99","#666666"]
rows=NamedTuple[]
for ds in datasets
 v=p[p.dataset.==ds,:]
 @assert nrow(unique(v[:,[:consumer,:resource]]))==nrow(v)
 for role in [:consumer,:resource]
  counts=Dict(r[role]=>r.degree for r in eachrow(combine(groupby(v,role),nrow=>:degree)))
  for r in eachrow(v)
   push!(rows,(dataset=ds,role=String(role),species=r[role],consumer=r.consumer,resource=r.resource,regional_partners=counts[r[role]],K=r.K))
  end
 end
end
t=DataFrame(rows)
CSV.write(joinpath(out,"points.csv"),t)
set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false)))
f=Figure(size=(1650,620))
for (j,role) in enumerate(["consumer","resource","both"])
 ax=Axis(f[1,j],title=["Consumers","Resources","Both combined"][j],xlabel="Number of regional partners",ylabel=j==1 ? "Sites supporting the interaction (K)" : "")
 for (ds,c) in zip(datasets,colours)
  v=t[(t.dataset.==ds).&((t.role.==role).|(role=="both")),:]
  scatter!(ax,v.regional_partners,v.K;color=(c,.35),markersize=5)
 end
 xlims!(ax,0,maximum(t.regional_partners)*1.03);ylims!(ax,0,maximum(t.K)*1.03)
end
Legend(f[2,1:3],[MarkerElement(color=c,marker=:circle) for c in colours],replace.(datasets,"_"=>" "),orientation=:horizontal,nbanks=2,framevisible=false,labelsize=16)
save(joinpath(out,"regional_partners_vs_link_support.png"),f;px_per_unit=2)
open(joinpath(out,"README.md"),"w") do io
 println(io,"# Regional partners versus interaction support\n\nEach point is one species–interaction combination. X is the species' number of regional partners; Y is the number of distinct sites supporting that interaction, K. All ten datasets are shown, coloured by dataset. Axes are linear and identical across panels; points are not jittered or averaged, so identical coordinates overlap.\n\nThe combined panel contains both endpoints of each link: once with its consumer's degree and once with its resource's degree. It is not a set of independent links.\n\nSource: ",input,"\n\nSHA256: ",bytes2hex(sha256(read(input))),"\n\nUses the same confirmed positive-support pair records as the preceding analysis. Existing sampling and site-universe flags remain applicable.")
end

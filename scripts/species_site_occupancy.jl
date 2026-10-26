using CSV,DataFrames,CairoMakie,SHA
root=ARGS[1];dest=ARGS[2];mkpath(dest)
datasets=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
sitepath=joinpath(root,"outputs/statistical description of local support/site_universe.csv")
sites=CSV.read(sitepath,DataFrame;types=Dict(:site=>String));rows=NamedTuple[];versions=NamedTuple[]
push!(versions,(path=relpath(sitepath,root),sha256=bytes2hex(sha256(read(sitepath)))))
for ds in datasets
 path=joinpath(root,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv")
 v=CSV.read(path,DataFrame;types=String);ss=Set(sites.site[sites.dataset.==ds]);N=length(ss)
 @assert all(in.(v.site,Ref(ss))) && all(in.(v.trophic_level,Ref(["consumer","resource"])))
 for g in groupby(v,[:trophic_level,:species])
  n=length(unique(g.site));@assert 1<=n<=N
  push!(rows,(dataset=ds,role=g.trophic_level[1],species=g.species[1],occupied_sites=n,N=N,occupied_percent=100n/N))
 end
 push!(versions,(path=relpath(path,root),sha256=bytes2hex(sha256(read(path)))))
end
t=DataFrame(rows);frequencies=NamedTuple[]
for ds in datasets,role in ["consumer","resource","combined"]
 v=t[(t.dataset.==ds).&(role=="combined" ? trues(nrow(t)) : t.role.==role),:];N=only(unique(v.N))
 for n in 1:N
  countn=count(==(n),v.occupied_sites)
  push!(frequencies,(dataset=ds,group=role,occupied_sites=n,occupied_percent=100n/N,N=N,species_count=countn,total_species=nrow(v),proportion=countn/nrow(v)))
 end
end
f=DataFrame(frequencies)
set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
for (ri,role) in enumerate(["consumer","resource","combined"]),percentage in [false,true]
 fig=Figure(size=(1800,850));color=["#0072B2","#D55E00","#009E73"][ri]
 ymax=min(1.,ceil(10maximum(f.proportion[f.group.==role]))/10+.05)
 for (i,ds) in enumerate(datasets)
  v=f[(f.dataset.==ds).&(f.group.==role),:];N=only(unique(v.N));maxn=maximum(v.occupied_sites[v.species_count.>0]);v=v[v.occupied_sites.<=maxn,:]
  @assert isapprox(sum(v.proportion),1;atol=1e-8)
  ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" ")*" (N = $N)",xlabel=percentage ? "Sampled sites occupied (%)" : "Occupied sites, n",ylabel=mod1(i,5)==1 ? "Proportion of species" : "")
  x=percentage ? v.occupied_percent : v.occupied_sites
  barplot!(ax,x,v.proportion;width=percentage ? .85*100/N : .85,color=color)
  xlims!(ax,percentage ? 0 : .4,percentage ? 100+50/N : maxn+.6);ylims!(ax,0,ymax)
  if percentage;ax.xticks=0:25:100
  else;step=maxn<=10 ? 1 : maxn<=30 ? 5 : maxn<=100 ? 20 : 50;ax.xticks=sort(unique(vcat([1],collect(step:step:maxn))));end
 end
 name="0$(2ri-(percentage ? 0 : 1))_$(role)_$(percentage ? "percent" : "absolute").png"
 save(joinpath(dest,name),fig;px_per_unit=2)
end
CSV.write(joinpath(dest,"species_occupancy.csv"),t);CSV.write(joinpath(dest,"frequency_distributions.csv"),f);CSV.write(joinpath(dest,"input_versions.csv"),DataFrame(versions))
write(joinpath(dest,"README.md"),"# Species site occupancy frequencies\n\nConfirmed species-occurrence records, including recorded presences without an observed interaction. Each species' occupancy is the number of distinct sampled sites where it is recorded. N comes from the established site universe. No undocumented or wholly unrecorded species are inferred. Y is the proportion of species in each group. Absolute and percentage plots show identical discrete frequencies: percentage=100*n/N, with no rebinning. Percentage plots span the full occupancy range; absolute ranges adapt to observed occupancy. Y scales match across datasets and both versions within each group. Combined panels pool species with equal species weight, retaining trophic-role identity (not an equal-weight average of consumer/resource distributions).\n\nThese describe site occupancy, not partner specialism. Sampling effort, detectability, geographic extent and interaction-derived occurrence records affect interpretation. Confirmed restorations and unresolved occurrence flags are unchanged: Gottin PP missing-identity records, Garraf PP2 zero columns and Montseny zero row remain excluded as before. Shared sites/collections limit independence.\n\nReproduce with: julia --startup-file=no --project=. scripts/species_site_occupancy.jl REPO OUTPUT_FOLDER\n")

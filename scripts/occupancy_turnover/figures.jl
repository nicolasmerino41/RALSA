function plots()
 st=CSV.read(joinpath(DEST,"singleton_stages.csv"),DataFrame);su=CSV.read(joinpath(DEST,"dataset_summary.csv"),DataFrame);d=CSV.read(joinpath(DEST,"occupancy_distributions.csv"),DataFrame)
 set_theme!(Theme(font="Arial",fontsize=18,Axis=(xgridvisible=false,ygridvisible=false,)))
 fig=Figure(size=(1800,850));colors=["#0072B2","#D55E00","#009E73"]
 for (i,ds) in enumerate(DATA)
  ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xticks=(1:4,["Species","Weighted\nspecies","CO","Links"]),ylabel=mod1(i,5)==1 ? "Singleton frequency (%)" : "")
  for (ri,role) in enumerate(["consumer","resource","combined"])
   v=only(eachrow(st[(st.dataset.==ds).&(st.role.==role),:]));ys=100 .* [v.singleton_species,v.singleton_degree_weighted,v.singleton_C,v.singleton_K]
   lines!(ax,1:4,ys;color=colors[ri],linewidth=2);scatter!(ax,1:4,ys;color=colors[ri],markersize=8)
  end
  ylims!(ax,0,80)
 end
 Legend(fig[3,1:5],[LineElement(color=c) for c in colors],["Consumers","Resources","Combined"],orientation=:horizontal,framevisible=false)
 save(joinpath(DEST,"01_singleton_stages.png"),fig;px_per_unit=2)
 fig=Figure(size=(1150,650));ax=Axis(fig[1,1],yticks=(1:10,replace.(DATA,"_"=>" ")),yreversed=true,xlabel="Share of singleton interactions (%)")
 cols=["#0072B2","#CC79A7","#E69F00"]
 for (i,ds) in enumerate(DATA)
  r=only(eachrow(su[su.dataset.==ds,:]));v=100 .* [r.single_site_species,r.single_site_overlap,r.repeated_opportunities]./r.singletons;left=0.
  for j in 1:3;barplot!(ax,[i],[left+v[j]];fillto=left,direction=:x,color=cols[j],width=.65);left+=v[j];end
 end
 xlims!(ax,0,100)
 Legend(fig[2,1],[PolyElement(color=c) for c in cols],["At least one single-site species","Both species occupy >1 site; overlap once","Repeated co-occurrence; interaction once"],orientation=:vertical,framevisible=false)
 rowsize!(fig.layout,1,Relative(.75))
 save(joinpath(DEST,"02_singleton_origins.png"),fig;px_per_unit=2)
 fig=Figure(size=(1350,650))
 ax=Axis(fig[1,1],yticks=(1:10,replace.(DATA,"_"=>" ")),yreversed=true,xlabel="Mean pairwise dissimilarity",title="Poisot decomposition")
 bx=Axis(fig[1,2],yticks=(1:10,replace.(DATA,"_"=>" ")),yreversed=true,xlabel="Share of differing links (%)",title="Direct classification of link differences")
 for (i,ds) in enumerate(DATA)
  r=only(eachrow(su[su.dataset.==ds,:]))
  barplot!(ax,[i],[r.mean_beta_ST];direction=:x,color=cols[1],width=.65)
  barplot!(ax,[i],[r.mean_beta_WN];fillto=r.mean_beta_ST,direction=:x,color=cols[3],width=.65)
  p=100r.turnover_fraction_of_differing_links
  barplot!(bx,[i],[p];direction=:x,color=cols[1],width=.65)
  barplot!(bx,[i],[100];fillto=p,direction=:x,color=cols[3],width=.65)
 end
 xlims!(ax,0,1);xlims!(bx,0,100)
 Legend(fig[2,1:2],[PolyElement(color=cols[1]),PolyElement(color=cols[3])],["Species-turnover component","Differences among shared species"],orientation=:horizontal,framevisible=false)
 save(joinpath(DEST,"03_network_dissimilarity.png"),fig;px_per_unit=2)
 fig=Figure(size=(1800,850))
 for (i,ds) in enumerate(DATA)
  ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Occupied sites",ylabel=mod1(i,5)==1 ? "Proportion with occupancy ≥ n" : "",xscale=log10,yscale=log10)
  v=d[(d.dataset.==ds).&(d.role.=="combined"),:]
  for (j,field) in enumerate([:unweighted,:linked_only,:degree_weighted])
   y=reverse(cumsum(reverse(v[!,field])));keep=v[!,field].>0;lines!(ax,v.n[keep],y[keep];color=colors[j],linewidth=2.5)
  end
 end
 Legend(fig[3,1:5],[LineElement(color=c) for c in colors],["Each recorded species once","Each linked species once","Weighted by regional degree"],orientation=:horizontal,framevisible=false)
 save(joinpath(DEST,"04_degree_weighted_occupancy.png"),fig;px_per_unit=2)
 lines=["# Findings","","High–low–high is defined strictly as species occupancy singleton fraction > known-pair C=1 fraction, and K=1 fraction > C=1 fraction. It compares different populations, not a probability-conserving cascade."]
 for role in ["consumer","resource","combined"]
  v=st[st.role.==role,:];yes=v.dataset[v.high_low_high]
  push!(lines,"- $role: high–low–high in $(length(yes))/10 datasets ($(join(yes,", "))). Among linked species only: $(count(v.high_low_high_linked))/10. Degree weighting reduces occupancy singleton frequency in $(count(v.singleton_degree_weighted.<v.singleton_species))/10.")
 end
 for r in eachrow(su)
  push!(lines,@sprintf("- %s: singleton origins %.1f%% single-site species, %.1f%% restricted overlap, %.1f%% repeated opportunities. Mean βWN=%.3f, βOS=%.3f, βST=%.3f; %.1f%% of differing-link occurrences involve species turnover. βOS undefined for %d/%d nonempty site pairs (zero convention).",r.dataset,100r.single_site_species/r.singletons,100r.single_site_overlap/r.singletons,100r.repeated_opportunities/r.singletons,r.mean_beta_WN,r.mean_beta_OS,r.mean_beta_ST,100r.turnover_fraction_of_differing_links,r.OS_undefined_pairs,r.nonempty_pairs))
 end
 append!(lines,["","For every known link C=1 forces K=1. Species occupancy constrains the origin of those forced singletons, whereas K=1,C>1 identifies singleton interactions despite repeated recorded opportunities. These are structural categories, not causal estimates.","","Poisot 2012 uses separately normalised dissimilarities: βWN=(r+t)/(2a+r+t); βOS=r/(2a+r); βST=βWN−βOS, where a is shared interactions, r is differing interactions whose endpoints occur at both sites, and t is differing interactions with at least one endpoint absent from the other site. The subtraction is not a literal fraction of differing links. The second panel classifies the actual differences as t/(r+t), pooled across site pairs. Contrasts between panels can therefore reflect their denominators, not contradictory results.","","When no interaction among shared species occurs in either site (2a+r=0), βOS is undefined; the all-pair display uses zero by explicit convention. Counts and means restricted to defined cases are saved. Pairs with no interactions at either site are excluded from mean interaction dissimilarity. βS uses all recorded species, including species not interacting there. Each unordered site pair has equal weight in mean β; pairs are dependent and no significance tests are applied.","","Known links and confirmed occurrence records are used throughout. Every reconstructed K and C is checked against authoritative records. Degree weighting repeats a species once per regional partner; combined weighting counts both endpoints, total weight 2L. It is not identical to the fraction of pairs containing a singleton species, which is counted directly in the origins table.","","Sampling/occupancy limitations remain: many presences derive from interactions; unconfirmed Gottin PP missing-identity records, Garraf PP2 zero columns and Montseny zero row remain excluded. Shared sites and collections limit independence. Low observed occupancy does not establish a small geographic range or habitat specialism. Static interaction differences are not proof of behavioural rewiring. These analyses do not identify an abundance/compatibility mechanism or add an independent predictor of random-removal gap."])
 write(joinpath(DEST,"FINDINGS.md"),join(lines,"\n"))
end

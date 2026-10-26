function plots()
 t=CSV.read(joinpath(DEST,"fit_summary.csv"),DataFrame);d=CSV.read(joinpath(DEST,"fitted_distributions.csv"),DataFrame);emp=CSV.read(joinpath(DEST,"empirical_distributions.csv"),DataFrame)
 set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
 for (gi,group) in enumerate(["CO","SP"])
  fig=Figure(size=(1850,900))
  for (i,ds) in enumerate(DATA)
   row=only(eachrow(t[(t.dataset.==ds).&(t.group.==group).&(t.model.==MODELS[1]),:]));maxk=row.maximum_sites
   ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel=group=="CO" ? "Co-occurring sites" : "Occupied sites",ylabel=mod1(i,5)==1 ? "Proportion with ≥ n sites" : "",xscale=log10,yscale=log10)
   for (mi,model) in enumerate(MODELS)
    v=d[(d.dataset.==ds).&(d.group.==group).&(d.model.==model),:];keep=v.fitted_ccdf.>0
    lines!(ax,v.n_sites[keep],v.fitted_ccdf[keep];color=COLORS[mi],linewidth=2,linestyle=mi==5 ? :dash : :solid)
    if mi==1
     keep=v.observed_frequency.>0;lines!(ax,v.n_sites[keep],v.observed_ccdf[keep];color=:black,linewidth=2.5);scatter!(ax,v.n_sites[keep],v.observed_ccdf[keep];color=:black,markersize=4)
    end
   end
   xlims!(ax,.95,min(row.N,maxk+3)+.5);ylims!(ax,.5/row.sample_size,1.1)
   ticks=filter(<=(maxk+3),[1,2,5,10,20,50,100,200,500]);ax.xticks=(ticks,string.(ticks))
  end
  Legend(fig[3,1:5],vcat([LineElement(color=COLORS[i],linestyle=i==5 ? :dash : :solid) for i in 1:5],[LineElement(color=:black)]),vcat(LABELS,["Empirical"]),orientation=:horizontal,framevisible=false,labelsize=16)
  save(joinpath(DEST,"0$(gi)_$(group)_fitted_ccdf.png"),fig;px_per_unit=2)
  fig=Figure(size=(1550,800))
  for (j,field,title) in [(1,:delta_AIC,"ΔAIC (lower is better)"),(2,:bootstrap_p,"Bootstrap goodness-of-fit p")]
   ax=Axis(fig[1,j],title=title,xticks=(1:5,["Exponential","Power law","Cutoff","Lognormal","Singleton"]),yticks=(1:10,replace.(DATA,"_"=>" ")),xticklabelrotation=pi/6)
   for (i,ds) in enumerate(DATA),(mi,m) in enumerate(MODELS)
    r=only(eachrow(t[(t.dataset.==ds).&(t.group.==group).&(t.model.==m),:]));v=r[field];shade=j==1 ? clamp(v/30,0,1) : 1-clamp(v/.5,0,1)
    poly!(ax,Rect2f(mi-.48,i-.46,.96,.92);color=RGBAf(.94-.35shade,.97-.32shade,.98-.25shade,1))
    text!(ax,mi,i;text=j==1 ? @sprintf("%.1f",v) : @sprintf("%.3f",v),align=(:center,:center),fontsize=17)
   end
   xlims!(ax,.5,5.5);ylims!(ax,.5,10.5);ax.yreversed=true;hidespines!(ax)
  end
  save(joinpath(DEST,"0$(gi+2)_$(group)_AIC_goodness_of_fit.png"),fig;px_per_unit=2)
 end
 fig=Figure(size=(1850,900));cols=["#0072B2","#D55E00","#009E73"]
 for (i,ds) in enumerate(DATA)
  ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Number of sites",ylabel=mod1(i,5)==1 ? "Proportion with ≥ n sites" : "",xscale=log10,yscale=log10)
  for (j,group) in enumerate(["K","CO","SP"])
   v=emp[(emp.dataset.==ds).&(emp.group.==group),:];keep=v.frequency.>0;lines!(ax,v.n_sites[keep],v.ccdf[keep];color=cols[j],linewidth=3)
  end
 end
 Legend(fig[3,1:5],[LineElement(color=c) for c in cols],["K: local interaction support","CO: co-occurrence of known pairs","SP: species occupancy"],orientation=:horizontal,framevisible=false)
 save(joinpath(DEST,"05_K_CO_SP_empirical_ccdf.png"),fig;px_per_unit=2)
 lines=["# Distribution-family results","","Fits describe the complete positive distribution on integers 1,...,N. CO uses known interacting pairs, matching the earlier K/CO plots. SP pools consumers and resources with equal species weight, including confirmed species never observed interacting. K is included empirically only, not refitted here.","","| Dataset | Quantity | Lowest AIC | AIC | Bootstrap p | Other models within 2 AIC |","|---|---|---|---:|---:|---|"]
 winners=NamedTuple[]
 for ds in DATA,group in ["CO","SP"]
  v=t[(t.dataset.==ds).&(t.group.==group),:];r=v[argmin(v.AIC),:];others=v.model[(v.delta_AIC.<=2).&(v.model.!=r.model)]
  push!(winners,(dataset=ds,group=group,best_AIC_model=r.model,AIC=r.AIC,bootstrap_p=r.bootstrap_p,best_AICc_model=v.model[argmin(v.AICc)],any_model_p_ge_005=any(v.bootstrap_p.>=.05)))
  push!(lines,"| $ds | $group | $(r.model) | $(round(r.AIC;digits=1)) | $(round(r.bootstrap_p;digits=3)) | $(join(others,", ")) |")
 end
 CSV.write(joinpath(DEST,"best_models.csv"),DataFrame(winners))
 append!(lines,["","AIC ranks alternatives, while refitted parametric-bootstrap KS goodness-of-fit tests assess detectable mismatch. A best-ranked model may still fit poorly. p≥0.05 does not establish the true family; compare curves, AIC differences, boundaries and Monte Carlo errors. AICc is also provided because species samples can be small. AIC values cannot be compared between SP and CO or between datasets with different samples.","","499 simulations per fit, each refitted. Candidate definitions match the previous K analysis: decreasing bounded exponential/geometric, decreasing bounded power law, power law with exponential cutoff, discrete lognormal density kernel, and singleton-plus-geometric. All are normalised over 1,...,N. This lognormal evaluates a density kernel at integers rather than integrating rounding bins. A zero-curvature lognormal boundary represents its limiting power-law kernel, not a finite lognormal. There is no tail cutoff selection or unbounded scale-free inference here.","","Observations share species/sites and sites have unequal sampling; independent-observation likelihoods and bootstrap tests are working diagnostics, not dependence-robust ecological tests. Shared campaigns further limit independence. Positive K/CO and all-species SP refer to distinct observational units. Similar winning families do not imply equal parameters, common mechanisms or identical shapes.","","All source hashes and reconstructed C were verified. Confirmed occurrence restorations and existing data flags remain unchanged: uncertain Gottin PP missing-identity records, Garraf PP2 zero columns and Montseny zero row are excluded as before; many presences are interaction-derived. Full N includes the established site universe but not undocumented sites. Existing analyses are preserved."])
 push!(lines,"\nBootstrap convergence flags: $(count(t.bootstrap_convergence_flag)) of $(nrow(t)) fits have at least one replicate KKT residual ≥1e-6; consult fit_summary.csv before interpreting borderline p-values.")
 write(joinpath(DEST,"FINDINGS.md"),join(lines,"\n"))
 CSV.write(joinpath(DEST,"code_versions.csv"),DataFrame(file=["run_analysis.jl","figures.jl","distribution_helpers.jl"],sha256=[bytes2hex(sha256(read(joinpath(@__DIR__,p)))) for p in ["run_analysis.jl","figures.jl","distribution_helpers.jl"]]))
end

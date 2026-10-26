function addlegend(f)
 Legend(f[3,1:5],[LineElement(color=COLORS[m],linestyle=m=="Empirical" ? :dash : :solid) for m in MODELS],["Empirical","Common, fixed regional M","Persistent pair, fixed regional M"];orientation=:horizontal,framevisible=false)
end
function plots()
 t=CSV.read(joinpath(DEST,"accumulation_summary.csv"),DataFrame);ex=CSV.read(joinpath(DEST,"exact_expectations.csv"),DataFrame);inp=CSV.read(joinpath(DEST,"dataset_inputs.csv"),DataFrame)
 set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
 for (field,label,file) in [("links","Regional links accumulated (%)","01_regional_link_accumulation"),("mean_support","Mean support per observed link","02_mean_support_accumulation"),("gini","Support Gini (%)","03_concentration_accumulation"),("singleton","Singleton links (%)","05_singleton_accumulation")]
  f=Figure(size=(1800,880))
  for (i,ds) in enumerate(DATA)
   ax=Axis(f[cld(i,5),mod1(i,5)];title=replace(ds,"_"=>" "),xlabel="Sites accumulated (%)",ylabel=mod1(i,5)==1 ? label : "")
   L=only(inp.L[inp.dataset.==ds]);scale=field=="links" ? 100/L : field=="mean_support" ? 1. : 100.
   for model in MODELS
    v=t[(t.dataset.==ds).&(t.model.==model).&(t.metric.==field),:];sort!(v,:n)
    band!(ax,100v.site_fraction,scale*v.lower,scale*v.upper;color=(COLORS[model],.10))
    if field=="links"
     e=ex[(ex.dataset.==ds).&(ex.model.==model),:];lines!(ax,100e.site_fraction,100e.expected_link_fraction;color=COLORS[model],linewidth=2.5,linestyle=model=="Empirical" ? :dash : :solid)
    else
     lines!(ax,100v.site_fraction,scale*v.mean;color=COLORS[model],linewidth=2.5,linestyle=model=="Empirical" ? :dash : :solid)
    end
   end
   xlims!(ax,0,100);field in ["links","singleton"] ? ylims!(ax,0,100) : field=="gini" ? ylims!(ax,0,65) : ylims!(ax,0,nothing)
  end
  addlegend(f);save(joinpath(DEST,file*".png"),f;px_per_unit=2)
 end
 f=Figure(size=(1800,880))
 for (i,ds) in enumerate(DATA)
  ax=Axis(f[cld(i,5),mod1(i,5)];title=replace(ds,"_"=>" "),xlabel="Regional links accumulated (%)",ylabel=mod1(i,5)==1 ? "Support Gini (%)" : "")
  for model in MODELS
   v=t[(t.dataset.==ds).&(t.model.==model).&(t.metric.=="gini"),:];e=ex[(ex.dataset.==ds).&(ex.model.==model).&(ex.n.>0),:];@assert v.n==e.n
   lines!(ax,100e.expected_link_fraction,100v.mean;color=COLORS[model],linewidth=2.5,linestyle=model=="Empirical" ? :dash : :solid)
   markers=[argmin(abs.(v.site_fraction.-q)) for q in [.25,.5,1.]]
   scatter!(ax,100e.expected_link_fraction[markers],100v.mean[markers];color=COLORS[model],markersize=[6,9,12])
  end
  xlims!(ax,0,100);ylims!(ax,0,65)
 end
 addlegend(f);save(joinpath(DEST,"04_diversity_and_concentration.png"),f;px_per_unit=2)
end
function report()
    t=CSV.read(joinpath(DEST,"accumulation_summary.csv"),DataFrame);ex=CSV.read(joinpath(DEST,"exact_expectations.csv"),DataFrame);inp=CSV.read(joinpath(DEST,"dataset_inputs.csv"),DataFrame)
 draws=CSV.read(joinpath(DEST,"accumulation_draws.csv"),DataFrame);validation=NamedTuple[]
 for v in groupby(draws,[:dataset,:model])
  ds=v.dataset[1];model=v.model[1];rep=combine(groupby(v,:replicate),:links=>mean=>:links,:support=>mean=>:support)
  truth=ex[(ex.dataset.==ds).&(ex.model.==model).&(ex.n.>0),:]
  @assert all((v.links.==0).|(abs.(v.mean_support.*v.links.-v.support).<1e-8))
  for (field,target) in [(:links,:expected_links),(:support,:expected_total_support)]
   error=abs(mean(rep[!,field])-mean(truth[!,target]));se=std(rep[!,field])/sqrt(nrow(rep));@assert error<=6se+1e-8
   push!(validation,(dataset=ds,model=model,quantity=String(field),integrated_mean_error=error,mcse=se,passed=true))
  end
 end
 CSV.write(joinpath(DEST,"analytical_checks.csv"),DataFrame(validation))
 stages=NamedTuple[]
 for ds in DATA,q in [.25,.5,1.]
  N=only(inp.N[inp.dataset.==ds]);n=clamp(round(Int,q*N),1,N)
  for model in MODELS
   getmean(field)=only(t.mean[(t.dataset.==ds).&(t.model.==model).&(t.n.==n).&(t.metric.==field)])
   er=only(eachrow(ex[(ex.dataset.==ds).&(ex.model.==model).&(ex.n.==n),:]))
   push!(stages,(dataset=ds,model=model,target_fraction=q,n=n,N=N,actual_fraction=n/N,expected_link_percent=100er.expected_link_fraction,mean_support=getmean("mean_support"),gini_percent=100getmean("gini"),singleton_percent=100getmean("singleton"),top10_support_percent=100getmean("top10")))
  end
 end
 st=DataFrame(stages);CSV.write(joinpath(DEST,"stage_summary.csv"),st)
 fmt(x)=@sprintf("%.2f",x)
 open(joinpath(DEST,"FINDINGS.md"),"w") do io
  println(io,"# Regional diversity and support accumulate together\n\nThree hundred paired random site orders per dataset compare observed interactions with Common and persistent-Pair allocations. Each model allocation preserves the observed full-region M, all known links and their opportunity sites. Model parameters are unchanged nine-dataset estimates. No parameters are fitted to accumulation curves.\n")
  println(io,"The recurring difference is visible during accumulation, not only in the completed regional network: empirical support is more concentrated than Common at the tested quarter, half and full site counts in all ten datasets. Pair is closer for Gini in all ten at each of those checkpoints. At the half-site checkpoint Common also discovers a larger fraction of regional links than observed in every dataset. With the same expected total occurrence budget, its support is spread across more distinct relationships; empirically, more support is concentrated in a subset while other relationships remain restricted.\n\nThis supports a scale-dependent departure from homogeneous realisation under the supplied opportunities. It does not show a universal family or prove a specific heterogeneity mechanism. Intermediate budgets are equal in expectation, not fixed separately in each subsample, and model differences include their implied placement among sites.\n")
  println(io,"| Dataset | Sites at half-sample | Links accumulated: empirical % | Common % | Pair % | Gini: empirical % | Common % | Pair % |\n|---|---:|---:|---:|---:|---:|---:|---:|")
  for ds in DATA
   v=st[(st.dataset.==ds).&(st.target_fraction.==.5),:];r=[only(eachrow(v[v.model.==m,:])) for m in MODELS]
   println(io,"| $(replace(ds,"_"=>" ")) | $(r[1].n)/$(r[1].N) | $(join(fmt.([x.expected_link_percent for x in r])," | ")) | $(join(fmt.([x.gini_percent for x in r])," | ")) |")
  end
  println(io,"\n## What follows automatically, and what is informative\n\nAt every realised step, mean support equals accumulated local support divided by accumulated regional links. Every link has support one when it first appears. At one nonempty site Gini is zero; full-region mean support is fixed to M/L in all models. Those endpoints and identities are not new ecological findings. Full-region concentration is not fixed.\n\nThe comparison asks how quickly links are discovered and how unequally support builds up between these endpoints, given identical opportunities and full-region budgets. Expected total accumulated support is nM/N for every model, but realised intermediate totals need not match: site totals are not constrained. Thus intermediate nonlinear concentration differences are not a pure fixed-intermediate-M allocation test. The paired difference table retains support totals alongside other metrics.\n")
  for q in [.25,.5,1.]
   v=st[st.target_fraction.==q,:];o=v[v.model.=="Empirical",:];c=v[v.model.=="Common",:];p=v[v.model.=="Pair",:]
   println(io,"At nominal $(round(Int,100q))% sites, empirical mean Gini exceeds Common in $(count(o.gini_percent.>c.gini_percent))/10 datasets. Pair is closer to empirical Gini than Common in $(count(abs.(p.gini_percent.-o.gini_percent).<abs.(c.gini_percent.-o.gini_percent)))/10. These are descriptive comparisons of means, not significance counts.")
  end
  println(io,"\n## How to read the figures\n\nFigures 1–3 follow regional link richness, mean positive support and Gini against the fraction of sites accumulated. Figure 4 plots the mean concentration trajectory against exact expected link richness; dots indicate nearest available 25%, 50% and 100% site counts, with increasing size. Figure 5 tracks singleton abundance. Link-richness means are analytical expectations; concentration and mean-support curves average replicate-level statistics, not ratios of averages.\n\nEmpirical bands are central 95% variation over random site orders. Model bands include both random allocation and site order. They are not confidence intervals and do not propagate parameter uncertainty. Empty prefixes have link richness and support zero but undefined mean support/concentration; these are excluded from those conditional averages and counted explicitly. Full curves contain every integer site count.\n\n## Scope\n\nRandom ordering is accumulation of sampled communities, not chronological assembly or geographic expansion. Expected link accumulation under uniform sampling is mathematically determined by the full K distribution; it is the complementary view of the existing random-removal calculation, not independent corroboration of that distribution. Site-level arrangements also affect joint fluctuations and nonlinear statistics. This exercise visualises and tests scale-dependent implications of the existing allocation models, rather than discovering an independent mechanism.\n\nKnown identities, confirmed occurrence inputs, observation limitations and shared-site/collection dependence remain. Opportunity-only corrections can alter permitted placement without changing original K. No spatial coordinates, new turnover model or directed habitat-loss scenario are introduced. See README and provenance files.")
 end
 CSV.write(joinpath(DEST,"code_versions.csv"),DataFrame(file=["run_analysis.jl","figures.jl"],sha256=[bytes2hex(sha256(read(joinpath(@__DIR__,f)))) for f in ["run_analysis.jl","figures.jl"]]))
end

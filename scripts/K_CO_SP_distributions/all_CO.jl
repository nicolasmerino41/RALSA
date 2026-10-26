module AllCoDistributions
using CSV,DataFrames,Statistics,Random,LinearAlgebra,SHA,CairoMakie,Printf
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"../..")) : abspath(ARGS[1])
const BASE=joinpath(ROOT,"outputs/statistical description of local support/K_CO_SP_distributions")
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : BASE
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const MODELS=["Exponential","Power_law","Cutoff_power_law","Discrete_lognormal","Singleton_geometric"]
const LABELS=["Exponential","Power law","Cutoff power law","Discrete lognormal","Singleton + geometric"]
const COLORS=["#888888","#CC8844","#4477AA","#228877","#882255"]
const B=499
include(joinpath(ROOT,"scripts/K_CO_SP_distributions/distribution_helpers.jl"))
mkpath(DEST)
function main()
 versions=CSV.read(joinpath(BASE,"input_versions.csv"),DataFrame)
 for r in eachrow(versions);@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256;end
 source=joinpath(ROOT,"outputs/support_prediction_zero/confirmed/pairs.csv");pairs=CSV.read(source,DataFrame)
 qc=CSV.read(joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory/quality_control.csv"),DataFrame)
 old=CSV.read(joinpath(BASE,"empirical_distributions.csv"),DataFrame)
 rows=NamedTuple[];curves=NamedTuple[];emp=NamedTuple[];checks=NamedTuple[]
 for (di,ds) in enumerate(DATA)
  v=pairs[(pairs.dataset.==ds).&(pairs.n.>=1),:];N=only(unique(v.N));x=Int.(v.n);n=length(x)
  @assert all(0 .<=v.K .<=v.n .<=N)
  known=v[v.K.>0,:];q=qc[qc.dataset.==ds,:]
  @assert Dict((r.consumer,r.resource)=>(r.K,r.n) for r in eachrow(known))==Dict((r.consumer,r.resource)=>(r.K,r.C) for r in eachrow(q))
  occ=CSV.read(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv"),DataFrame;types=String)
  pres=Dict((g.trophic_level[1],g.species[1])=>Set(g.site) for g in groupby(occ,[:trophic_level,:species]))
  @assert all(length(intersect(pres[("consumer",r.consumer)],pres[("resource",r.resource)]))==r.n for r in eachrow(v))
  for group in ["K","CO"]
   prior=old[(old.dataset.==ds).&(old.group.==group),:];sort!(prior,:n_sites)
   vals=group=="K" ? known.K : known.n
   @assert prior.N[1]==N && prior.count==[count(==(k),vals) for k in 1:N]
  end
  counts=[count(==(k),x) for k in 1:N];obs=counts/n;cc=reverse(cumsum(reverse(obs)))
  for k in 1:N;push!(emp,(dataset=ds,group="All_CO",N=N,sample_size=n,n_sites=k,count=counts[k],frequency=obs[k],ccdf=cc[k]));end
  push!(checks,(dataset=ds,N=N,all_cooccurring_pairs=n,known_pairs=nrow(known),never_interacting_pairs=count(v.K.==0),known_K_CO_match=true,occurrence_C_match=true))
  for (mi,m) in enumerate(MODELS)
   f=fit(counts,m);@assert f.kkt<1e-6
   boot=bootstrap(counts,m,f,MersenneTwister(910000+100di+mi))
   push!(rows,(dataset=ds,group="All_CO",model=m,N=N,sample_size=n,maximum_sites=maximum(x),parameters=f.parameters,log_likelihood=f.ll,AIC=f.aic,AICc=f.aicc,KS=f.D,bootstrap_p=boot.p,bootstrap_mcse=boot.mcse,bootstrap_replicates=B,alpha=f.alpha,lambda=f.lambda,mu=f.mu,sigma=f.sigma,beta=f.beta,s=f.s,r=f.r,boundary=f.boundary,kkt_residual=f.kkt,bootstrap_max_kkt=boot.maxkkt,bootstrap_convergence_flag=boot.maxkkt>=1e-6))
   fcc=reverse(cumsum(reverse(f.p)))
   for k in 1:N;push!(curves,(dataset=ds,model=m,n_sites=k,observed_frequency=obs[k],observed_ccdf=cc[k],fitted_frequency=f.p[k],fitted_ccdf=fcc[k]));end
   println(ds," ",m," AIC=",round(f.aic;digits=1)," p=",boot.p);flush(stdout)
  end
  CSV.write(joinpath(DEST,"all_CO_fit_summary.csv"),DataFrame(rows));CSV.write(joinpath(DEST,"all_CO_fitted_distributions.csv"),DataFrame(curves));CSV.write(joinpath(DEST,"all_CO_empirical_distributions.csv"),DataFrame(emp))
 end
 t=DataFrame(rows);t.delta_AIC=[r.AIC-minimum(t.AIC[t.dataset.==r.dataset]) for r in eachrow(t)];t.delta_AICc=[r.AICc-minimum(t.AICc[t.dataset.==r.dataset]) for r in eachrow(t)]
 CSV.write(joinpath(DEST,"all_CO_fit_summary.csv"),t);CSV.write(joinpath(DEST,"all_CO_checks.csv"),DataFrame(checks))
 paths=[source,joinpath(BASE,"empirical_distributions.csv"),joinpath(ROOT,"scripts/K_CO_SP_distributions/distribution_helpers.jl")]
 CSV.write(joinpath(DEST,"all_CO_source_versions.csv"),DataFrame(path=relpath.(paths,ROOT),sha256=[bytes2hex(sha256(read(p))) for p in paths]))
end
function plots()
 t=CSV.read(joinpath(DEST,"all_CO_fit_summary.csv"),DataFrame);d=CSV.read(joinpath(DEST,"all_CO_fitted_distributions.csv"),DataFrame)
 set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
 fig=Figure(size=(1850,900))
 for (i,ds) in enumerate(DATA)
  row=only(eachrow(t[(t.dataset.==ds).&(t.model.==MODELS[1]),:]))
  ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Co-occurring sites",ylabel=mod1(i,5)==1 ? "Proportion of all pairs with ≥ n sites" : "",xscale=log10,yscale=log10)
  for (mi,m) in enumerate(MODELS)
   v=d[(d.dataset.==ds).&(d.model.==m),:];keep=v.fitted_ccdf.>0
   lines!(ax,v.n_sites[keep],v.fitted_ccdf[keep];color=COLORS[mi],linewidth=2,linestyle=mi==5 ? :dash : :solid)
   if mi==1;keep=v.observed_frequency.>0;lines!(ax,v.n_sites[keep],v.observed_ccdf[keep];color=:black,linewidth=2.5);end
  end
  xlims!(ax,.95,min(row.N,row.maximum_sites+3)+.5);ylims!(ax,.5/row.sample_size,1.1)
 end
 Legend(fig[3,1:5],vcat([LineElement(color=COLORS[i],linestyle=i==5 ? :dash : :solid) for i in 1:5],[LineElement(color=:black)]),vcat(LABELS,["Empirical"]),orientation=:horizontal,framevisible=false,labelsize=16)
 save(joinpath(DEST,"06_All_CO_fitted_ccdf.png"),fig;px_per_unit=2)
 fig=Figure(size=(1550,800))
 for (j,field,title) in [(1,:delta_AIC,"ΔAIC (lower is better)"),(2,:bootstrap_p,"Bootstrap goodness-of-fit p")]
  ax=Axis(fig[1,j],title=title,xticks=(1:5,["Exponential","Power law","Cutoff","Lognormal","Singleton"]),yticks=(1:10,replace.(DATA,"_"=>" ")),xticklabelrotation=pi/6)
  for (i,ds) in enumerate(DATA),(mi,m) in enumerate(MODELS)
   r=only(eachrow(t[(t.dataset.==ds).&(t.model.==m),:]));v=r[field];shade=j==1 ? clamp(v/30,0,1) : 1-clamp(v/.5,0,1)
   poly!(ax,Rect2f(mi-.48,i-.46,.96,.92);color=RGBAf(.94-.35shade,.97-.32shade,.98-.25shade,1))
   text!(ax,mi,i;text=j==1 ? @sprintf("%.1f",v) : @sprintf("%.3f",v),align=(:center,:center),fontsize=17)
  end
  xlims!(ax,.5,5.5);ylims!(ax,.5,10.5);ax.yreversed=true;hidespines!(ax)
 end
 save(joinpath(DEST,"07_All_CO_AIC_goodness_of_fit.png"),fig;px_per_unit=2)
 old=CSV.read(joinpath(BASE,"empirical_distributions.csv"),DataFrame);allco=CSV.read(joinpath(DEST,"all_CO_empirical_distributions.csv"),DataFrame);emp=vcat(old,allco)
 CSV.write(joinpath(DEST,"four_curve_empirical_distributions.csv"),emp)
 fig=Figure(size=(1850,900));cols=["#0072B2","#D55E00","#009E73","#CC79A7"]
 for (i,ds) in enumerate(DATA)
  ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Number of sites",ylabel=mod1(i,5)==1 ? "Proportion with ≥ n sites" : "",xscale=log10,yscale=log10)
  for (j,group) in enumerate(["K","CO","SP","All_CO"])
   v=emp[(emp.dataset.==ds).&(emp.group.==group),:];keep=v.frequency.>0;lines!(ax,v.n_sites[keep],v.ccdf[keep];color=cols[j],linewidth=3,linestyle=j==4 ? :dash : :solid)
  end
 end
 Legend(fig[3,1:5],[LineElement(color=c,linestyle=j==4 ? :dash : :solid) for (j,c) in enumerate(cols)],["K: interaction support","CO: known interacting pairs","SP: species occupancy","All CO: every co-occurring pair"],orientation=:horizontal,framevisible=false,labelsize=17)
 save(joinpath(DEST,"08_K_CO_SP_All_CO_empirical_ccdf.png"),fig;px_per_unit=2)
 lines=["# All-pair co-occurrence distribution fits","","Every pair with C≥1, including K=0, from confirmed/pairs.csv. Same five finite-support families, maximum-likelihood methods and 499 refitted parametric-bootstrap KS replicates as the known-pair analysis. Original figures and tables are unchanged; 08 adds All CO to the empirical comparison.","","| Dataset | Best AIC model | AIC | Bootstrap p | Any candidate p≥0.05 |","|---|---|---:|---:|---|"]
 winners=NamedTuple[]
 for ds in DATA
  v=t[t.dataset.==ds,:];r=v[argmin(v.AIC),:];ok=any(v.bootstrap_p.>=.05)
  push!(winners,(dataset=ds,model=r.model,AIC=r.AIC,AICc=r.AICc,bootstrap_p=r.bootstrap_p,any_model_p_ge_005=ok))
  push!(lines,"| $ds | $(r.model) | $(round(r.AIC;digits=1)) | $(r.bootstrap_p) | $ok |")
 end
 append!(lines,["","AIC only ranks models for this particular sample. Do not compare absolute AIC with known-pair AIC, because sample sizes and populations differ. Goodness-of-fit p<0.05 flags mismatch under the independent-pair working model; p≥0.05 does not prove a true family. Dependence through shared species/sites, larger all-pair sample sizes, and sampling affect these diagnostics. Near-threshold p-values have Monte Carlo uncertainty. All current records and original known K/CO counts were verified against confirmed occurrences. Existing occurrence corrections and exclusions are unchanged; see the original README.","","K and known-pair CO refer only to observed regional links; All CO includes co-occurring pairs with no observed interaction; SP pools consumers/resources once per species. Each curve has its own denominator. Positive-C pairs only: C=0 combinations are not included. No original fits were rerun.","","Bootstrap convergence flags: $(count(t.bootstrap_convergence_flag)) of $(nrow(t)) fits."])
 write(joinpath(DEST,"ALL_CO_FINDINGS.md"),join(lines,"\n"));CSV.write(joinpath(DEST,"all_CO_best_models.csv"),DataFrame(winners))
end
if !(length(ARGS)>2 && ARGS[3]=="plot");main();end
plots()
end

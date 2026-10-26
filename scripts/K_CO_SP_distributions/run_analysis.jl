module KCoSpDistributions
using CSV,DataFrames,Statistics,Random,LinearAlgebra,SHA,CairoMakie,Printf
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"../..")) : abspath(ARGS[1])
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(ROOT,"outputs/statistical description of local support/K_CO_SP_distributions")
const SOURCE=joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const MODELS=["Exponential","Power_law","Cutoff_power_law","Discrete_lognormal","Singleton_geometric"]
const LABELS=["Exponential","Power law","Cutoff power law","Discrete lognormal","Singleton + geometric"]
const COLORS=["#888888","#CC8844","#4477AA","#228877","#882255"]
const B=499
include("distribution_helpers.jl")
mkpath(DEST)
function main()
 versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
 for r in eachrow(versions);@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256;end
 used=NamedTuple[]
 function load(path;kwargs...)
  push!(used,(path=relpath(path,ROOT),sha256=bytes2hex(sha256(read(path)))));CSV.read(path,DataFrame;kwargs...)
 end
 qc=load(joinpath(SOURCE,"quality_control.csv"));sites=load(joinpath(ROOT,"outputs/statistical description of local support/site_universe.csv");types=Dict(:site=>String))
 rows=NamedTuple[];curves=NamedTuple[];emp=NamedTuple[];observations=NamedTuple[]
 for (di,ds) in enumerate(DATA)
  v=qc[qc.dataset.==ds,:];N=only(unique(v.N));ss=Set(sites.site[sites.dataset.==ds]);@assert length(ss)==N
  occ=load(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv");types=String)
  pres=Dict((g.trophic_level[1],g.species[1])=>Set(g.site) for g in groupby(occ,[:trophic_level,:species]))
  @assert all(s in ss for locs in Base.values(pres) for s in locs)
  @assert [length(intersect(pres[("consumer",r.consumer)],pres[("resource",r.resource)])) for r in eachrow(v)]==v.C
  for (gi,group) in enumerate(["K","CO","SP"])
   values=group=="K" ? Int.(v.K) : group=="CO" ? Int.(v.C) : [length(pres[key]) for key in sort(collect(keys(pres)))]
   @assert all(1 .<=values .<=N)
   counts=[count(==(k),values) for k in 1:N];n=length(values);p=counts/n;cc=reverse(cumsum(reverse(p)))
   for k in 1:N;push!(emp,(dataset=ds,group=group,N=N,sample_size=n,n_sites=k,count=counts[k],frequency=p[k],ccdf=cc[k]));end
   for (i,x) in enumerate(values);push!(observations,(dataset=ds,group=group,index=i,n_sites=x,N=N));end
   group=="K" && continue
   for (mi,model) in enumerate(MODELS)
    f=fit(counts,model);@assert f.kkt<1e-6 "Fit convergence $ds $group $model $(f.kkt)"
    boot=bootstrap(counts,model,f,MersenneTwister(884000+1000di+100gi+mi))
    push!(rows,(dataset=ds,group=group,model=model,N=N,sample_size=n,maximum_sites=maximum(values),parameters=f.parameters,log_likelihood=f.ll,AIC=f.aic,AICc=f.aicc,KS=f.D,bootstrap_p=boot.p,bootstrap_mcse=boot.mcse,bootstrap_replicates=B,alpha=f.alpha,lambda=f.lambda,mu=f.mu,sigma=f.sigma,beta=f.beta,s=f.s,r=f.r,boundary=f.boundary,kkt_residual=f.kkt,bootstrap_max_kkt=boot.maxkkt,bootstrap_convergence_flag=boot.maxkkt>=1e-6))
    fcc=reverse(cumsum(reverse(f.p)))
    for k in 1:N;push!(curves,(dataset=ds,group=group,model=model,n_sites=k,observed_frequency=p[k],observed_ccdf=cc[k],fitted_frequency=f.p[k],fitted_ccdf=fcc[k]));end
    println(ds," ",group," ",model," AIC=",round(f.aic;digits=1)," p=",boot.p);flush(stdout)
   end
  end
  CSV.write(joinpath(DEST,"fit_summary.csv"),DataFrame(rows));CSV.write(joinpath(DEST,"fitted_distributions.csv"),DataFrame(curves));CSV.write(joinpath(DEST,"empirical_distributions.csv"),DataFrame(emp));CSV.write(joinpath(DEST,"observations.csv"),DataFrame(observations))
 end
 t=DataFrame(rows)
 t.delta_AIC=[r.AIC-minimum(t.AIC[(t.dataset.==r.dataset).&(t.group.==r.group)]) for r in eachrow(t)]
 t.delta_AICc=[r.AICc-minimum(t.AICc[(t.dataset.==r.dataset).&(t.group.==r.group)]) for r in eachrow(t)]
 CSV.write(joinpath(DEST,"fit_summary.csv"),t);CSV.write(joinpath(DEST,"input_versions.csv"),versions)
 original=joinpath(ROOT,"scripts/support_distribution_families/fit_distributions.jl");push!(used,(path=relpath(original,ROOT),sha256=bytes2hex(sha256(read(original)))))
 CSV.write(joinpath(DEST,"source_versions.csv"),DataFrame(used))
 for file in ["metadata_relationships.csv","metadata_source_versions.csv"];cp(joinpath(SOURCE,file),joinpath(DEST,file);force=true);end
end
include("figures.jl")
if !(length(ARGS)>2 && ARGS[3]=="plot");main();end
plots()
end

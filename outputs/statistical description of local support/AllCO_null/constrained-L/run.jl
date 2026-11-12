using CSV, DataFrames, Statistics, Random, SHA, CairoMakie, Printf
const ROOT=normpath(joinpath(@__DIR__,"../../../.."))
const OUT=@__DIR__
const SOURCE=joinpath(ROOT,"outputs/support_prediction_zero/confirmed")
const B=5000
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const BLUE="#0072B2"
const RED="#D55E00"
function gapweights(N)
 [mean(m/N-(m<k ? 0. : prod((m-j)/(N-j) for j in 0:k-1)) for m in 0:fld(4N,5)) for k in 1:N]
end
include(joinpath(OUT,"sampler.jl"))
function main()
 sampler_test();diagnosticrows=NamedTuple[]
 pairs=CSV.read(joinpath(SOURCE,"pairs.csv"),DataFrame)

 versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
 paths=[joinpath(SOURCE,"pairs.csv"),joinpath(SOURCE,"input_versions.csv"),@__FILE__,joinpath(OUT,"sampler.jl")]
 summaries=NamedTuple[]; draws=NamedTuple[]; histrows=NamedTuple[]; checks=NamedTuple[]
 for (di,ds) in enumerate(DATA)
  g=pairs[(pairs.dataset.==ds).&(pairs.n.>0),:]
  N=only(unique(g.N)); C=Int.(g.n); K=Int.(g.K); P=length(C); T=sum(C); M=sum(K); L=count(>(0),K)
  @assert all(0 .<=K .<=C .<=N) && M>0
  @assert nrow(unique(g,[:consumer,:resource]))==P
  occpath=joinpath(SOURCE,"inputs",ds*"_occupancy.csv");push!(paths,occpath)
  vr=only(eachrow(versions[(versions.dataset.==ds).&(versions.kind.=="occupancy"),:]))
  @assert bytes2hex(sha256(read(occpath)))==vr.sha256
  occ=CSV.read(occpath,DataFrame;types=String)
  pres=Dict((x.trophic_level[1],x.species[1])=>Set(x.site) for x in groupby(occ,[:trophic_level,:species]))
  expected=Dict((a[2],b[2])=>length(intersect(sa,sb)) for (a,sa) in pres for (b,sb) in pres if a[1]=="consumer" && b[1]=="resource" && !isempty(intersect(sa,sb)))
  @assert expected==Dict((r.consumer,r.resource)=>r.n for r in eachrow(g))
  known=g[g.K.>0,:]
  ip=joinpath(SOURCE,"inputs",ds*"_interactions.csv");push!(paths,ip)
  iv=only(eachrow(versions[(versions.dataset.==ds).&(versions.kind.=="interactions"),:]))
  @assert bytes2hex(sha256(read(ip)))==iv.sha256
  ints=CSV.read(ip,DataFrame;types=String)
  observed=Dict((x.consumer[1],x.resource[1])=>length(unique(x.site)) for x in groupby(ints,[:consumer,:resource]))
  @assert Dict((r.consumer,r.resource)=>r.K for r in eachrow(known))==observed
  owners=reduce(vcat,[fill(i,C[i]) for i in 1:P]);weights=gapweights(N)
  obsfreq=[count(==(k),K)/L for k in 1:N]
  emp=[Float64(L),100(P-L)/P,100count(==(1),K)/L,100sum(weights[k] for k in K if k>0)/L]
  sims=zeros(B,4); freq=zeros(B,N);ccdf=zeros(B,N);rng=MersenneTwister(202611120+di)
  lc=logcombinations(C); perchain=B÷4; monitor=zeros(perchain,4,4)
  for chain in 1:4
   sim=copy(K)
   if chain>1
    chosen=sortperm(C.*(.5 .+rand(rng,P));rev=true)[1:L]
    if sum(C[chosen])>=M
     fill!(sim,0);sim[chosen].=1;remaining=M-L
     for i in shuffle(rng,chosen);take=min(C[i]-1,remaining);sim[i]+=take;remaining-=take;end
    end
   end
   active=findall(>(0),sim)
   for _ in 1:500;sweep!(rng,sim,C,active,lc);end
   for draw in 1:perchain
    for _ in 1:5;sweep!(rng,sim,C,active,lc);end
    b=(chain-1)*perchain+draw
    @assert sum(sim)==M && count(>(0),sim)==L && all(0 .<=sim .<=C)
    for k in sim;k>0 && (freq[b,k]+=1/L);end
    ccdf[b,:]=reverse(cumsum(reverse(freq[b,:])))
    sims[b,:]=[L,100(P-L)/P,100freq[b,1],100sum(freq[b,:].*weights)]
    overlap=count((sim.>0).&(K.>0))/L
    monitor[draw,chain,:]=[sims[b,3],sims[b,4],mean(C[sim.>0]),overlap]
    push!(draws,(dataset=ds,chain=chain,replicate=draw,L=L,zero_pct=sims[b,2],singleton_pct=sims[b,3],gap_pp=sims[b,4],known_link_overlap=overlap))
   end
  end
  for (j,metric) in enumerate(["singleton_pct","gap_pp","selected_mean_C","known_link_overlap"])
   z=diagnostics(monitor[:,:,j]);push!(diagnosticrows,(dataset=ds,metric=metric,split_Rhat=z.rhat,ESS=z.ess,mean_mcse=z.mcse,passes=z.rhat<1.01 && z.ess>=1000))
  end
  push!(checks,(dataset=ds,N=N,pairs=P,total_opportunities=T,M=M,L=L,all_pairs_verified=true,known_support_verified=true,all_draws_M_L_fixed=true))
  for (j,metric) in enumerate(["regional_links","zero_pct","singleton_pct","gap_pp"])
   v=sims[:,j];lo,hi=quantile(v,[.025,.975])
   push!(summaries,(dataset=ds,metric=metric,empirical=emp[j],null_mean=mean(v),null_lower=lo,null_upper=hi,difference=emp[j]-mean(v),position=emp[j]<lo ? "below" : emp[j]>hi ? "above" : "inside",mean_mcse=j<=2 ? 0. : diagnostics(monitor[:,:,j-2]).mcse,replicates=B))
  end
  obscc=reverse(cumsum(reverse(obsfreq)))
  for k in 1:N
   push!(histrows,(dataset=ds,K=k,empirical_frequency=obsfreq[k],null_mean_frequency=mean(freq[:,k]),null_lower_frequency=quantile(freq[:,k],.025),null_upper_frequency=quantile(freq[:,k],.975),empirical_ccdf=obscc[k],null_mean_ccdf=mean(ccdf[:,k]),null_lower_ccdf=quantile(ccdf[:,k],.025),null_upper_ccdf=quantile(ccdf[:,k],.975)))
  end
  println(ds," L observed/null: ",L," / ",round(mean(sims[:,1]);digits=1),"; gap: ",round(emp[4];digits=2)," / ",round(mean(sims[:,4]);digits=2));flush(stdout)
 end
 CSV.write(joinpath(OUT,"diagnostics.csv"),DataFrame(diagnosticrows))
 CSV.write(joinpath(OUT,"summary.csv"),DataFrame(summaries));CSV.write(joinpath(OUT,"draws.csv"),DataFrame(draws));CSV.write(joinpath(OUT,"support_distributions.csv"),DataFrame(histrows));CSV.write(joinpath(OUT,"checks.csv"),DataFrame(checks))
 for p in [joinpath(SOURCE,"occurrence_audit.csv"),joinpath(ROOT,"outputs/statistical description of local support/K_CO_SP_distributions/metadata_relationships.csv")]
  push!(paths,p);cp(p,joinpath(OUT,basename(p));force=true)
 end
 CSV.write(joinpath(OUT,"input_versions.csv"),DataFrame(path=relpath.(paths,ROOT),sha256=[bytes2hex(sha256(read(p))) for p in paths]))
end
function plots()
 s=CSV.read(joinpath(OUT,"summary.csv"),DataFrame);h=CSV.read(joinpath(OUT,"support_distributions.csv"),DataFrame)
 set_theme!(Theme(font="Arial",fontsize=18,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false)))
 f=Figure(size=(1550,650))
 for (j,(metric,label)) in enumerate(zip(["singleton_pct","gap_pp"],["Singletons among realised links (%)","Mean retention gap (percentage points)"]))
  ax=Axis(f[cld(j,2),mod1(j,2)],xlabel=label,yticks=(1:10,replace.(DATA,"_"=>" ")),yreversed=true)
  for (i,ds) in enumerate(DATA)
   r=only(eachrow(s[(s.dataset.==ds).&(s.metric.==metric),:]))
   lines!(ax,[r.null_lower,r.null_upper],[i,i],color=BLUE,linewidth=5)
   scatter!(ax,[r.null_mean],[i],color=BLUE,markersize=11)
   scatter!(ax,[r.empirical],[i],color=RED,marker=:diamond,markersize=14)
  end
 end
 Legend(f[2,1:2],[LineElement(color=BLUE,linewidth=5),MarkerElement(color=RED,marker=:diamond)],["Allocation mean and central 95% interval","Empirical"],orientation=:horizontal,framevisible=false)
 save(joinpath(OUT,"01_null_comparison.png"),f,px_per_unit=2)
 for (number,kind,ylabel) in [(2,"frequency","Proportion of realised links"),(3,"ccdf","P(support ≥ K | K > 0)")]
  fig=Figure(size=(1800,820))
  for (i,ds) in enumerate(DATA)
   v=h[h.dataset.==ds,:];ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Supporting sites (K)",ylabel=mod1(i,5)==1 ? ylabel : "")
   band!(ax,v.K,v[!,Symbol("null_lower_"*kind)],v[!,Symbol("null_upper_"*kind)],color=(BLUE,.2))
   lines!(ax,v.K,v[!,Symbol("null_mean_"*kind)],color=BLUE,linewidth=2.5)
   lines!(ax,v.K,v[!,Symbol("empirical_"*kind)],color=RED,linewidth=2.5)
   xmax=maximum(v.K[(v.empirical_frequency.>0).|(v.null_mean_frequency.>0)]);xlims!(ax,.8,xmax+.3);ylims!(ax,0,nothing)
  end
  Legend(fig[3,1:5],[LineElement(color=RED),LineElement(color=BLUE),PolyElement(color=(BLUE,.2))],["Empirical","Allocation mean","Pointwise central 95% interval"],orientation=:horizontal,framevisible=false)
  save(joinpath(OUT,lpad(number,2,'0')*"_support_"*kind*".png"),fig,px_per_unit=2)
 end
 open(joinpath(OUT,"FINDINGS.md"),"w") do io
  println(io,"# All CO null with constrained L\n\nAll positive co-occurring pairs are eligible. Observed M and L are fixed; realised pair identities vary. The target is uniform over opportunity-cell subsets conditional on both M and L, giving support vectors weights proportional to product binomial(C,K). This differs from the original known-link null, where pair identities are fixed.\n")
  println(io,"Four heat-bath Markov chains, 500 burn-in sweeps and 1,250 retained allocations per chain, separated by five sweeps. Each sweep attempts P pair-block updates. These are correlated MCMC samples, not independent exact draws. A finite-state enumeration checks detailed balance; split Rhat, ESS and Monte Carlo standard errors assess mixing for singleton percentage, gap, selected opportunity counts and overlap with observed links. Interpret intervals only alongside diagnostics.csv.\n")
  println(io,"| Dataset | Singleton % empirical / null [95%] | Gap pp empirical / null [95%] |\n|---|---:|---:|")
  for ds in DATA
   a=only(eachrow(s[(s.dataset.==ds).&(s.metric.=="singleton_pct"),:]));b=only(eachrow(s[(s.dataset.==ds).&(s.metric.=="gap_pp"),:]))
   println(io,@sprintf("| %s | %.2f / %.2f [%.2f, %.2f] | %.2f / %.2f [%.2f, %.2f] |",ds,a.empirical,a.null_mean,a.null_lower,a.null_upper,b.empirical,b.null_mean,b.null_lower,b.null_upper))
  end
  for metric in ["singleton_pct","gap_pp"]
   v=s[s.metric.==metric,:];println(io,"\n",metric,": empirical below / inside / above central 95% allocation interval: ",join([count(v.position.==p) for p in ["below","inside","above"]]," / ")," datasets.")
  end
  diag=CSV.read(joinpath(OUT,"diagnostics.csv"),DataFrame)
  println(io,"\nDiagnostics: maximum split Rhat = ",maximum(diag.split_Rhat),"; minimum ESS = ",minimum(diag.ESS),"; failing Rhat<1.01 and ESS≥1000: ",count(.!diag.passes)," / ",nrow(diag),".")
  println(io,"\nLink richness, zero-pair percentage and mean positive support M/L match observations by construction. Remaining differences concern allocation and identities, not total support or link richness. Gap is the arithmetic mean of exact expected regional retention minus 1-m/N over integer m/N≤0.8, including zero. Singleton percentages and distribution curves use positive-support links. Intervals describe conditional allocations, not confidence intervals or removal-sequence variation.\n\nThe same confirmed inputs are hash-checked and all C and K counts reconstructed as in free-L. Sampling/detection limitations, interaction-derived occurrence records and shared collections remain; see copied occurrence_audit.csv and metadata_relationships.csv. No mechanism is established by rejecting equal allocation. Existing analyses are preserved.\n\nRun with Julia using the RALSA project and this run.jl; append plot to redraw saved outputs.")
 end
end
if !("plot" in ARGS);main();end;plots()

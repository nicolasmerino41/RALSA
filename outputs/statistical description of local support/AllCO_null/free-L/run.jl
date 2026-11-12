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
function main()
 pairs=CSV.read(joinpath(SOURCE,"pairs.csv"),DataFrame)

 versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
 paths=[joinpath(SOURCE,"pairs.csv"),joinpath(SOURCE,"input_versions.csv"),@__FILE__]
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
  perm=collect(1:T);sim=zeros(Int,P)
  for b in 1:B
   fill!(sim,0)
   # Partial Fisher-Yates gives a uniform fixed-size subset of opportunity cells.
   for j in 1:M
    z=rand(rng,j:T);perm[j],perm[z]=perm[z],perm[j];sim[owners[perm[j]]]+=1
   end
   @assert sum(sim)==M && all(sim.<=C)
   l=count(>(0),sim)
   for k in sim;k>0 && (freq[b,k]+=1/l);end
   ccdf[b,:]=reverse(cumsum(reverse(freq[b,:])))
   sims[b,:]=[l,100(P-l)/P,100freq[b,1],100sum(freq[b,:].*weights)]
   push!(draws,(dataset=ds,replicate=b,L=l,zero_pct=sims[b,2],singleton_pct=sims[b,3],gap_pp=sims[b,4]))
  end
  exactL=sum(1-(c>T-M ? 0. : prod((T-M-j)/(T-j) for j in 0:c-1)) for c in C)
  @assert abs(mean(sims[:,1])-exactL)<=6std(sims[:,1])/sqrt(B)+1e-8
  push!(checks,(dataset=ds,N=N,pairs=P,total_opportunities=T,M=M,L=L,all_pairs_verified=true,known_support_verified=true,exact_expected_L=exactL,simulated_mean_L=mean(sims[:,1]),L_mcse=std(sims[:,1])/sqrt(B)))
  for (j,metric) in enumerate(["regional_links","zero_pct","singleton_pct","gap_pp"])
   v=sims[:,j];lo,hi=quantile(v,[.025,.975])
   push!(summaries,(dataset=ds,metric=metric,empirical=emp[j],null_mean=mean(v),null_lower=lo,null_upper=hi,difference=emp[j]-mean(v),position=emp[j]<lo ? "below" : emp[j]>hi ? "above" : "inside",mean_mcse=std(v)/sqrt(B),replicates=B))
  end
  obscc=reverse(cumsum(reverse(obsfreq)))
  for k in 1:N
   push!(histrows,(dataset=ds,K=k,empirical_frequency=obsfreq[k],null_mean_frequency=mean(freq[:,k]),null_lower_frequency=quantile(freq[:,k],.025),null_upper_frequency=quantile(freq[:,k],.975),empirical_ccdf=obscc[k],null_mean_ccdf=mean(ccdf[:,k]),null_lower_ccdf=quantile(ccdf[:,k],.025),null_upper_ccdf=quantile(ccdf[:,k],.975)))
  end
  println(ds," L observed/null: ",L," / ",round(mean(sims[:,1]);digits=1),"; gap: ",round(emp[4];digits=2)," / ",round(mean(sims[:,4]);digits=2));flush(stdout)
 end
 CSV.write(joinpath(OUT,"summary.csv"),DataFrame(summaries));CSV.write(joinpath(OUT,"draws.csv"),DataFrame(draws));CSV.write(joinpath(OUT,"support_distributions.csv"),DataFrame(histrows));CSV.write(joinpath(OUT,"checks.csv"),DataFrame(checks))
 for p in [joinpath(SOURCE,"occurrence_audit.csv"),joinpath(ROOT,"outputs/statistical description of local support/K_CO_SP_distributions/metadata_relationships.csv")]
  push!(paths,p);cp(p,joinpath(OUT,basename(p));force=true)
 end
 CSV.write(joinpath(OUT,"input_versions.csv"),DataFrame(path=relpath.(paths,ROOT),sha256=[bytes2hex(sha256(read(p))) for p in paths]))
end
function plots()
 s=CSV.read(joinpath(OUT,"summary.csv"),DataFrame);h=CSV.read(joinpath(OUT,"support_distributions.csv"),DataFrame)
 set_theme!(Theme(font="Arial",fontsize=18,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false)))
 f=Figure(size=(1550,1000))
 for (j,(metric,label)) in enumerate(zip(["regional_links","zero_pct","singleton_pct","gap_pp"],["Realised regional links","Pairs with K = 0 (%)","Singletons among realised links (%)","Mean retention gap (percentage points)"]))
  ax=Axis(f[cld(j,2),mod1(j,2)],xlabel=label,yticks=(1:10,replace.(DATA,"_"=>" ")),yreversed=true)
  for (i,ds) in enumerate(DATA)
   r=only(eachrow(s[(s.dataset.==ds).&(s.metric.==metric),:]))
   lines!(ax,[r.null_lower,r.null_upper],[i,i],color=BLUE,linewidth=5)
   scatter!(ax,[r.null_mean],[i],color=BLUE,markersize=11)
   scatter!(ax,[r.empirical],[i],color=RED,marker=:diamond,markersize=14)
  end
 end
 Legend(f[3,1:2],[LineElement(color=BLUE,linewidth=5),MarkerElement(color=RED,marker=:diamond)],["Allocation mean and central 95% interval","Empirical"],orientation=:horizontal,framevisible=false)
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
  println(io,"# All-co-occurrence null, free regional link richness\n")
  println(io,"Each allocation selects exactly M distinct pair-site opportunity cells uniformly without replacement from all positive co-occurring pairs. Regional link identities and richness are free. No parameter fitting. 5,000 independent draws per dataset, fixed seeds. Cell selection is uniform, not uniform over support-count vectors: vectors receive probability proportional to the product of binomial(C,K).\n")
  println(io,"K=0 fraction uses all positive-C pairs. Singleton fractions and support distributions use only K>0 links within each allocation. Gap uses that allocation's own original realised link count, with empirical and null local retention 1-m/N. Mean gaps average exact expected retention over integer removal levels m/N≤0.8, including zero. Thus gap differences include changed regional richness and allocation, not allocation alone.\n")
  println(io,"| Dataset | Regional links empirical / null | Zero pairs % empirical / null | Singleton % empirical / null | Gap pp empirical / null [95%] |\n|---|---:|---:|---:|---:|")
  for ds in DATA
   rows=[only(eachrow(s[(s.dataset.==ds).&(s.metric.==m),:])) for m in ["regional_links","zero_pct","singleton_pct","gap_pp"]]
   a,b,c,d=rows
   println(io,@sprintf("| %s | %.0f / %.1f | %.1f / %.1f | %.1f / %.1f | %.2f / %.2f [%.2f, %.2f] |",ds,a.empirical,a.null_mean,b.empirical,b.null_mean,c.empirical,c.null_mean,d.empirical,d.null_mean,d.null_lower,d.null_upper))
  end
  for metric in ["regional_links","zero_pct","singleton_pct","gap_pp"]
   v=s[s.metric.==metric,:];println(io,"\n",metric,": empirical below / inside / above allocation interval: ",join([count(v.position.==p) for p in ["below","inside","above"]]," / ")," datasets.")
  end
  println(io,"\nAcross all ten systems, empirical interactions occupy fewer pairs and repeat more often within the realised subset than equal random allocation across all opportunities predicts. Both the lower singleton fraction and higher gap lie outside the null central 95% intervals in every dataset. This differs from a known-link null because regional richness is no longer fixed. At fixed M, fewer realised links necessarily raise mean support M/L; the gap comparison cannot isolate allocation shape from this change in mean support. A future constrained-L test would separate those contributions.\n")
  println(io,"\nThese are conditional allocation intervals, not confidence intervals or removal-sequence variation. They test equal allocation among recorded opportunities, not ecological neutrality generally. Departures do not identify abundance, compatibility or detection mechanisms. No unseen species/pairs are estimated. Restored confirmed occurrences, interaction-derived presences, incomplete detection, sampling differences and shared collections retain their existing limitations (see copied audit and relationship tables). All opportunities were reconstructed from confirmed occupancy; known K was reconstructed from confirmed interaction records. Simulation mean link richness is checked against its exact hypergeometric expectation. Existing analyses are preserved.\n\nRun: julia --startup-file=no --project=RALSA_PATH \"PATH_TO_THIS_FOLDER/run.jl\".")
 end
end
if !("plot" in ARGS);main();end;plots()



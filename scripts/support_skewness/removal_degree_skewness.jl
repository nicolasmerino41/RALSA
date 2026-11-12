using CSV,DataFrames,Statistics,Random,SHA,CairoMakie,Printf
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const OUT=joinpath(ROOT,"outputs/statistical description of local support/support_skewness")
const SOURCE=joinpath(ROOT,"outputs/support_prediction_zero/confirmed")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const COLORS=["#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499","#661100"]
const B=2000
function skewness(x)
 length(x)<3 && return NaN
 d=x.-mean(x);v=mean(d.^2)
 v>0 ? mean(d.^3)/v^1.5 : NaN
end
function main()
 versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
 pairs=CSV.read(joinpath(SOURCE,"pairs.csv"),DataFrame)
 paths=[@__FILE__,joinpath(SOURCE,"pairs.csv"),joinpath(SOURCE,"input_versions.csv")]
 records=NamedTuple[]
 for (di,ds) in enumerate(DATA)
  p=pairs[(pairs.dataset.==ds).&(pairs.K.>0),:];N=only(unique(p.N));L=nrow(p)
  for kind in ["occupancy","interactions"]
   file=joinpath(SOURCE,"inputs",ds*"_"*kind*".csv");push!(paths,file)
   expected=only(versions.sha256[(versions.dataset.==ds).&(versions.kind.==kind)])
   @assert bytes2hex(sha256(read(file)))==expected
  end
  occ=CSV.read(paths[end-1],DataFrame;types=String);ints=CSV.read(paths[end],DataFrame;types=String)
  sites=sort(unique(occ.site));@assert length(sites)==N
  si=Dict(s=>i for (i,s) in enumerate(sites));li=Dict((r.consumer,r.resource)=>i for (i,r) in enumerate(eachrow(p)))
  support=[Int[] for _ in 1:L]
  for r in eachrow(ints);push!(support[li[(r.consumer,r.resource)]],si[r.site]);end
  foreach(unique!,support);@assert length.(support)==p.K
  consumers=unique(p.consumer);resources=unique(p.resource)
  ci=Dict(s=>i for (i,s) in enumerate(consumers));ri=Dict(s=>i for (i,s) in enumerate(resources))
  cids=[ci[s] for s in p.consumer];rids=[ri[s] for s in p.resource]
  rng=MersenneTwister(202611180+di)
  for b in 0:B
   rank=b==0 ? collect(1:N) : invperm(randperm(rng,N))
   last=[maximum(rank[s]) for s in support]
   levels=b==0 ? [0] : [20,40,60,80]
   for level in levels
    m=round(Int,N*level/100);cd=zeros(Int,length(consumers));rd=zeros(Int,length(resources))
    for link in 1:L
     if last[link]>m;cd[cids[link]]+=1;rd[rids[link]]+=1;end
    end
    for (guild,degrees) in [("consumer",cd),("resource",rd)]
     active=degrees[degrees.>0];skew=skewness(active)
     push!(records,(dataset=ds,guild=guild,removal_pct=level,removed_sites=m,N=N,actual_removal_pct=100m/N,replicate=b,connected_species=length(active),skewness=skew,valid=isfinite(skew)))
    end
   end
  end
  println(ds," complete");flush(stdout)
 end
 d=DataFrame(records);CSV.write(joinpath(OUT,"degree_skewness_removal_draws.csv"),d)
 rows=NamedTuple[]
 for g in groupby(d,[:dataset,:guild,:removal_pct])
  x=g.skewness[g.valid];n=length(x)
  push!(rows,(dataset=g.dataset[1],guild=g.guild[1],removal_pct=g.removal_pct[1],removed_sites=g.removed_sites[1],N=g.N[1],actual_removal_pct=g.actual_removal_pct[1],mean_skewness=n>0 ? mean(x) : NaN,mcse=n>1 ? std(x)/sqrt(n) : 0.,valid_iterations=n,undefined_iterations=nrow(g)-n,mean_connected_species=mean(g.connected_species)))
 end
 summary=DataFrame(rows);CSV.write(joinpath(OUT,"degree_skewness_removal_summary.csv"),summary)
 CSV.write(joinpath(OUT,"degree_skewness_removal_input_versions.csv"),DataFrame(path=relpath.(paths,ROOT),sha256=[bytes2hex(sha256(read(p))) for p in paths]))
 set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false)))
 for guild in ["consumer","resource"]
  f=Figure(size=(1700,780));v=summary[summary.guild.==guild,:];finite=filter(isfinite,v.mean_skewness);lo=min(0.,minimum(finite));hi=maximum(finite);pad=.08max(hi-lo,1.)
  for (i,ds) in enumerate(DATA)
   q=v[v.dataset.==ds,:];sort!(q,:removal_pct)
   ax=Axis(f[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Sites removed (%)",ylabel=mod1(i,5)==1 ? "Degree skewness" : "",xticks=([0,20,40,60,80],string.([0,20,40,60,80])))
   hlines!(ax,[0.],color=:gray80,linewidth=1)
   lines!(ax,q.removal_pct,q.mean_skewness;color=COLORS[i],linewidth=3)
   scatter!(ax,q.removal_pct,q.mean_skewness;color=COLORS[i],markersize=13)
   xlims!(ax,-4,84);ylims!(ax,lo-pad,hi+pad)
  end
  save(joinpath(OUT,"degree_skewness_removal_"*guild*"s.png"),f,px_per_unit=2)
 end
 open(joinpath(OUT,"degree_skewness_removal_README.md"),"w") do io
  println(io,"# Degree skewness across site removal\n\nSeparate consumer and resource figures, ten datasets each. Each species contributes its retained regional partner count once; only species with positive degree are included, matching the earlier degree-distribution analyses. These are degree distributions, not species occupancy or pair-support distributions.\n\nUniform random site permutations, 2,000 independent sequences per dataset; nested 20/40/60/80% removals within a sequence. Zero removal is calculated once. Site counts are rounded to the nearest integer, matching previous degree analyses. Plot positions label requested percentages; exact realised fractions are recorded in the CSV.\n\nSkewness is the population third central moment divided by variance to the power 3/2, matching support_skewness. Calculate it separately for each retained network, then average; no pooling of degrees across replicates. States with fewer than three connected species or zero variance are undefined and excluded explicitly. Means therefore condition on defined skewness; counts are reported. No uncertainty bands are plotted, as requested. Monte Carlo standard errors are saved.\n\nConfirmed occurrence/interaction inputs are hash-checked; site counts and link supports are verified. Existing sampling, detection, interaction-derived occurrence and shared-campaign limitations remain. No original files are replaced.\n")
  println(io,"| Guild | Datasets with lower mean skewness at 80% than baseline | Undefined draws | Largest Monte Carlo SE |\n|---|---:|---:|---:|")
  for guild in ["consumer","resource"]
   v=summary[summary.guild.==guild,:];n=count(ds->only(v.mean_skewness[(v.dataset.==ds).&(v.removal_pct.==80)])<only(v.mean_skewness[(v.dataset.==ds).&(v.removal_pct.==0)]),DATA)
   println(io,"| ",guild," | ",n," / 10 | ",sum(v.undefined_iterations)," | ",round(maximum(v.mcse);digits=4)," |")
  end
 end
end
main()

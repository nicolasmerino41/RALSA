using CSV, DataFrames, Statistics, CairoMakie, SHA, Printf
root=normpath(joinpath(@__DIR__,"../.."))
out=joinpath(root,"outputs/statistical description of local support/support_skewness")
mkpath(out)
source=joinpath(root,"outputs/statistical description of local support/K_CO_SP_distributions")
versions=NamedTuple[]
function load(p;kwargs...)
 push!(versions,(path=relpath(p,root),sha256=bytes2hex(sha256(read(p)))))
 CSV.read(p,DataFrame;kwargs...)
end
obs=load(joinpath(source,"observations.csv"))
datasets=unique(obs.dataset)
function moments(x,w=ones(length(x)))
 p=w./sum(w); mu=sum(p.*x); v=sum(p.*(x.-mu).^2)
 (mean=mu,sd=sqrt(v),skewness=v>0 ? sum(p.*(x.-mu).^3)/v^1.5 : NaN)
end
rows=NamedTuple[]; raw=NamedTuple[]
for ds in datasets
 occ=load(joinpath(root,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv");types=String)
 ints=load(joinpath(root,"data/processed",ds*"_interactions.csv");types=String)
 pres=Dict((g.trophic_level[1],g.species[1])=>Set(g.site) for g in groupby(occ,[:trophic_level,:species]))
 deg=Dict(k=>0 for k in keys(pres)); cs=Int[];ks=Int[]
 for g in groupby(ints,[:consumer,:resource])
  a=("consumer",g.consumer[1]);b=("resource",g.resource[1]);deg[a]+=1;deg[b]+=1
  c=intersect(pres[a],pres[b]);support=Set(g.site)
  @assert issubset(support,c)
  push!(cs,length(c));push!(ks,length(support))
 end
 saved=obs[obs.dataset.==ds,:]
 @assert sort(cs)==sort(saved.n_sites[saved.group.=="CO"])
 @assert sort(ks)==sort(saved.n_sites[saved.group.=="K"])
 keys_=collect(keys(pres));x=[length(pres[k]) for k in keys_];w=[deg[k] for k in keys_]
 @assert sort(x)==sort(saved.n_sites[saved.group.=="SP"])
 @assert sum(w)==2length(ks)
 for (stage,xx,ww) in [("Species weighted",x,w),("Co-occurrence",cs,ones(Int,length(cs))),("Local support",ks,ones(Int,length(ks))),("Species unweighted",x,ones(Int,length(x)))]
  z=moments(xx,ww)
  push!(rows,(dataset=ds,stage=stage,N=first(saved.N),n=length(xx),weight=sum(ww),mean=z.mean,sd=z.sd,skewness=z.skewness,singleton_fraction=sum(ww[xx.==1])/sum(ww)))
  for i in eachindex(xx);push!(raw,(dataset=ds,stage=stage,value=xx[i],weight=ww[i]));end
 end
end
summary=DataFrame(rows);CSV.write(joinpath(out,"skewness_summary.csv"),summary)
CSV.write(joinpath(out,"observations_and_weights.csv"),DataFrame(raw))
changes=NamedTuple[]
for ds in datasets
 g=summary[summary.dataset.==ds,:];s=g.skewness
 push!(changes,(dataset=ds,species_weighted=s[1],cooccurrence=s[2],local_support=s[3],delta_species_to_CO=s[2]-s[1],delta_CO_to_K=s[3]-s[2],both_increase=s[1]<s[2]<s[3]))
end
d=DataFrame(changes);CSV.write(joinpath(out,"stage_changes.csv"),d)
CSV.write(joinpath(out,"input_versions.csv"),DataFrame(versions))
colors=["#0072B2","#D55E00","#009E73"]
set_theme!(Theme(fontsize=18))
f=Figure(size=(1500,680))
lo=min(0.,minimum(d[:,[:species_weighted,:cooccurrence,:local_support]] |> Matrix))-.2
hi=maximum(d[:,[:species_weighted,:cooccurrence,:local_support]] |> Matrix)+.4
for (i,ds) in enumerate(datasets)
 ax=Axis(f[div(i-1,5)+1,mod(i-1,5)+1],title=replace(ds,"_"=>" "),xticks=(1:3,["Species","CO","Support"]),ylabel=mod(i-1,5)==0 ? "Skewness" : "")
 vals=[d.species_weighted[i],d.cooccurrence[i],d.local_support[i]]
 hlines!(ax,[0.],color=:gray80,linewidth=1)
 lines!(ax,1:3,vals,color=:gray40,linewidth=2)
 scatter!(ax,1:3,vals,color=colors,markersize=13)
 for j in 1:3;text!(ax,j,vals[j]+.10,text=@sprintf("%.2f",vals[j]),align=(:center,:bottom),fontsize=14);end
 ylims!(ax,lo,hi);xlims!(ax,.65,3.35)
 hidexdecorations!(ax,label=false,ticklabels=false,ticks=false);hideydecorations!(ax,label=false,ticklabels=false,ticks=false);hidespines!(ax,:t,:r)
end
save(joinpath(out,"01_skewness_stages.png"),f,px_per_unit=2)
f2=Figure(size=(1100,650))
for (j,col,label) in [(1,:delta_species_to_CO,"Species → co-occurrence"),(2,:delta_CO_to_K,"Co-occurrence → local support")]
 ax=Axis(f2[1,j],xlabel="Change in skewness",yticks=(1:10,replace.(datasets,"_"=>" ")),yreversed=true,title=label)
 v=d[!,col];barplot!(ax,1:10,v,direction=:x,color=[z>0 ? colors[1] : colors[2] for z in v]);vlines!(ax,[0.],color=:black,linewidth=1)
 hidespines!(ax,:t,:r)
end
save(joinpath(out,"02_skewness_changes.png"),f2,px_per_unit=2)
open(joinpath(out,"FINDINGS.md"),"w") do io
 println(io,"# Skewness across organisational stages\n")
 println(io,"Species occupancy is weighted by initial regional degree, pooling both guilds. CO and K follow the same known regional pairs. Confirmed occurrence-only species have zero degree weight. Unweighted species results are retained separately in the summary CSV.\n")
 println(io,"Skewness = weighted third central moment / weighted variance^(3/2). Weights sum to one; this is descriptive population-moment skewness, without sample-size correction. Positive values indicate a right tail. No prior scaling is required. Zero variance would return NaN.\n")
 println(io,"Species → CO increases: $(count(d.delta_species_to_CO.>0))/10. CO → support increases: $(count(d.delta_CO_to_K.>0))/10. Both steps increase: $(count(d.both_increase))/10.\n")
 println(io,"| Dataset | Species weighted | CO | Support | Both increase |\n|---|---:|---:|---:|---|")
 for r in eachrow(d);println(io,@sprintf("| %s | %.3f | %.3f | %.3f | %s |",r.dataset,r.species_weighted,r.cooccurrence,r.local_support,string(r.both_increase)));end
 println(io,"\nThis is a direct descriptive comparison, not a significance test or evidence of causal assembly. Shared species and sites prevent treating observations as independent replicates. Skewness is sensitive to rare high values; singleton fraction and mean remain in the CSV for context. Finite site count, sampling effort and detection affect the distributions. Degree weighting uses interaction information. Existing confirmed-occurrence exclusions and shared-collection limitations remain. Reconstructed C, K and unweighted occupancy distributions must exactly match the saved family-fitting observations before results are written. Other analyses are preserved.")
end
show(stdout,d;allcols=true);println()

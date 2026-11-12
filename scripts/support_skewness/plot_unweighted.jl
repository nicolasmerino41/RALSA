using CSV,DataFrames,Statistics,CairoMakie,Printf,SHA
root=normpath(joinpath(@__DIR__,"../.."))
out=joinpath(root,"outputs/statistical description of local support/support_skewness")
summary_path=joinpath(out,"skewness_summary.csv")
raw_path=joinpath(out,"observations_and_weights.csv")
summary=CSV.read(summary_path,DataFrame)
raw=CSV.read(raw_path,DataFrame)
datasets=unique(summary.dataset)
changes=NamedTuple[]
for ds in datasets
 vals=Float64[]
 for stage in ["Species unweighted","Co-occurrence","Local support"]
  x=raw.value[(raw.dataset.==ds).&(raw.stage.==stage)]
  mu=mean(x);v=mean((x.-mu).^2)
  skew=v>0 ? mean((x.-mu).^3)/v^1.5 : NaN
  saved=only(summary.skewness[(summary.dataset.==ds).&(summary.stage.==stage)])
  @assert isapprox(skew,saved;atol=1e-10)
  push!(vals,skew)
 end
 a,b,c=vals
 push!(changes,(dataset=ds,species_unweighted=a,cooccurrence=b,local_support=c,delta_species_to_CO=b-a,delta_CO_to_K=c-b,both_increase=a<b<c))
end
d=DataFrame(changes)
CSV.write(joinpath(out,"stage_changes_unweighted.csv"),d)
CSV.write(joinpath(out,"input_versions_unweighted.csv"),DataFrame([(path=relpath(p,root),sha256=bytes2hex(sha256(read(p)))) for p in [summary_path,raw_path]]))
colors=["#0072B2","#D55E00","#009E73"]
set_theme!(Theme(fontsize=18))
f=Figure(size=(1500,680))
lo=min(0.,minimum(d[:,[:species_unweighted,:cooccurrence,:local_support]] |> Matrix))-.2
hi=maximum(d[:,[:species_unweighted,:cooccurrence,:local_support]] |> Matrix)+.4
for (i,ds) in enumerate(datasets)
 ax=Axis(f[div(i-1,5)+1,mod(i-1,5)+1],title=replace(ds,"_"=>" "),xticks=(1:3,["Species","CO","Support"]),ylabel=mod(i-1,5)==0 ? "Skewness" : "")
 vals=[d.species_unweighted[i],d.cooccurrence[i],d.local_support[i]]
 hlines!(ax,[0.],color=:gray80,linewidth=1)
 lines!(ax,1:3,vals,color=:gray40,linewidth=2)
 scatter!(ax,1:3,vals,color=colors,markersize=13)
 for j in 1:3;text!(ax,j,vals[j]+.10,text=@sprintf("%.2f",vals[j]),align=(:center,:bottom),fontsize=14);end
 ylims!(ax,lo,hi);xlims!(ax,.65,3.35)
 hidexdecorations!(ax,label=false,ticklabels=false,ticks=false);hideydecorations!(ax,label=false,ticklabels=false,ticks=false);hidespines!(ax,:t,:r)
end
save(joinpath(out,"01_skewness_stages_unweighted.png"),f,px_per_unit=2)
f2=Figure(size=(1100,650))
for (j,col,label) in [(1,:delta_species_to_CO,"Species → co-occurrence"),(2,:delta_CO_to_K,"Co-occurrence → local support")]
 ax=Axis(f2[1,j],xlabel="Change in skewness",yticks=(1:10,replace.(datasets,"_"=>" ")),yreversed=true,title=label)
 v=d[!,col];barplot!(ax,1:10,v,direction=:x,color=[z>0 ? colors[1] : colors[2] for z in v]);vlines!(ax,[0.],color=:black,linewidth=1)
 hidespines!(ax,:t,:r)
end
save(joinpath(out,"02_skewness_changes_unweighted.png"),f2,px_per_unit=2)

open(joinpath(out,"FINDINGS_unweighted.md"),"w") do io
 println(io,"# Unweighted species occupancy to co-occurrence to support\n")
 println(io,"All recorded species from both guilds receive equal weight. CO and support use the same known regional interacting pairs as the weighted analysis. No other input changes. Skewness is the third central moment divided by variance to the power 3/2, with no sample-size correction. Values are recalculated from the saved observations and verified against the existing summary.\n")
 println(io,"Species to CO increases in $(count(d.delta_species_to_CO.>0))/10 datasets. CO to support increases in $(count(d.delta_CO_to_K.>0))/10. Both steps increase in $(count(d.both_increase))/10.\n")
 println(io,"| Dataset | Species unweighted | CO | Support | Both increase |\n|---|---:|---:|---:|---|")
 for r in eachrow(d);println(io,@sprintf("| %s | %.3f | %.3f | %.3f | %s |",r.dataset,r.species_unweighted,r.cooccurrence,r.local_support,string(r.both_increase)));end
 println(io,"\nDescriptive comparisons, not significance tests or a causal assembly test. Species and pairs are different observational units; skewness need not increase when per-pair counts decrease. Existing sampling, detection, occurrence and shared-collection limitations apply. Weighted plots and all previous analyses are preserved.")
end
show(stdout,d;allcols=true);println()

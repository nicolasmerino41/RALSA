using CSV,DataFrames,Statistics,Random,CairoMakie,LinearAlgebra,SHA
root=ARGS[1];dest=ARGS[2];mkpath(dest)
source=joinpath(root,"outputs/statistical description of local support/opportunity_theory/quality_control.csv")
d=CSV.read(source,DataFrame);datasets=unique(d.dataset)
function stats(x)
 v=sort(x);n=length(v);mu=mean(v)
 [std(v;corrected=false)/mu,2dot(1:n,v)/(n*sum(v))-(n+1)/n,quantile(v,.9)/mu]
end
rows=NamedTuple[];B=1999
for (di,ds) in enumerate(datasets)
 v=d[d.dataset.==ds,:];k=Float64.(v.K);c=Float64.(v.C);n=length(k);a=stats(k);b=stats(c);rng=MersenneTwister(7821+di);boot=zeros(B,3)
 for rep in 1:B
  ix=rand(rng,1:n,n);boot[rep,:]=stats(k[ix]).-stats(c[ix])
 end
 for (j,metric) in enumerate(["CV","Gini","p90_over_mean"])
  lo,hi=quantile(boot[:,j],[.025,.975]);push!(rows,(dataset=ds,metric=metric,support=a[j],cooccurrence=b[j],difference=a[j]-b[j],lower=lo,upper=hi))
 end
end
CSV.write(joinpath(dest,"shape_comparison.csv"),DataFrame(rows))
set_theme!(Theme(fontsize=18,Axis=(xgridvisible=false,ygridvisible=false,)))
fig=Figure(size=(1800,850))
for (i,ds) in enumerate(datasets)
 v=d[d.dataset.==ds,:]
 ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Sites / mean number of sites",ylabel=mod1(i,5)==1 ? "Proportion of pairs ≥ x" : "",xscale=log10,yscale=log10)
 for (j,field) in enumerate([:K,:C])
  values=v[!,field];x=sort(unique(values));y=[mean(values.>=t) for t in x]
  lines!(ax,x./mean(values),y;color=j==1 ? "#0072B2" : "#D55E00",linewidth=2.5)
 end
end
Legend(fig[3,1:5],[LineElement(color="#0072B2"),LineElement(color="#D55E00")],["Local support / mean support","Co-occurrences / mean co-occurrences"],orientation=:horizontal,framevisible=false)
save(joinpath(dest,"01_mean_scaled_ccdf.png"),fig;px_per_unit=2)
write(joinpath(dest,"README.md"),"# Quick empirical shape comparison\n\nSame known interacting pairs for K and C. Mean rescaling removes a multiplicative scale difference. Compare CV, Gini and 90th percentile divided by mean, which are unchanged by multiplication by a constant. Differences are K minus C. Paired link bootstrap, 1999 replicates, percentile 95% intervals. These are exploratory working intervals assuming independent pairs; shared species/sites mean they are not definitive ecological significance tests. No multiple-testing claim or distribution-family identification. Equal summary values do not establish equal distributions. Discrete values and finite support complicate collapse after scaling. No models fitted. Existing confirmed-occurrence flags remain unchanged.\n\nSource SHA256: "*bytes2hex(sha256(read(source)))*"\n")

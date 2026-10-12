module LocalSupportDescription
using CSV, DataFrames, Statistics, SHA, Printf, CairoMakie
include(joinpath(@__DIR__,"..","src","config.jl"))
const DEST=joinpath(OUT,"statistical description of local support")
const PALETTE=["#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499","#661100"]
const COLOR=Dict(zip(DATA,PALETTE))
mkpath(DEST)
summary=NamedTuple[];distributions=NamedTuple[];concentration=NamedTuple[];links=NamedTuple[];versions=NamedTuple[];sites=NamedTuple[]
expected=Dict(r.file=>r.sha256 for r in eachrow(CSV.read(joinpath(ROOT,"data/checksums.csv"),DataFrame)))
for ds in DATA
    tables=Dict{String,DataFrame}()
    for kind in ["interactions","occupancy","cooccurrences"]
        file=ds*"_"*kind*".csv";path=joinpath(INPUT,file);raw=joinpath(RAW,file)
        @assert bytes2hex(sha256(read(raw)))==expected[file]
        d=CSV.read(path,DataFrame;types=String)
        @assert !any(ismissing,Matrix(d)) && nrow(unique(d))==nrow(d)
        @assert all(x->!isempty(strip(x)),Matrix(d))
        @assert Set(Tuple.(eachrow(d)))==Set(Tuple.(eachrow(CSV.read(raw,DataFrame;types=String))))
        tables[kind]=d
        push!(versions,(dataset=ds,kind=kind,path=relpath(path,ROOT),sha256=bytes2hex(sha256(read(path))),raw_sha256=expected[file],rows=nrow(d)))
    end
    ints=tables["interactions"];occ=tables["occupancy"];co=tables["cooccurrences"]
    universe=sort(unique(vcat(ints.site,occ.site,co.site)));N=length(universe)
    for site in universe;push!(sites,(dataset=ds,site=site,has_interactions=site in ints.site));end
    p=combine(groupby(ints,[:consumer,:resource]),:site=>(x->length(unique(x)))=>:K)
    K=sort(p.K);L=length(K);T=sum(K);@assert all(1 .<= K .<= N) && T==nrow(ints)
    for r in eachrow(p);push!(links,(dataset=ds,consumer=r.consumer,resource=r.resource,K=r.K,N=N,relative_support=r.K/N));end
    x=vcat(0.,collect(1:L)./L);y=vcat(0.,cumsum(K)./T)
    gini=1-sum(diff(x).*(y[1:end-1].+y[2:end]))
    # Interpolation gives exactly 10% of links, with fractional boundary weight.
    pos=.9L;i=floor(Int,pos);bottom=i==L ? 1. : y[i+1]+(pos-i)*(y[i+2]-y[i+1])
    push!(summary,(dataset=ds,N=N,regional_links=L,total_support=T,mean=mean(K),median=median(K),
        standard_deviation=std(K),coefficient_of_variation=std(K)/mean(K),q25=quantile(K,.25),q75=quantile(K,.75),
        p90=quantile(K,.9),singleton_fraction=count(==(1),K)/L,maximum_support=maximum(K),mean_relative_support=mean(K)/N,
        support_gini=gini,top_10_percent_support_share=1-bottom))
    for k in 1:N
        push!(distributions,(dataset=ds,N=N,K=k,relative_support=k/N,link_count=count(==(k),K),frequency=count(==(k),K)/L,ccdf=count(>=(k),K)/L))
    end
    for j in eachindex(x);push!(concentration,(dataset=ds,cumulative_link_fraction=x[j],cumulative_support_share=y[j]));end
    @assert isapprox(last(y),1) && all(y .<= x .+ 1e-12)
end
s=DataFrame(summary);d=DataFrame(distributions);c=DataFrame(concentration)
for (name,rows) in [("summary",s),("distributions",d),("concentration",c),("link_support",DataFrame(links)),("input_versions",DataFrame(versions)),("site_universe",DataFrame(sites))]
    CSV.write(joinpath(DEST,name*".csv"),rows)
end
set_theme!(Theme(font="Arial",fontsize=20,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
function legend(slot)
    Legend(slot,[LineElement(color=COLOR[ds],linewidth=3) for ds in DATA],[DISPLAY[ds] for ds in DATA];nbanks=2,orientation=:horizontal,framevisible=false,labelsize=17)
end
savefig(f,name)=save(joinpath(DEST,name*".png"),f;px_per_unit=2)
f=Figure(size=(1750,820))
for (j,ds) in enumerate(DATA)
    v=d[(d.dataset.==ds).&(d.link_count.>0),:];row=only(eachrow(s[s.dataset.==ds,:]))
    ax=Axis(f[cld(j,5),mod1(j,5)],title=DISPLAY[ds],titlesize=19,xlabel="Supporting sites, K",ylabel=mod1(j,5)==1 ? "Proportion of regional links" : "")
    for r in eachrow(v);lines!(ax,[r.K,r.K],[0,r.frequency],color=COLOR[ds],linewidth=2);end
    scatter!(ax,v.K,v.frequency,color=COLOR[ds],markersize=6)
    xlims!(ax,.5,row.maximum_support+.5);ylims!(ax,0,.85)
end
savefig(f,"01_support_frequencies")
f=Figure(size=(1400,660))
for (j,label) in [(1,"Supporting sites, K"),(2,"Relative support, K / N")]
    ax=Axis(f[1,j],xlabel=label,ylabel="Proportion of links with support ≥ threshold",xscale=j==1 ? log10 : identity)
    for ds in DATA
        v=d[d.dataset.==ds,:];maxk=only(s.maximum_support[s.dataset.==ds]);v=v[v.K.<=min(maxk+1,first(v.N)),:]
        stairs!(ax,j==1 ? v.K : v.relative_support,v.ccdf;color=COLOR[ds],linewidth=2.5,step=:pre)
    end
    ylims!(ax,0,1.02)
    if j==1
        ticks=filter(x->x<=maximum(s.maximum_support),[1,2,5,10,20,50,100,200]);ax.xticks=(ticks,string.(ticks))
    else
        xlims!(ax,0,1)
    end
end
legend(f[2,1:2]);savefig(f,"02_support_ccdf_comparison")
f=Figure(size=(1400,700))
ax=Axis(f[1,1],xlabel="Cumulative proportion of links",ylabel="Cumulative share of total support",aspect=1)
lines!(ax,[0,1],[0,1],color=:gray65,linestyle=:dash)
for ds in DATA
    v=c[c.dataset.==ds,:];lines!(ax,v.cumulative_link_fraction,v.cumulative_support_share;color=COLOR[ds],linewidth=2.5)
end
xlims!(ax,0,1);ylims!(ax,0,1)
sorted=sort(s,:top_10_percent_support_share)
ax=Axis(f[1,2],yticks=(1:10,[DISPLAY[x] for x in sorted.dataset]),xlabel="Support share held by the top 10% of links")
vlines!(ax,[.1],color=:gray65,linestyle=:dash)
scatter!(ax,sorted.top_10_percent_support_share,1:10,color=[COLOR[x] for x in sorted.dataset],markersize=14)
for (i,r) in enumerate(eachrow(sorted));text!(ax,r.top_10_percent_support_share,i;text=@sprintf("%.1f%%",100r.top_10_percent_support_share),offset=(10,0),align=(:left,:center),fontsize=16);end
xlims!(ax,0,min(1,maximum(sorted.top_10_percent_support_share)+.12))
legend(f[2,1:2]);savefig(f,"03_support_concentration")
println(s)
end

module ParameterTransfer
using CSV, DataFrames, Statistics, SHA, Printf, CairoMakie
include(joinpath(@__DIR__,"..","src","config.jl"))
const SOURCE=joinpath(OUT,"statistical description of local support")
const FIT=joinpath(SOURCE,"two_parameter_fit")
const PREVIOUS=joinpath(SOURCE,"gap_reconstruction")
const DEST=joinpath(SOURCE,"parameter_transfer")
const PALETTE=["#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499","#661100"]
mkpath(DEST)
versions=CSV.read(joinpath(FIT,"input_versions.csv"),DataFrame)
for file in ["input_versions.csv","reconstruction_input_versions.csv"]
    for row in eachrow(CSV.read(joinpath(PREVIOUS,file),DataFrame))
        @assert bytes2hex(sha256(read(joinpath(ROOT,row.path))))==row.sha256 "Input changed since reconstruction"
    end
end
paths=[joinpath(FIT,"fit_summary.csv"),joinpath(SOURCE,"link_support.csv"),joinpath(SOURCE,"site_universe.csv"),joinpath(PREVIOUS,"curves.csv"),joinpath(PREVIOUS,"summary.csv")]
CSV.write(joinpath(DEST,"input_versions.csv"),versions)
CSV.write(joinpath(DEST,"transfer_input_versions.csv"),DataFrame(path=relpath.(paths,ROOT),sha256=[bytes2hex(sha256(read(p))) for p in paths]))
allparams=CSV.read(paths[1],DataFrame);params=allparams[allparams.model.=="two_parameter",:]
@assert nrow(params)==10 && Set(params.dataset)==Set(DATA) && all(isfinite,params.s) && all(isfinite,params.r)
links=CSV.read(paths[2],DataFrame);sites=CSV.read(paths[3],DataFrame;types=Dict(:site=>String))
previous=CSV.read(paths[4],DataFrame);oldsummary=CSV.read(paths[5],DataFrame)
function pmf(s,r,N)
    @assert N>=2 && 0<=s<=1 && 0<=r<=1
    t=r.^(0:N-2);vcat(s,(1-s).*t./sum(t))
end
function metrics(v,column)
    e=v[!,column];j=argmax(abs.(e))
    (mae=mean(abs.(e)),maxae=maximum(abs.(e)),bias=mean(e),worst=v.removal_fraction[j],signed_at_worst=e[j])
end
rows=NamedTuple[];summaries=NamedTuple[];training=NamedTuple[];probabilities=NamedTuple[]
for ds in DATA
    train=params[params.dataset.!=ds,:];@assert nrow(train)==9 && !(ds in train.dataset)
    sbar=mean(train.s);rbar=mean(train.r)
    N=nrow(sites[sites.dataset.==ds,:])
    # The held-out support is not used to calculate these probabilities.
    transferred=pmf(sbar,rbar,N)
    @assert isapprox(sum(transferred),1;atol=1e-12)
    for tr in eachrow(train);push!(training,(held_out_dataset=ds,training_dataset=tr.dataset,weight=1/9,s=tr.s,r=tr.r));end
    own=only(eachrow(params[params.dataset.==ds,:]));ownp=pmf(own.s,own.r,N)
    v=links[links.dataset.==ds,:];K=Int.(v.K);L=length(K);@assert only(unique(v.N))==N
    ints=CSV.read(joinpath(INPUT,ds*"_interactions.csv"),DataFrame;types=String)
    raw=combine(groupby(ints,[:consumer,:resource]),:site=>(x->length(unique(x)))=>:K)
    check=innerjoin(raw,v;on=[:consumer,:resource],makeunique=true)
    @assert nrow(check)==L==nrow(raw) && all(check.K.==check.K_1)
    emp=[count(==(k),K)/L for k in 1:N]
    for k in 1:N;push!(probabilities,(dataset=ds,N=N,K=k,transferred_probability=transferred[k],own_fit_probability=ownp[k],empirical_probability=emp[k]));end
    for m in 0:N
        survival=ones(N);logext=0.
        for k in 1:m
            logext+=log((m-k+1)/(N-k+1));survival[k]=-expm1(logext)
        end
        re=sum(emp.*survival);rt=sum(transferred.*survival);ro=sum(ownp.*survival);localret=1-m/N
        @assert min(re,rt,ro)>=localret-1e-12
        push!(rows,(dataset=ds,N=N,removed_sites=m,removal_fraction=m/N,in_main_range=5m<=4N,
            empirical_regional_retention=re,transferred_regional_retention=rt,own_fit_regional_retention=ro,local_support_retention=localret,
            empirical_gap=re-localret,transferred_gap=rt-localret,own_fit_gap=ro-localret,
            transferred_error_pp=100(rt-re),own_fit_error_pp=100(ro-re)))
    end
    alldf=DataFrame(filter(x->x.dataset==ds,rows));main=alldf[alldf.in_main_range,:]
    prior=sort(previous[previous.dataset.==ds,:],:removed_sites)
    @assert alldf.removed_sites==prior.removed_sites
    @assert maximum(abs.(alldf.empirical_gap.-prior.empirical_gap))<1e-12
    @assert maximum(abs.(alldf.own_fit_gap.-prior.fitted_gap))<1e-12
    for col in [:empirical_regional_retention,:transferred_regional_retention,:own_fit_regional_retention]
        @assert isapprox(first(alldf[!,col]),1) && last(alldf[!,col])==0 && maximum(diff(alldf[!,col]))<=1e-12
    end
    t=metrics(main,:transferred_error_pp);o=metrics(main,:own_fit_error_pp)
    @assert isapprox(o.mae,only(oldsummary.mean_absolute_error_pp[oldsummary.dataset.==ds]);atol=1e-10)
    push!(summaries,(dataset=ds,N=N,regional_links=L,training_datasets=join(train.dataset,";"),training_dataset_count=9,
        transferred_s=sbar,transferred_r=rbar,own_s=own.s,own_r=own.r,
        comparison_levels=nrow(main),maximum_removal_fraction=last(main.removal_fraction),
        transferred_mae_pp=t.mae,transferred_maximum_absolute_error_pp=t.maxae,transferred_signed_mean_error_pp=t.bias,
        transferred_worst_removal_fraction=t.worst,transferred_signed_error_at_worst_pp=t.signed_at_worst,
        own_fit_mae_pp=o.mae,own_fit_maximum_absolute_error_pp=o.maxae,own_fit_signed_mean_error_pp=o.bias,
        own_fit_worst_removal_fraction=o.worst,mae_deterioration_pp=t.mae-o.mae))
end
d=DataFrame(rows);s=DataFrame(summaries)
for (file,table) in [("summary",s),("curves",d),("training_parameters",DataFrame(training)),("support_probabilities",DataFrame(probabilities))]
    CSV.write(joinpath(DEST,file*".csv"),table)
end
aggregate=NamedTuple[]
for (label,mae,maxae,bias) in [("Transferred",s.transferred_mae_pp,s.transferred_maximum_absolute_error_pp,s.transferred_signed_mean_error_pp),
    ("Own fit (in-sample)",s.own_fit_mae_pp,s.own_fit_maximum_absolute_error_pp,s.own_fit_signed_mean_error_pp)]
    push!(aggregate,(comparison=label,datasets=10,mean_dataset_mae_pp=mean(mae),median_dataset_mae_pp=median(mae),
        mean_dataset_maximum_error_pp=mean(maxae),largest_dataset_maximum_error_pp=maximum(maxae),
        mean_dataset_signed_error_pp=mean(bias),mean_absolute_dataset_signed_error_pp=mean(abs.(bias)),
        mean_mae_deterioration_pp=label=="Transferred" ? mean(s.mae_deterioration_pp) : 0.))
end
CSV.write(joinpath(DEST,"across_dataset_summary.csv"),DataFrame(aggregate))
set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=15,yticklabelsize=15)))
main=d[d.in_main_range,:];ymax=5ceil(100maximum(vcat(main.empirical_gap,main.transferred_gap,main.own_fit_gap))/5)
f=Figure(size=(1850,850))
for (j,ds) in enumerate(DATA)
    v=main[main.dataset.==ds,:]
    ax=Axis(f[cld(j,5),mod1(j,5)],title=DISPLAY[ds],titlesize=19,xlabel="Sites removed (%)",ylabel=mod1(j,5)==1 ? "Regional − local retention (pp)" : "",xticks=0:20:80)
    lines!(ax,100v.removal_fraction,100v.empirical_gap;color=PALETTE[j],linewidth=3)
    lines!(ax,100v.removal_fraction,100v.transferred_gap;color=:black,linewidth=2.5,linestyle=:dash)
    lines!(ax,100v.removal_fraction,100v.own_fit_gap;color=:gray55,linewidth=2,linestyle=:dot)
    xlims!(ax,0,80);ylims!(ax,0,ymax)
end
Legend(f[3,1:5],[LineElement(color=:black,linewidth=3),LineElement(color=:black,linestyle=:dash),LineElement(color=:gray55,linestyle=:dot)],
    ["Empirical","Transferred parameters (held out)","Own-fit reconstruction (in-sample)"];orientation=:horizontal,framevisible=false)
save(joinpath(DEST,"01_parameter_transfer_gaps.png"),f;px_per_unit=2)
f=Figure(size=(1400,700))
for (j,field,label) in [(1,:transferred_mae_pp,"Mean absolute gap error (pp)"),(2,:transferred_signed_mean_error_pp,"Signed mean gap error (pp)")]
    ax=Axis(f[1,j],yticks=(1:10,[DISPLAY[x] for x in DATA]),yreversed=true,xlabel=label)
    ownfield=j==1 ? :own_fit_mae_pp : :own_fit_signed_mean_error_pp
    vlines!(ax,[0],color=:gray65,linestyle=:dash)
    for (i,ds) in enumerate(DATA)
        r=only(eachrow(s[s.dataset.==ds,:]));lines!(ax,[r[ownfield],r[field]],[i,i],color=PALETTE[i],linewidth=2)
        scatter!(ax,[r[field]],[i],color=PALETTE[i],marker=:circle,markersize=14)
        scatter!(ax,[r[ownfield]],[i],color=:gray40,marker=:diamond,markersize=11)
    end
end
Legend(f[2,1:2],[MarkerElement(color=:black,marker=:circle),MarkerElement(color=:gray40,marker=:diamond)],
    ["Transferred parameters","Own fit (in-sample)"];orientation=:horizontal,framevisible=false)
save(joinpath(DEST,"02_transfer_errors.png"),f;px_per_unit=2)
println(s);println(DataFrame(aggregate))
end

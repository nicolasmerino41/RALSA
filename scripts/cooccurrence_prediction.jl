module CooccurrencePrediction
using CSV, DataFrames, Statistics, LinearAlgebra, SHA, CairoMakie
include(joinpath(@__DIR__,"..","src","config.jl"))
const SOURCE=joinpath(OUT,"statistical description of local support")
const BASE=joinpath(SOURCE,"parameter_transfer")
const ZERO=joinpath(OUT,"support_prediction_zero")
const DEST=joinpath(SOURCE,"cooccurrence_prediction")
const COLORS=["#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499","#661100"]
mkpath(DEST)
logit(x)=log(x/(1-x))
logistic(x)=x>=0 ? 1/(1+exp(-x)) : exp(x)/(1+exp(x))
for file in [joinpath(BASE,"input_versions.csv"),joinpath(BASE,"transfer_input_versions.csv"),joinpath(ZERO,"confirmed/input_versions.csv")]
    for row in eachrow(CSV.read(file,DataFrame))
        path=hasproperty(row,:path) ? joinpath(ROOT,row.path) : joinpath(ZERO,"confirmed/inputs",row.dataset*"_"*row.kind*".csv")
        @assert bytes2hex(sha256(read(path)))==row.sha256 "Input version changed"
    end
end
links=CSV.read(joinpath(SOURCE,"link_support.csv"),DataFrame)
sites=CSV.read(joinpath(SOURCE,"site_universe.csv"),DataFrame;types=Dict(:site=>String))
par=CSV.read(joinpath(SOURCE,"two_parameter_fit/fit_summary.csv"),DataFrame);par=par[par.model.=="two_parameter",:]
@assert nrow(par)==10 && all(0 .< par.s .< 1) && all(0 .< par.r .< 1)
audit=NamedTuple[];features=NamedTuple[];versions=NamedTuple[];occupancy_audit=NamedTuple[]
for ds in DATA
    path=joinpath(ZERO,"confirmed/inputs",ds*"_occupancy.csv");occ=CSV.read(path,DataFrame;types=String)
    @assert nrow(occ)==nrow(unique(occ)) && !any(ismissing,Matrix(occ))
    universe=Set(sites.site[sites.dataset.==ds]);N=length(universe);@assert Set(occ.site)⊆universe
    presence=Dict((g.trophic_level[1],g.species[1])=>Set(g.site) for g in groupby(occ,[:trophic_level,:species]))
    ids=select(links[links.dataset.==ds,:],:consumer,:resource);Cs=Int[]
    for pair in eachrow(ids)
        C=length(intersect(get(presence,("consumer",pair.consumer),Set{String}()),get(presence,("resource",pair.resource),Set{String}())))
        push!(Cs,C)
        K=only(links.K[(links.dataset.==ds).&(links.consumer.==pair.consumer).&(links.resource.==pair.resource)])
        push!(audit,(dataset=ds,consumer=pair.consumer,resource=pair.resource,C=C,K_qc_only=K,N=N,valid=1<=K<=C<=N))
    end
    if all(>(0),Cs)
        push!(features,(dataset=ds,N=N,regional_links=length(Cs),x1=mean(Cs.==1),x2=mean(log.(Cs))))
    end
    original=CSV.read(joinpath(INPUT,ds*"_occupancy.csv"),DataFrame;types=String)
    @assert Set(Tuple.(eachrow(original)))⊆Set(Tuple.(eachrow(occ)))
    push!(occupancy_audit,(dataset=ds,confirmed_added_presences=nrow(occ)-nrow(original),occupancy_version="confirmed restorations; uncertain additions excluded",independently_sampled_both_guilds=false))
    push!(versions,(path=relpath(path,ROOT),sha256=bytes2hex(sha256(read(path)))))
end
qc=DataFrame(audit);CSV.write(joinpath(DEST,"link_cooccurrence_qc.csv"),qc)
CSV.write(joinpath(DEST,"unresolved_inconsistencies.csv"),qc[.!qc.valid,:])
CSV.write(joinpath(DEST,"occupancy_audit.csv"),DataFrame(occupancy_audit))
@assert all(qc.valid) "K/C/N inconsistency: see unresolved_inconsistencies.csv; no links altered or dropped"
predictors=DataFrame(features);CSV.write(joinpath(DEST,"predictors.csv"),predictors)
for path in [joinpath(SOURCE,"link_support.csv"),joinpath(SOURCE,"site_universe.csv"),joinpath(SOURCE,"two_parameter_fit/fit_summary.csv"),joinpath(BASE,"curves.csv"),joinpath(BASE,"summary.csv"),joinpath(ZERO,"presence_record_provenance.csv"),joinpath(ZERO,"source_checksums.csv")]
    push!(versions,(path=relpath(path,ROOT),sha256=bytes2hex(sha256(read(path)))))
end
CSV.write(joinpath(DEST,"input_versions.csv"),DataFrame(versions))
dataset=innerjoin(predictors,select(par,:dataset,:s,:r);on=:dataset)
baseline=CSV.read(joinpath(BASE,"curves.csv"),DataFrame);basemetrics=CSV.read(joinpath(BASE,"summary.csv"),DataFrame)
rows=NamedTuple[];metrics=NamedTuple[];fits=NamedTuple[];parameter_predictions=NamedTuple[]
for ds in DATA
    train=dataset[dataset.dataset.!=ds,:];target=only(eachrow(predictors[predictors.dataset.==ds,:]))
    @assert nrow(train)==9 && !(ds in train.dataset)
    predictions=Float64[]
    for (response,predictor) in [(:s,:x1),(:r,:x2)]
        X=hcat(ones(9),train[!,predictor]);@assert rank(X)==2
        coef=X\logit.(train[!,response]);pred=logistic(coef[1]+coef[2]*target[predictor]);push!(predictions,pred)
        push!(fits,(held_out_dataset=ds,response=String(response),predictor=String(predictor),intercept=coef[1],slope=coef[2],
            training_datasets=join(train.dataset,";"),training_weight=1/9,predictor_min=minimum(train[!,predictor]),predictor_max=maximum(train[!,predictor]),
            target_predictor=target[predictor],extrapolation=!(minimum(train[!,predictor])<=target[predictor]<=maximum(train[!,predictor]))))
    end
    shat,rhat=predictions;N=target.N;tail=rhat.^(0:N-2);p=vcat(shat,(1-shat).*tail./sum(tail));@assert isapprox(sum(p),1;atol=1e-12)
    # Held-out fitted parameters enter evaluation only after prediction.
    own=only(eachrow(par[par.dataset.==ds,:]));b=only(eachrow(basemetrics[basemetrics.dataset.==ds,:]))
    @assert isapprox(b.transferred_s,mean(train.s)) && isapprox(b.transferred_r,mean(train.r))
    push!(parameter_predictions,(dataset=ds,x1=target.x1,x2=target.x2,observed_s=own.s,predicted_s=shat,baseline_s=b.transferred_s,
        observed_r=own.r,predicted_r=rhat,baseline_r=b.transferred_r))
    prior=sort(baseline[baseline.dataset.==ds,:],:removed_sites);@assert prior.removed_sites==collect(0:N)
    for z in eachrow(prior)
        m=z.removed_sites;survival=ones(N);logext=0.
        for k in 1:m;logext+=log((m-k+1)/(N-k+1));survival[k]=-expm1(logext);end
        retention=sum(p.*survival);gap=retention-(1-m/N)
        push!(rows,(dataset=ds,N=N,removed_sites=m,removal_fraction=m/N,in_main_range=5m<=4N,local_support_retention=1-m/N,
            empirical_gap=z.empirical_gap,cooccurrence_gap=gap,baseline_gap=z.transferred_gap,own_fit_gap=z.own_fit_gap,
            predicted_regional_retention=retention,cooccurrence_error_pp=100(gap-z.empirical_gap),baseline_error_pp=z.transferred_error_pp,own_fit_error_pp=z.own_fit_error_pp))
    end
    v=DataFrame(filter(x->x.dataset==ds && x.in_main_range,rows));mean_gap=100mean(v.empirical_gap);@assert mean_gap>0
    base_mae=mean(abs.(v.baseline_error_pp));@assert isapprox(base_mae,b.transferred_mae_pp;atol=1e-10)
    for (model,field) in [("Cooccurrence",:cooccurrence_error_pp),("Common baseline",:baseline_error_pp),("Own fit (in-sample)",:own_fit_error_pp)]
        err=v[!,field];mae=mean(abs.(err));worst=argmax(abs.(err))
        push!(metrics,(dataset=ds,model=model,N=N,comparison_levels=nrow(v),mean_empirical_gap_pp=mean_gap,mae_pp=mae,
            maximum_absolute_error_pp=maximum(abs.(err)),signed_mean_error_pp=mean(err),worst_removal_fraction=v.removal_fraction[worst],
            change_in_mae_vs_baseline_pp=mae-base_mae,relative_error_percent=100mae/mean_gap,improved_vs_baseline=mae<base_mae-1e-10))
    end
end
d=DataFrame(rows);s=DataFrame(metrics);pred=DataFrame(parameter_predictions)
for (name,table) in [("curves",d),("performance",s),("parameter_predictions",pred),("regressions",DataFrame(fits))]
    CSV.write(joinpath(DEST,name*".csv"),table)
end
aggregate=combine(groupby(s,:model),:mae_pp=>mean=>:mean_mae_pp,:mae_pp=>median=>:median_mae_pp,
    :maximum_absolute_error_pp=>mean=>:mean_maximum_error_pp,:signed_mean_error_pp=>mean=>:mean_signed_error_pp,
    :relative_error_percent=>mean=>:mean_relative_error_percent,:relative_error_percent=>median=>:median_relative_error_percent,
    :change_in_mae_vs_baseline_pp=>mean=>:mean_change_in_mae_pp,:improved_vs_baseline=>sum=>:datasets_improved)
CSV.write(joinpath(DEST,"across_dataset_summary.csv"),aggregate)
set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=15,yticklabelsize=15)))
f=Figure(size=(1250,700))
for (j,response) in enumerate(["s","r"])
    ax=Axis(f[1,j],xlabel="Observed "*response,ylabel="Held-out predicted "*response,aspect=1)
    lines!(ax,[0,1],[0,1],color=:gray60,linestyle=:dash)
    for (i,ds) in enumerate(DATA)
        z=only(eachrow(pred[pred.dataset.==ds,:]));x=z[Symbol("observed_"*response)]
        scatter!(ax,[x],[z[Symbol("predicted_"*response)]],color=COLORS[i],markersize=13)
        scatter!(ax,[x],[z[Symbol("baseline_"*response)]],color=COLORS[i],marker=:xcross,markersize=12)
    end
    xlims!(ax,.25,.95);ylims!(ax,.25,.95)
end
Legend(f[2,1:2],[MarkerElement(color=:black,marker=:circle),MarkerElement(color=:black,marker=:xcross)],["Co-occurrence regressions","Common parameter baseline"];orientation=:horizontal,framevisible=false)
Legend(f[3,1:2],[MarkerElement(color=c,marker=:circle) for c in COLORS],[DISPLAY[x] for x in DATA];orientation=:horizontal,nbanks=2,framevisible=false)
save(joinpath(DEST,"01_parameter_predictions.png"),f;px_per_unit=2)
main=d[d.in_main_range,:];f=Figure(size=(1850,870));ymax=5ceil(100maximum(vcat(main.empirical_gap,main.cooccurrence_gap,main.baseline_gap,main.own_fit_gap))/5)
for (j,ds) in enumerate(DATA)
    v=main[main.dataset.==ds,:];ax=Axis(f[cld(j,5),mod1(j,5)],title=DISPLAY[ds],titlesize=19,xlabel="Sites removed (%)",ylabel=mod1(j,5)==1 ? "Regional − local retention (pp)" : "",xticks=0:20:80)
    for (field,color,style,width) in [(:empirical_gap,COLORS[j],:solid,3),(:cooccurrence_gap,:black,:dash,2.5),(:baseline_gap,:gray50,:dashdot,2),(:own_fit_gap,:gray65,:dot,2)]
        lines!(ax,100v.removal_fraction,100v[!,field],color=color,linestyle=style,linewidth=width)
    end
    xlims!(ax,0,80);ylims!(ax,0,ymax)
end
Legend(f[3,1:5],[LineElement(color=:black),LineElement(color=:black,linestyle=:dash),LineElement(color=:gray50,linestyle=:dashdot),LineElement(color=:gray65,linestyle=:dot)],
    ["Empirical","Co-occurrence regressions (held out)","Common baseline (held out)","Own fit (in-sample)"];orientation=:horizontal,framevisible=false)
save(joinpath(DEST,"02_gap_predictions.png"),f;px_per_unit=2)
f=Figure(size=(1200,690))
for (j,field,label) in [(1,:mae_pp,"Gap MAE (pp)"),(2,:relative_error_percent,"MAE / mean empirical gap (%)")]
    ax=Axis(f[1,j],yticks=(1:10,[DISPLAY[x] for x in DATA]),yreversed=true,xlabel=label)
    for (i,ds) in enumerate(DATA)
        v=s[s.dataset.==ds,:];b=only(v[v.model.=="Common baseline",field]);p=only(v[v.model.=="Cooccurrence",field])
        lines!(ax,[b,p],[i,i],color=COLORS[i],linewidth=2);scatter!(ax,[b],[i],color=:gray45,marker=:diamond,markersize=11);scatter!(ax,[p],[i],color=COLORS[i],markersize=13)
    end
end
Legend(f[2,1:2],[MarkerElement(color=:black,marker=:circle),MarkerElement(color=:gray45,marker=:diamond)],["Co-occurrence regressions","Common baseline"];orientation=:horizontal,framevisible=false)
save(joinpath(DEST,"03_performance_comparison.png"),f;px_per_unit=2)
println(aggregate);println(s[s.model.=="Cooccurrence",:])
end

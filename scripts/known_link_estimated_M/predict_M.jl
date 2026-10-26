module KnownLinkEstimatedM
using CSV, DataFrames, Statistics, LinearAlgebra, SHA, CairoMakie, Printf
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"../..")) : abspath(ARGS[1])
const SOURCE=joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory")
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(SOURCE,"known_link_estimated_M")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP2","Olot","Montseny"]
const COLORS=["#332288","#88CCEE","#44AA99","#117733","#999933","#882255","#AA4499","#661100"]
logz(x)=log(-expm1(x))
lse(x)=maximum(x)+log(sum(exp.(x.-maximum(x))))
logistic(x)=1/(1+exp(-x))
include(joinpath(ROOT,"scripts/beta_pair_conditioned_M/conditional_support.jl"))
mkpath(DEST)
function regression(x,y,intercept)
    b=intercept ? dot(x.-mean(x),y.-mean(y))/sum(abs2,x.-mean(x)) : dot(x,y)/dot(x,x)
    a=intercept ? mean(y)-b*mean(x) : 0.
    a,b
end
function main()
    versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
    for r in eachrow(versions);@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256 "Changed input: $(r.path)";end
    used=NamedTuple[]
    function readtable(path;kwargs...)
        push!(used,(path=relpath(path,ROOT),sha256=bytes2hex(sha256(read(path)))))
        CSV.read(path,DataFrame;kwargs...)
    end
    pairs=readtable(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/pairs.csv"))
    qc=readtable(joinpath(SOURCE,"quality_control.csv"))
    fits=readtable(joinpath(SOURCE,"fitted_parameters.csv"))
    old=readtable(joinpath(SOURCE,"conditioned_M/summary.csv"))
    oldcurves=readtable(joinpath(SOURCE,"conditioned_M/curves.csv"))
    sites=readtable(joinpath(ROOT,"outputs/statistical description of local support/site_universe.csv");types=Dict(:site=>String))
    @assert Set(pairs.dataset)==Set(DATA)
    totals=NamedTuple[];checks=NamedTuple[]
    for ds in DATA
        v=pairs[pairs.dataset.==ds,:];g=qc[qc.dataset.==ds,:];N=only(unique(g.N))
        @assert all(v.N.==N) && N==count(==(ds),sites.dataset)
        @assert all(0 .<=v.K .<=v.n .<=N)
        positive=v[v.K.>0,:]
        @assert Dict((r.consumer,r.resource)=>(r.K,r.n) for r in eachrow(positive))==Dict((r.consumer,r.resource)=>(r.K,r.C) for r in eachrow(g))
        @assert nrow(unique(v,[:consumer,:resource]))==nrow(v)
        occ=readtable(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv");types=String)
        pres=Dict((z.trophic_level[1],z.species[1])=>Set(z.site) for z in groupby(occ,[:trophic_level,:species]))
        @assert all(r.n==length(intersect(pres[("consumer",r.consumer)],pres[("resource",r.resource)])) for r in eachrow(v))
        total_from_sites=sum(length(unique(z.species[z.trophic_level.=="consumer"]))*length(unique(z.species[z.trophic_level.=="resource"])) for z in groupby(occ,:site))
        @assert sum(v.n)==total_from_sites
        M=sum(v.K);L=nrow(g);T=sum(v.n);capacity=sum(g.C)
        @assert all(old.M[old.dataset.==ds].==M) && all(old.L[old.dataset.==ds].==L)
        push!(totals,(dataset=ds,N=N,L=L,M=M,total_cooccurrences=T,known_link_opportunities=capacity,zero_K_pairs=count(v.K.==0),zero_K_cooccurrences=sum(v.n[v.K.==0]),observed_realisation=M/T))
        push!(checks,(dataset=ds,check="All-pair opportunities match occupancy; positive pairs and M match original",maximum_error=0.))
    end
    total=DataFrame(totals);CSV.write(joinpath(DEST,"dataset_totals.csv"),total)
    x=Float64.(total.known_link_opportunities);y=Float64.(total.M)
    full=NamedTuple[]
    for intercept in [false,true]
        a,b=regression(x,y,intercept)
        push!(full,(model=intercept ? "Intercept" : "Proportional",intercept=a,slope=b,slope_percent=100b,pearson_r=cor(x,y),scope="All ten datasets; descriptive only"))
    end
    CSV.write(joinpath(DEST,"full_data_regressions.csv"),DataFrame(full))
    rows=NamedTuple[];curves=NamedTuple[];regrows=NamedTuple[];freqrows=NamedTuple[]
    for (di,ds) in enumerate(DATA)
        t=total[di,:];g=qc[qc.dataset.==ds,:];C=Int.(g.C);K=Int.(g.K);N=t.N;L=t.L
        H=gapmatrix(N);idx=1:fld(4N,5)+1;emp=100vec(mean(H[:,K];dims=2));training=filter(!=(ds),DATA);train=total.dataset.!=ds
        for (scenario,intercept) in [("Proportional",false),("Intercept",true)]
            a,b=regression(x[train],y[train],intercept);raw=a+b*t.known_link_opportunities
            rounded=round(Int,raw);M=clamp(rounded,L,t.known_link_opportunities)
            flag=rounded<L ? "Below known-link minimum; projected to L" : rounded>t.known_link_opportunities ? "Above known-link capacity; projected to sum(C)" : "Feasible"
            push!(regrows,(dataset=ds,scenario=scenario,intercept=a,slope=b,slope_percent=100b,training_datasets=join(training,";"),predictor_cooccurrences=t.known_link_opportunities,predicted_M_raw=raw,predicted_M_rounded=rounded,predicted_M_used=M,feasibility=flag,observed_M=t.M,relative_M_error_percent=100(raw/t.M-1)))
            for model in ["Pair","Common"]
                par=only(eachrow(fits[(fits.dataset.==ds).&(fits.scope.=="LODO").&(fits.model.==model),:]))
                @assert Set(split(par.training_datasets,";"))==Set(training)
                if M==L
                    f=zeros(N);f[1]=1.
                elseif M==t.known_link_opportunities
                    f=[count(==(k),C)/L for k in 1:N]
                else
                    D,P=dp(pairweights(C,par.a,par.h),M);f=exactfreq(M,D,P,N)
                end
                @assert isapprox(sum(f),1;atol=1e-8) && isapprox(dot(1:N,f),M/L;atol=1e-8)
                pred=100H*f;err=pred[idx].-emp[idx]
                reference=only(eachrow(old[(old.dataset.==ds).&(old.model.==model),:]))
                push!(rows,(dataset=ds,scenario=scenario,model=model,N=N,L=L,observed_M=t.M,predicted_M_used=M,feasibility=flag,a=par.a,h=par.h,training_datasets=par.training_datasets,mean_empirical_gap_pp=mean(emp[idx]),mean_predicted_gap_pp=mean(pred[idx]),mae_pp=mean(abs.(err)),maximum_error_pp=maximum(abs.(err)),signed_error_pp=mean(err),removal_at_maximum_error=(argmax(abs.(err))-1)/N,relative_gap_error_percent=100mean(abs.(err))/mean(emp[idx]),observed_M_mae_pp=reference.mae_pp,mae_change_pp=mean(abs.(err))-reference.mae_pp))
                ref=oldcurves[(oldcurves.dataset.==ds).&(oldcurves.model.==model),:];@assert ref.m==collect(0:N)
                @assert maximum(abs.(ref.empirical_gap_pp.-emp))<1e-9
                for m in 0:N;push!(curves,(dataset=ds,scenario=scenario,model=model,m=m,removal_fraction=m/N,empirical_gap_pp=emp[m+1],predicted_gap_pp=pred[m+1],observed_M_gap_pp=ref.predicted_gap_pp[m+1]));end
                for k in 1:N;push!(freqrows,(dataset=ds,scenario=scenario,model=model,K=k,predicted_frequency=f[k]));end
                push!(checks,(dataset=ds,check="$scenario $model: exact M, normalized frequencies, original empirical curve",maximum_error=abs(dot(1:N,f)-M/L)))
                println(ds," ",scenario," ",model,": M=",M," (observed ",t.M,"); MAE=",round(mean(abs.(err));digits=3));flush(stdout)
            end
        end
        CSV.write(joinpath(DEST,"summary.csv"),DataFrame(rows));CSV.write(joinpath(DEST,"curves.csv"),DataFrame(curves));CSV.write(joinpath(DEST,"regression_folds.csv"),DataFrame(regrows))
    end
    tab=DataFrame(rows);agg=NamedTuple[]
    for scenario in ["Proportional","Intercept"],model in ["Pair","Common"]
        v=tab[(tab.scenario.==scenario).&(tab.model.==model),:]
        push!(agg,(scenario=scenario,model=model,mean_mae_pp=mean(v.mae_pp),median_mae_pp=median(v.mae_pp),mean_signed_error_pp=mean(v.signed_error_pp),mean_relative_error_percent=mean(v.relative_gap_error_percent),observed_M_mean_mae_pp=mean(v.observed_M_mae_pp),mean_mae_change_pp=mean(v.mae_change_pp),datasets_improved=count(v.mae_change_pp.<0),datasets_projected=count(v.feasibility.!="Feasible")))
    end
    CSV.write(joinpath(DEST,"across_dataset_summary.csv"),DataFrame(agg));CSV.write(joinpath(DEST,"support_frequencies.csv"),DataFrame(freqrows));CSV.write(joinpath(DEST,"result_checks.csv"),DataFrame(checks))
    CSV.write(joinpath(DEST,"input_versions.csv"),versions)
    helper=joinpath(ROOT,"scripts/beta_pair_conditioned_M/conditional_support.jl");push!(used,(path=relpath(helper,ROOT),sha256=bytes2hex(sha256(read(helper)))))
    CSV.write(joinpath(DEST,"source_versions.csv"),DataFrame(used))
    CSV.write(joinpath(DEST,"parameter_estimates_used.csv"),fits[(fits.scope.=="LODO").&in.(fits.model,Ref(["Pair","Common"])),:])
    for file in ["metadata_relationships.csv","metadata_source_versions.csv"];cp(joinpath(SOURCE,file),joinpath(DEST,file);force=true);end
end

include("compare.jl")
if !(length(ARGS)>=3 && ARGS[3]=="plot");main();compare();end
plots();report()
end

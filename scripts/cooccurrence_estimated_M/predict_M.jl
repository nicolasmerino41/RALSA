module CooccurrenceEstimatedM
using CSV, DataFrames, Statistics, LinearAlgebra, SHA, CairoMakie, Printf
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"../..")) : abspath(ARGS[1])
const SOURCE=joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory")
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(SOURCE,"cooccurrence_estimated_M")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const COLORS=["#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499","#661100"]
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
    x=Float64.(total.total_cooccurrences);y=Float64.(total.M)
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
            a,b=regression(x[train],y[train],intercept);raw=a+b*t.total_cooccurrences
            rounded=round(Int,raw);M=clamp(rounded,L,t.known_link_opportunities)
            flag=rounded<L ? "Below known-link minimum; projected to L" : rounded>t.known_link_opportunities ? "Above known-link capacity; projected to sum(C)" : "Feasible"
            push!(regrows,(dataset=ds,scenario=scenario,intercept=a,slope=b,slope_percent=100b,training_datasets=join(training,";"),total_cooccurrences=t.total_cooccurrences,predicted_M_raw=raw,predicted_M_rounded=rounded,predicted_M_used=M,feasibility=flag,observed_M=t.M,relative_M_error_percent=100(raw/t.M-1)))
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
function plots()
    t=CSV.read(joinpath(DEST,"dataset_totals.csv"),DataFrame);r=CSV.read(joinpath(DEST,"regression_folds.csv"),DataFrame);f=CSV.read(joinpath(DEST,"full_data_regressions.csv"),DataFrame);c=CSV.read(joinpath(DEST,"curves.csv"),DataFrame);s=CSV.read(joinpath(DEST,"summary.csv"),DataFrame)
    set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
    fig=Figure(size=(1420,650));ax=Axis(fig[1,1];xlabel="Total co-occurrences (including K = 0)",ylabel="Total local support, M");bx=Axis(fig[1,2];xlabel="Observed total support, M",ylabel="Held-out estimated support, M")
    xx=range(0,maximum(t.total_cooccurrences)*1.05;length=100)
    for z in eachrow(f);lines!(ax,xx,z.intercept .+z.slope.*xx;color=z.model=="Proportional" ? :black : :gray,linestyle=z.model=="Proportional" ? :solid : :dash,linewidth=2);end
    upper=max(maximum(t.M),maximum(r.predicted_M_raw))*1.1;lines!(bx,[0,upper],[0,upper];color=:gray,linestyle=:dash)
    for (i,ds) in enumerate(DATA)
        scatter!(ax,[t.total_cooccurrences[i]],[t.M[i]];color=COLORS[i],markersize=12)
        for (scenario,marker) in [("Proportional",:circle),("Intercept",:utriangle)]
            z=only(eachrow(r[(r.dataset.==ds).&(r.scenario.==scenario),:]));scatter!(bx,[z.observed_M],[z.predicted_M_raw];color=COLORS[i],marker=marker,markersize=12)
        end
    end
    Legend(fig[2,1:2],[MarkerElement(color=COLORS[i],marker=:circle) for i in 1:10],replace.(DATA,"_"=>" ");orientation=:horizontal,nbanks=2,framevisible=false,labelsize=15)
    Legend(fig[3,1:2],[LineElement(color=:black),LineElement(color=:gray,linestyle=:dash),MarkerElement(color=:black,marker=:circle),MarkerElement(color=:black,marker=:utriangle)],["Proportional fit (all ten)","Intercept fit (all ten)","Proportional holdout","Intercept holdout"];orientation=:horizontal,framevisible=false,labelsize=15)
    save(joinpath(DEST,"01_support_estimation.png"),fig;px_per_unit=2)
    fig=Figure(size=(1400,760))
    for (j,scenario) in enumerate(["Proportional","Intercept"])
        ax=Axis(fig[1,j];xlabel="Gap MAE (percentage points)",title=scenario=="Proportional" ? "Proportional support estimate" : "Regression with intercept",yticks=(1:10,replace.(DATA,"_"=>" ")))
        for (i,ds) in enumerate(DATA)
            z=only(eachrow(s[(s.dataset.==ds).&(s.scenario.==scenario).&(s.model.=="Pair"),:]))
            lines!(ax,[z.observed_M_mae_pp,z.mae_pp],[i,i];color=:gray,linewidth=1.5)
            scatter!(ax,[z.observed_M_mae_pp],[i];color="#CC6677",markersize=11);scatter!(ax,[z.mae_pp],[i];color="#4477AA",markersize=11)
        end
        ylims!(ax,.5,10.5);ax.yreversed=true;xlims!(ax,0,17)
    end
    Legend(fig[2,1:2],[MarkerElement(color="#CC6677",marker=:circle),MarkerElement(color="#4477AA",marker=:circle)],["Pair: observed M","Pair: held-out estimated M"];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"02_prediction_errors.png"),fig;px_per_unit=2)
    fig=Figure(size=(1800,880))
    for (i,ds) in enumerate(DATA)
        ax=Axis(fig[cld(i,5),mod1(i,5)];title=replace(ds,"_"=>" "),xlabel="Sites removed (%)",ylabel=mod1(i,5)==1 ? "Gap (pp)" : "")
        for (scenario,model,color,style) in [("Proportional","Pair","#4477AA",:solid),("Intercept","Pair","#228833",:dashdot),("Proportional","Common","#999999",:dot)]
            v=c[(c.dataset.==ds).&(c.scenario.==scenario).&(c.model.==model).&(c.removal_fraction.<=.8),:]
            lines!(ax,100v.removal_fraction,v.predicted_gap_pp;color=color,linestyle=style,linewidth=2.5)
            if scenario=="Proportional" && model=="Pair"
                lines!(ax,100v.removal_fraction,v.observed_M_gap_pp;color="#CC6677",linewidth=2.5)
                lines!(ax,100v.removal_fraction,v.empirical_gap_pp;color=:black,linestyle=:dash,linewidth=2.5)
            end
        end
        xlims!(ax,0,80);ylims!(ax,0,50)
    end
    Legend(fig[3,1:5],[LineElement(color=:black,linestyle=:dash),LineElement(color="#CC6677"),LineElement(color="#4477AA"),LineElement(color="#228833",linestyle=:dashdot),LineElement(color="#999999",linestyle=:dot)],["Empirical","Pair: observed M","Pair: proportional M estimate","Pair: intercept M estimate","Common: proportional M estimate"];orientation=:horizontal,nbanks=2,framevisible=false,labelsize=17)
    save(joinpath(DEST,"03_gap_curves.png"),fig;px_per_unit=2)
end
function report()
    a=CSV.read(joinpath(DEST,"across_dataset_summary.csv"),DataFrame);r=CSV.read(joinpath(DEST,"regression_folds.csv"),DataFrame);s=CSV.read(joinpath(DEST,"summary.csv"),DataFrame);full=CSV.read(joinpath(DEST,"full_data_regressions.csv"),DataFrame)
    curves=CSV.read(joinpath(DEST,"curves.csv"),DataFrame);freq=CSV.read(joinpath(DEST,"support_frequencies.csv"),DataFrame);verified=NamedTuple[]
    for z in eachrow(s)
        v=curves[(curves.dataset.==z.dataset).&(curves.scenario.==z.scenario).&(curves.model.==z.model),:]
        f=freq[(freq.dataset.==z.dataset).&(freq.scenario.==z.scenario).&(freq.model.==z.model),:]
        @assert v.m==collect(0:z.N) && all(v.removal_fraction.==v.m./z.N)
        @assert maximum(abs.(100gapmatrix(z.N)*f.predicted_frequency.-v.predicted_gap_pp))<1e-8
        @assert isapprox(sum(f.predicted_frequency),1;atol=1e-8) && isapprox(dot(f.K,f.predicted_frequency),z.predicted_M_used/z.L;atol=1e-8)
        keep=v.removal_fraction.<=.8;e=(v.predicted_gap_pp.-v.empirical_gap_pp)[keep]
        @assert isapprox(mean(abs.(e)),z.mae_pp;atol=1e-9) && isapprox(mean(e),z.signed_error_pp;atol=1e-9)
        @assert isapprox(100mean(abs.(e))/mean(v.empirical_gap_pp[keep]),z.relative_gap_error_percent;atol=1e-9)
        @assert maximum(abs.(v.predicted_gap_pp[[1,end]]))<1e-8
        push!(verified,(dataset=z.dataset,scenario=z.scenario,model=z.model,passed=true))
    end
    for z in eachrow(a)
        v=s[(s.scenario.==z.scenario).&(s.model.==z.model),:]
        @assert isapprox(mean(v.mae_pp),z.mean_mae_pp;atol=1e-10)
    end
    CSV.write(joinpath(DEST,"saved_result_checks.csv"),DataFrame(verified))
    fmt(x)=@sprintf("%.3f",x)
    open(joinpath(DEST,"FINDINGS.md"),"w") do io
        println(io,"# Predicting total support from all-pair co-occurrences\n\nThe earlier persistent-pair parameters were learned on nine datasets, but its target total support M was observed. Here both the M regression and the persistent-pair parameters exclude the target. Known regional interaction identities, N and occurrence records are still supplied.\n")
        println(io,"The co-occurrence proxy does not preserve the earlier accuracy. The primary proportional estimate increases equal-dataset Pair MAE from 0.976 to 5.739 percentage points; the intercept regression gives 8.410 pp. The proportional approach remains close for Montseny, Olot and Nahuel, but is substantially worse for Quercus, Gottin HP/PP, Salix Galpar and Garraf PP2. This identifies total-support estimation as a major source of error in this particular prediction pipeline, not a failure of co-occurrence information generally.\n\nObserved all-pair realisation fractions range from 8.3% to 34.2%. The regression overestimates support in several low-realisation systems and underestimates it in Salix Galpar and Garraf PP2. A positive across-dataset association therefore does not make one percentage sufficiently accurate for every target. Pair still improves mean gap error over Common when both receive the proportional M estimate (5.739 versus 7.146 pp), but neither matches the observed-M Pair reference.\n")
        for z in eachrow(full);println(io,"All-data $(z.model) regression: M = $(fmt(z.intercept)) + $(fmt(z.slope)) T; slope $(fmt(z.slope_percent))%. This is descriptive only and is not used for held-out curves.\n");end
        println(io,"T sums pair-site co-occurrences for all pairs, including K=0. The existing totalCOVstotalM_acrossDatasets.jl filters to K>0, so its slope is not the requested all-pair percentage. The primary model is proportional (zero intercept); the intercept regression is a prespecified sensitivity check, not selected according to validation performance.\n")
        println(io,"| Support estimate | Allocation | Equal-dataset mean MAE (pp) | Median (pp) | Signed error (pp) | Observed-M mean MAE (pp) | Datasets improved |\n|---|---|---:|---:|---:|---:|---:|")
        for z in eachrow(a);println(io,"| $(z.scenario) | $(z.model) | $(fmt(z.mean_mae_pp)) | $(fmt(z.median_mae_pp)) | $(fmt(z.mean_signed_error_pp)) | $(fmt(z.observed_M_mean_mae_pp)) | $(z.datasets_improved)/10 |");end
        println(io,"\n| Dataset | Observed M | Proportional M estimate, raw | Pair MAE: observed M | Pair MAE: proportional estimate | Pair MAE: intercept estimate |\n|---|---:|---:|---:|---:|---:|")
        for ds in DATA
            z=only(eachrow(r[(r.dataset.==ds).&(r.scenario.=="Proportional"),:]));v=only(eachrow(s[(s.dataset.==ds).&(s.scenario.=="Proportional").&(s.model.=="Pair"),:]));w=only(eachrow(s[(s.dataset.==ds).&(s.scenario.=="Intercept").&(s.model.=="Pair"),:]))
            println(io,"| $(replace(ds,"_"=>" ")) | $(z.observed_M) | $(fmt(z.predicted_M_raw)) | $(fmt(v.observed_M_mae_pp)) | $(fmt(v.mae_pp)) | $(fmt(w.mae_pp)) |")
        end
        println(io,"\nPrimary fold slopes range from $(fmt(minimum(r.slope_percent[r.scenario.=="Proportional"])))% to $(fmt(maximum(r.slope_percent[r.scenario.=="Proportional"])))%. No target M or individual K calibrates these slopes. OLS weights datasets equally in its squared-error objective, although high-T datasets have greater leverage. This is not the arithmetic mean of dataset percentages.\n")
        flagged=r[r.feasibility.!="Feasible",:]
        println(io,"Predictions are rounded to integer M, then projected onto the known-link feasible range L≤M≤sum(C). These constraints use known identities and occurrences, not held-out K. There are $(nrow(flagged)) projected dataset/scenario cases. Raw, rounded and used predictions are all retained.\n")
        for z in eachrow(flagged);println(io,"- $(z.dataset), $(z.scenario): raw $(fmt(z.predicted_M_raw)), used $(z.predicted_M_used); $(z.feasibility).");end
        println(io,"\nMean curves are exact conditional dynamic-programming expectations, with no allocation or removal simulation. Gap errors average all integer removal levels m/N≤0.8, including zero. Aggregate results weight each dataset equally. The local-retention fraction is still 1−m/N; the predicted M changes the support allocation and hence regional retention.\n\nThis is a plug-in prediction: regression and parameter uncertainty are not propagated. The primary comparison isolates replacing observed M with estimated M; Common is an allocation control. The site and beta models are not refitted or extended. Known links exclude wholly undetected interactions from the generated network, even though zero-K pairs contribute to the total-co-occurrence predictor.\n\nShared Gottin sites and related collections prevent interpreting leave-one-dataset-out as independent landscape validation. Much occurrence information is interaction-derived; unresolved occupancy corrections and the existing interaction flags remain. See README and copied metadata. Similar curves do not demonstrate an ecological mechanism.")
    end
    CSV.write(joinpath(DEST,"code_versions.csv"),DataFrame(file=["predict_M.jl"],sha256=[bytes2hex(sha256(read(@__FILE__)))]))
end
if !(length(ARGS)>=3 && ARGS[3]=="plot");main();end
plots();report()
end

module SupportDistributionFamilies
using CSV, DataFrames, Statistics, Random, LinearAlgebra, SHA, CairoMakie, Printf
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"../..")) : abspath(ARGS[1])
const SOURCE=joinpath(ROOT,"outputs/statistical description of local support")
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(SOURCE,"distribution_families")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const MODELS=["Exponential","Power_law","Cutoff_power_law","Discrete_lognormal","Singleton_geometric"]
const LABELS=["Exponential","Power law","Cutoff power law","Discrete lognormal","Singleton + geometric"]
const COLORS=["#888888","#CC8844","#4477AA","#228877","#882255"]
const B=499
const MIN_TAIL=50
const MIN_DISTINCT=5
mkpath(DEST)
function normprob(logw)
    mx=maximum(logw);w=exp.(logw.-mx);z=sum(w)
    w/z,mx+log(z)
end
function onefit(x,target)
    p,_=normprob(zeros(length(x)))
    target>=dot(p,x)-1e-13 && return 0.
    target<=minimum(x)+1e-13 && return 1000.
    lo=0.;hi=1.
    while dot(first(normprob(-hi.*x)),x)>target && hi<1e5;hi*=2;end
    for _ in 1:48
        mid=(lo+hi)/2
        if dot(first(normprob(-mid.*x)),x)>target;lo=mid;else;hi=mid;end
    end
    (lo+hi)/2
end
function objective(theta,X,target)
    p,z=normprob(-X*theta);mu=transpose(X)*p;g=target.-mu
    H=transpose(X)*(X.*p).-mu*transpose(mu)
    dot(theta,target)+z,g,H,p
end
function optimise(X,target,start,positive)
    theta=copy(start);best=Inf
    for iter in 1:150
        value,g,H,p=objective(theta,X,target)
        free=[j for j in eachindex(theta) if !positive[j] || theta[j]>1e-10 || g[j]<0]
        isempty(free) && break
        maximum(abs.(g[free]))<1e-9 && break
        delta=zeros(length(theta));delta[free]=-(H[free,free]+1e-12I)\g[free]
        dot(delta,g)>=0 && (delta[free]=-g[free])
        step=1.;accepted=false
        for _ in 1:45
            proposal=theta.+step.*delta
            for j in eachindex(theta);positive[j] && (proposal[j]=max(0.,proposal[j]));end
            nv=objective(proposal,X,target)[1]
            if nv<=value+1e-4dot(g,proposal.-theta)+1e-14
                theta=proposal;best=nv;accepted=true;break
            end
            step/=2
        end
        !accepted && break
    end
    val,g,_,p=objective(theta,X,target)
    residual=maximum([positive[j] && theta[j]<1e-8 ? max(0.,-g[j]) : abs(g[j]) for j in eachindex(theta)])
    (theta=theta,value=val,p=p,kkt=residual)
end
function fit(counts,model,xmin=1)
    N=length(counts);k=Float64.(xmin:N);n=sum(counts[xmin:end]);obs=counts[xmin:end]./n
    alpha=NaN;lambda=NaN;mu=NaN;sigma=NaN;s=NaN;r=NaN;beta=NaN;residual=0.;boundary=false
    if model=="Exponential"
        lambda=onefit(k,dot(obs,k));p,z=normprob(-lambda.*k);d=1;boundary=lambda==0
    elseif model=="Power_law"
        alpha=onefit(log.(k),dot(obs,log.(k)));p,z=normprob(-alpha.*log.(k));d=1;boundary=alpha==0
    elseif model in ["Cutoff_power_law","Discrete_lognormal"]
        lk=log.(k);alpha0=onefit(lk,dot(obs,lk))
        X=model=="Cutoff_power_law" ? hcat(lk,k) : hcat(lk,lk.^2)
        scale=[max(std(X[:,j]),1e-3) for j in 1:2];Xs=X./transpose(scale);target=transpose(Xs)*obs
        positive=model=="Cutoff_power_law" ? [true,true] : [false,true]
        starts=[scale.*[alpha0,0.],scale.*[alpha0,.1]]
        if model=="Cutoff_power_law";push!(starts,scale.*[0.,onefit(k,dot(obs,k))]);else;push!(starts,scale.*[-1.,1.]);end
        candidates=[optimise(Xs,target,start,positive) for start in starts];best=candidates[argmin([x.value for x in candidates])]
        theta=best.theta./scale;p=best.p;residual=best.kkt;alpha=theta[1];d=2
        if model=="Cutoff_power_law";lambda=theta[2];boundary=alpha<1e-8 || lambda<1e-8
        else;beta=theta[2];boundary=beta<1e-8;sigma=boundary ? Inf : sqrt(1/(2beta));mu=boundary ? NaN : (1-alpha)/(2beta);end
    else
        @assert xmin==1
        s=counts[1]/n;tail=counts[2:end];nt=sum(tail)
        lambda=nt==0 ? 1000. : onefit(Float64.(2:N),dot(tail./nt,Float64.(2:N)))
        pt,_=normprob(-lambda.*Float64.(2:N));p=vcat(s,(1-s).*pt);r=exp(-lambda);d=2;boundary=s==1 || lambda==0
    end
    @assert abs(sum(p)-1)<1e-10 && minimum(p)>=0
    ll=sum(counts[xmin+i-1]*log(p[i]) for i in eachindex(p) if counts[xmin+i-1]>0)
    D=maximum(abs.(cumsum(obs).-cumsum(p)));aic=2d-2ll;aicc=n>d+1 ? aic+2d*(d+1)/(n-d-1) : Inf
    (p=p,ll=ll,aic=aic,aicc=aicc,D=D,parameters=d,alpha=alpha,lambda=lambda,mu=mu,sigma=sigma,beta=beta,s=s,r=r,kkt=residual,boundary=boundary)
end
function samplecounts(rng,p,n)
    cdf=cumsum(p);cdf[end]=1.;counts=zeros(Int,length(p))
    for _ in 1:n;counts[searchsortedfirst(cdf,rand(rng))]+=1;end
    counts
end
function bootstrap(counts,model,f,rng)
    vals=zeros(B);maxkkt=0.
    for b in 1:B
        sim=samplecounts(rng,f.p,sum(counts));sf=fit(sim,model);vals[b]=sf.D;maxkkt=max(maxkkt,sf.kkt)
    end
    p=(1+count(>=(f.D-1e-12),vals))/(B+1)
    (p=p,mcse=sqrt(p*(1-p)/(B+1)),maxkkt=maxkkt)
end
function tailselect(counts)
    choices=NamedTuple[]
    for xmin in 2:length(counts)
        tail=counts[xmin:end]
        sum(tail)>=MIN_TAIL && count(tail.>0)>=MIN_DISTINCT || continue
        f=fit(counts,"Power_law",xmin);push!(choices,(xmin=xmin,n=sum(tail),f=f))
    end
    isempty(choices) ? nothing : choices[argmin([r.f.D for r in choices])]
end
function tailbootstrap(counts,selected,rng)
    n=sum(counts);xmin=selected.xmin;nt=selected.n;p=zeros(length(counts))
    p[1:xmin-1]=counts[1:xmin-1]./n;p[xmin:end]=(nt/n).*selected.f.p
    vals=Float64[]
    for b in 1:B
        sim=samplecounts(rng,p,n);f=tailselect(sim);isnothing(f) || push!(vals,f.f.D)
    end
    isempty(vals) && return (p=NaN,valid=0)
    (p=(1+count(>=(selected.f.D-1e-12),vals))/(length(vals)+1),valid=length(vals))
end
function main()
    links=CSV.read(joinpath(SOURCE,"link_support.csv"),DataFrame)
    sites=CSV.read(joinpath(SOURCE,"site_universe.csv"),DataFrame;types=Dict(:site=>String))
    versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
    for r in eachrow(versions);@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256;end
    used=DataFrame(path=String[],sha256=String[])
    for file in [joinpath(SOURCE,"link_support.csv"),joinpath(SOURCE,"site_universe.csv"),joinpath(SOURCE,"two_parameter_fit/fit_summary.csv"),joinpath(ROOT,"scripts/statistical_description_local_support.jl"),joinpath(ROOT,"scripts/two_parameter_support_fit.jl")]
        push!(used,(relpath(file,ROOT),bytes2hex(sha256(read(file)))))
    end
    CSV.write(joinpath(DEST,"input_versions.csv"),versions);CSV.write(joinpath(DEST,"source_versions.csv"),used)
    prior=CSV.read(joinpath(SOURCE,"two_parameter_fit/fit_summary.csv"),DataFrame)
    summaries=NamedTuple[];curves=NamedTuple[];tailrows=NamedTuple[];tailcurves=NamedTuple[];tailstatus=NamedTuple[];audit=NamedTuple[]
    for (di,ds) in enumerate(DATA)
        v=links[links.dataset.==ds,:];N=only(unique(v.N));K=Int.(v.K);L=length(K);counts=[count(==(k),K) for k in 1:N]
        @assert N==count(==(ds),sites.dataset) && all(1 .<=K .<=N)
        ints=CSV.read(joinpath(ROOT,"data/processed",ds*"_interactions.csv"),DataFrame;types=String)
        raw=Dict((g.consumer[1],g.resource[1])=>length(unique(g.site)) for g in groupby(ints,[:consumer,:resource]))
        @assert raw==Dict((r.consumer,r.resource)=>r.K for r in eachrow(v))
        for (mi,model) in enumerate(MODELS)
            f=fit(counts,model);@assert f.kkt<1e-6 "Fit convergence $ds $model $(f.kkt)"
            if model in ["Exponential","Singleton_geometric"]
                oldmodel=model=="Exponential" ? "one_parameter" : "two_parameter"
                old=only(eachrow(prior[(prior.dataset.==ds).&(prior.model.==oldmodel),:]))
                @assert abs(f.ll-old.log_likelihood)<1e-7
                push!(audit,(dataset=ds,check="Existing $model fit reproduced",error=abs(f.ll-old.log_likelihood),passed=true))
            end
            boot=bootstrap(counts,model,f,MersenneTwister(51000+100di+mi))
            push!(summaries,(dataset=ds,model=model,N=N,L=L,maximum_K=maximum(K),alpha=f.alpha,lambda=f.lambda,mu=f.mu,sigma=f.sigma,beta=f.beta,s=f.s,r=f.r,parameters=f.parameters,log_likelihood=f.ll,AIC=f.aic,AICc=f.aicc,KS=f.D,bootstrap_p=boot.p,bootstrap_p_mcse=boot.mcse,bootstrap_replicates=B,boundary=f.boundary,kkt_residual=f.kkt,bootstrap_max_kkt=boot.maxkkt,observed_singleton=counts[1]/L,fitted_singleton=f.p[1],observed_mean=mean(K),fitted_mean=dot(1:N,f.p),observed_p90=sort(K)[ceil(Int,.9L)],fitted_p90=findfirst(>=(.9),cumsum(f.p))))
            obs=counts/L;occ=reverse(cumsum(reverse(obs)));fcc=reverse(cumsum(reverse(f.p)))
            for k in 1:N;push!(curves,(dataset=ds,model=model,K=k,observed_frequency=obs[k],fitted_frequency=f.p[k],observed_ccdf=occ[k],fitted_ccdf=fcc[k]));end
            println(ds," ",model,": AICc=",round(f.aicc;digits=2)," KS p=",boot.p);flush(stdout)
        end
        ts=tailselect(counts)
        if isnothing(ts)
            push!(tailstatus,(dataset=ds,status="Insufficient tail: need 50 links and 5 distinct K at xmin>=2",xmin=0,tail_links=0,tail_fraction=0.,observed_support_ratio=NaN,bootstrap_p=NaN,valid_bootstrap=0))
        else
            tb=tailbootstrap(counts,ts,MersenneTwister(83000+di));xmin=ts.xmin
            push!(tailstatus,(dataset=ds,status="Exploratory finite-support tail fit",xmin=xmin,tail_links=ts.n,tail_fraction=ts.n/L,observed_support_ratio=maximum(K)/xmin,bootstrap_p=tb.p,valid_bootstrap=tb.valid))
            for model in MODELS[1:4]
                f=fit(counts,model,xmin);@assert f.kkt<1e-6
                push!(tailrows,(dataset=ds,model=model,N=N,xmin=xmin,tail_links=ts.n,alpha=f.alpha,lambda=f.lambda,mu=f.mu,sigma=f.sigma,log_likelihood=f.ll,AICc=f.aicc,KS=f.D,boundary=f.boundary,kkt_residual=f.kkt))
                occ=reverse(cumsum(reverse(counts[xmin:end]./ts.n)));fcc=reverse(cumsum(reverse(f.p)))
                for (i,k) in enumerate(xmin:N);push!(tailcurves,(dataset=ds,model=model,K=k,observed_ccdf=occ[i],fitted_ccdf=fcc[i]));end
            end
        end
        CSV.write(joinpath(DEST,"fit_summary.csv"),DataFrame(summaries));CSV.write(joinpath(DEST,"fitted_distributions.csv"),DataFrame(curves))
        CSV.write(joinpath(DEST,"tail_selection.csv"),DataFrame(tailstatus));CSV.write(joinpath(DEST,"tail_fits.csv"),DataFrame(tailrows));CSV.write(joinpath(DEST,"tail_curves.csv"),DataFrame(tailcurves))
    end
    tab=DataFrame(summaries);tab.delta_AICc=[r.AICc-minimum(tab.AICc[tab.dataset.==r.dataset]) for r in eachrow(tab)]
    CSV.write(joinpath(DEST,"fit_summary.csv"),tab)
    tail=DataFrame(tailrows);tail.delta_AICc=[r.AICc-minimum(tail.AICc[tail.dataset.==r.dataset]) for r in eachrow(tail)]
    CSV.write(joinpath(DEST,"tail_fits.csv"),tail);CSV.write(joinpath(DEST,"fit_checks.csv"),DataFrame(audit))
end
function plots()
    t=CSV.read(joinpath(DEST,"fit_summary.csv"),DataFrame);d=CSV.read(joinpath(DEST,"fitted_distributions.csv"),DataFrame)
    set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
    for (kind,file) in [("frequency","01_frequency_fits"),("semilog","02_ccdf_semilog"),("loglog","03_ccdf_loglog")]
        f=Figure(size=(1830,890))
        for (i,ds) in enumerate(DATA)
            row=only(eachrow(t[(t.dataset.==ds).&(t.model.=="Exponential"),:]));maxk=row.maximum_K
            ax=Axis(f[cld(i,5),mod1(i,5)];title=replace(ds,"_"=>" "),xlabel="Supporting sites, K",ylabel=mod1(i,5)==1 ? (kind=="frequency" ? "Proportion of regional links" : "Proportion with support ≥ K") : "",xscale=kind=="loglog" ? log10 : identity,yscale=kind=="frequency" ? identity : log10)
            for (mi,model) in enumerate(MODELS)
                v=d[(d.dataset.==ds).&(d.model.==model),:];field=kind=="frequency" ? :fitted_frequency : :fitted_ccdf;mask=v[!,field].>0
                lines!(ax,v.K[mask],v[mask,field];color=COLORS[mi],linewidth=2,linestyle=model=="Singleton_geometric" ? :dash : :solid)
                if mi==1
                    field=kind=="frequency" ? :observed_frequency : :observed_ccdf;mask=v[!,field].>0
                    if kind=="frequency";scatter!(ax,v.K[mask],v[mask,field];color=:black,markersize=6)
                    else;stairs!(ax,v.K[mask],v[mask,field];color=:black,linewidth=1.5,step=:post);scatter!(ax,v.K[mask],v[mask,field];color=:black,markersize=4);end
                end
            end
            xlims!(ax,kind=="loglog" ? .95 : .5,min(row.N,maxk+3)+.5)
            ylims!(ax,kind=="frequency" ? 0. : .5/row.L,kind=="frequency" ? .85 : 1.1)
            if kind=="loglog";ticks=filter(x->x<=maxk+3,[1,2,5,10,20,50,100]);ax.xticks=(ticks,string.(ticks));end
        end
        Legend(f[3,1:5],vcat([LineElement(color=COLORS[i],linestyle=i==5 ? :dash : :solid) for i in 1:5],[MarkerElement(color=:black,marker=:circle)]),vcat(LABELS,["Empirical"]);orientation=:horizontal,framevisible=false,labelsize=16)
        save(joinpath(DEST,file*".png"),f;px_per_unit=2)
    end
    f=Figure(size=(1500,790))
    for (j,field,title) in [(1,:delta_AICc,"ΔAICc (lower is better)"),(2,:bootstrap_p,"Bootstrap goodness-of-fit p")]
        ax=Axis(f[1,j];title=title,xticks=(1:5,["Exp.","Power","Cutoff","Lognormal","Singleton"]),yticks=(1:10,replace.(DATA,"_"=>" ")),xticklabelrotation=pi/6)
        for (i,ds) in enumerate(DATA),(mi,m) in enumerate(MODELS)
            r=only(eachrow(t[(t.dataset.==ds).&(t.model.==m),:]));v=r[field]
            shade=j==1 ? clamp(v/30,0,1) : 1-clamp(v/.5,0,1)
            poly!(ax,Rect2f(mi-.48,i-.46,.96,.92);color=RGBAf(.94-.35shade,.97-.32shade,.98-.25shade,1))
            label=j==1 ? @sprintf("%.1f",v) : @sprintf("%.3f",v)
            text!(ax,mi,i;text=label,align=(:center,:center),fontsize=17,color=:black)
        end
        xlims!(ax,.5,5.5);ylims!(ax,.5,10.5);ax.yreversed=true;hidespines!(ax)
    end
    save(joinpath(DEST,"04_model_comparison.png"),f;px_per_unit=2)
    ts=CSV.read(joinpath(DEST,"tail_selection.csv"),DataFrame);tc=CSV.read(joinpath(DEST,"tail_curves.csv"),DataFrame)
    f=Figure(size=(1830,890))
    for (i,ds) in enumerate(DATA)
        r=only(eachrow(ts[ts.dataset.==ds,:]))
        ax=Axis(f[cld(i,5),mod1(i,5)];title=replace(ds,"_"=>" "),xlabel="Supporting sites, K",ylabel=mod1(i,5)==1 ? "Conditional tail CCDF" : "",xscale=log10,yscale=log10)
        if r.xmin==0
            text!(ax,2,.5;text="Insufficient tail data",fontsize=15);xlims!(ax,1,10);ylims!(ax,.01,1.1);continue
        end
        for (mi,m) in enumerate(MODELS[1:4])
            v=tc[(tc.dataset.==ds).&(tc.model.==m),:];mask=v.fitted_ccdf.>0
            lines!(ax,v.K[mask],v.fitted_ccdf[mask];color=COLORS[mi],linewidth=2)
            if mi==1;keep=v.observed_ccdf.>0;stairs!(ax,v.K[keep],v.observed_ccdf[keep];color=:black,linewidth=1.5,step=:post);end
        end
        maxk=only(t.maximum_K[(t.dataset.==ds).&(t.model.=="Exponential")]);xlims!(ax,.95r.xmin,maxk+2);ylims!(ax,.5/r.tail_links,1.1)
        ticks=sort(unique(vcat([r.xmin],filter(x->r.xmin<=x<=maxk+2,[2,5,10,20,50,100]))));ax.xticks=(ticks,string.(ticks))
        text!(ax,r.xmin,.8/r.tail_links;text="Kmin = $(r.xmin); links = $(r.tail_links)",fontsize=13,align=(:left,:bottom))
    end
    Legend(f[3,1:5],vcat([LineElement(color=COLORS[i]) for i in 1:4],[LineElement(color=:black)]),vcat(LABELS[1:4],["Empirical tail"]);orientation=:horizontal,framevisible=false,labelsize=16)
    save(joinpath(DEST,"05_tail_comparison.png"),f;px_per_unit=2)
end
function report()
    t=CSV.read(joinpath(DEST,"fit_summary.csv"),DataFrame);ts=CSV.read(joinpath(DEST,"tail_selection.csv"),DataFrame)
    fmt(x)=@sprintf("%.3f",x)
    open(joinpath(DEST,"FINDINGS.md"),"w") do io
        println(io,"# What distribution describes local support?\n")
        println(io,"K is the number of distinct sites supporting a known regional link (K≥1). These are fits to the empirical support distribution, not degree distributions, co-occurrence models, or gap predictions. All candidate laws are discrete and normalised on 1,...,N, including unobserved support levels.\n")
        println(io,"A single exponential shows detectable mismatch in all ten datasets. A bounded power law is compatible with eight, but it is not uniquely supported: several alternatives have similar AICc. Salix Galpar clearly favours an exponential cutoff, Quercus favours the singleton-plus-geometric model, and no tested family adequately describes Garraf PP2. These results support recurrent broad distributions with dataset-specific details, rather than one established universal law.\n")
        println(io,"| Dataset | Lowest-AICc family | KS bootstrap p for winner | Other families within 2 AICc units |\n|---|---|---:|---|")
        for ds in DATA
            v=t[t.dataset.==ds,:];r=v[argmin(v.AICc),:];other=v.model[(v.delta_AICc.<=2).&(v.model.!=r.model)]
            println(io,"| $(replace(ds,"_"=>" ")) | $(replace(r.model,"_"=>" ")) | $(fmt(r.bootstrap_p)) | $(isempty(other) ? "—" : join(replace.(other,"_"=>" "),", ")) |")
        end
        println(io,"\nAICc ranks candidates relatively. A low bootstrap p flags systematic distributional mismatch even if a model has the lowest AICc. Conversely, p≥0.05 means this check did not detect mismatch; it does not establish a true family or universality. Values near 0.05 should be read with their Monte Carlo errors.\n")
        println(io,"| Family | Datasets with bootstrap p ≥ 0.05 | Lowest-AICc datasets |\n|---|---:|---:|")
        for m in MODELS
            wins=count(ds->t.model[findall(t.dataset.==ds)][argmin(t.AICc[t.dataset.==ds])]==m,DATA)
            println(io,"| $(replace(m,"_"=>" ")) | $(count(t.bootstrap_p[t.model.==m].>=.05))/10 | $wins/10 |")
        end
        println(io,"\n## Tail-only check\n\nThe lower cutoff is selected by minimum power-law KS distance among cutoffs Kmin≥2 retaining at least 50 links and five distinct observed support values. The power-law bootstrap regenerates the tail, resamples the empirical body, and reselects the cutoff. Alternatives are compared on that same selected tail. Tail AICc values must not be compared with whole-distribution AICc.\n")
        println(io,"| Dataset | Kmin | Tail links | Observed Kmax/Kmin | Power-law bootstrap p |\n|---|---:|---:|---:|---:|")
        for r in eachrow(ts)
            println(io,r.xmin==0 ? "| $(replace(r.dataset,"_"=>" ")) | Insufficient data | — | — | — |" : "| $(replace(r.dataset,"_"=>" ")) | $(r.xmin) | $(r.tail_links) | $(round(r.observed_support_ratio;digits=2)) | $(fmt(r.bootstrap_p)) |")
        end
        println(io,"\nAny compatible tail is a limited finite-range description. K is bounded by N and observed tails cover a restricted range; these data cannot establish an unbounded, scale-free asymptotic law. A cutoff power law here adds an exponential decay below the hard upper bound N; the hard N bound applies to every model.\n")
        println(io,"Tail bootstrap samples without an eligible tail are excluded from the conditional bootstrap calculation. Gottin HP retains 420 of 499 samples and Quercus 498; other eligible datasets retain all 499. This selection particularly limits interpretation for Gottin HP. Valid replicate counts are recorded in tail_selection.csv.\n")
        println(io,"## Interpretation limits\n\nThe tests treat regional-link support counts as independent draws. Links share species and sites, so bootstrap p-values and AICc are working-model diagnostics, not independence-robust ecological tests. Neither counts of non-rejections nor common rank winners establish universality. Site counts, sampling effort, spatial extent and observation can affect shape; K/N alone would not eliminate those effects. Shared campaigns further reduce independent replication.\n")
        println(io,"Existing geometric and singleton-geometric likelihoods are reproduced. Input interaction triples and site counts are checked against the descriptive analysis. Occurrence restorations do not alter these K distributions because neither interactions nor N changed; outstanding occurrence and interaction flags remain in README. No original analysis is overwritten.\n")
        println(io,"Methodological references: Clauset, Shalizi & Newman (2009), https://doi.org/10.1137/070710111; Broido & Clauset (2019), https://doi.org/10.1038/s41467-019-08746-5. Our finite-support, discrete bootstrap is implemented for this bounded site-support problem, not copied as an unbounded power-law test.")
    end
    CSV.write(joinpath(DEST,"code_versions.csv"),DataFrame(file=["fit_distributions.jl"],sha256=[bytes2hex(sha256(read(@__FILE__)))]))
end
if !(length(ARGS)>=3 && ARGS[3] in ["plot","report"]);main();end
if !(length(ARGS)>=3 && ARGS[3]=="report");plots();end
report()
end

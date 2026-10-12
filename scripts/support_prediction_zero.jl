module SupportPredictionZero
using CSV, DataFrames, Statistics, LinearAlgebra, Printf, SHA, CairoMakie
include(joinpath(@__DIR__, "..", "src", "config.jl"))
include(joinpath(ROOT, "src", "site_removal.jl"))

const SCENARIO = "--inclusive" in ARGS ? "inclusive" : "confirmed"
const DEST = joinpath(OUT, "support_prediction_zero", SCENARIO)
const ZERO_INPUT = joinpath(DEST,"inputs")
const H = collect(0.0:0.02:0.8)
const COLORS = ["#D55E00", "#0072B2"]
const MODELS = ["M0", "M1"]
const BLOCK = Dict(d => (startswith(d,"Gottin") ? "Gottin" :
    (startswith(d,"Garraf") || d in ["Olot","Montseny"]) ? "Bosch_collection" : d) for d in DATA)
const LOGFACT = cumsum(vcat(0.0, log.(Float64.(1:2000))))
const LOWER = [-12.0, -10.0]
const UPPER = [12.0, log(1e8)]
logchoose(n,k) = LOGFACT[n+1]-LOGFACT[k+1]-LOGFACT[n-k+1]
logistic(x) = 1/(1+exp(-x))
logit(x) = log(x/(1-x))
area(y) = sum(diff(H).*(y[1:end-1].+y[2:end])./2)/last(H)
savecsv(name,x) = CSV.write(joinpath(DEST,name*".csv"),x isa DataFrame ? x : DataFrame(x))

# Rising products avoid subtracting nearly equal log-gamma terms at large phi.
function logprob(n,k,mu,phi)
    v = logchoose(n,k)+k*log(mu)+(n-k)*log1p(-mu)
    if isfinite(phi)
        v += sum(log1p(i/(mu*phi)) for i in 0:k-1;init=0.0)
        v += sum(log1p(i/((1-mu)*phi)) for i in 0:n-k-1;init=0.0)
        v -= sum(log1p(i/phi) for i in 0:n-1;init=0.0)
    end
    v
end

function read_pairs()
    allpairs, versions, audits = DataFrame(), NamedTuple[], NamedTuple[]
    expected = Dict(r.file=>r.sha256 for r in eachrow(CSV.read(joinpath(ROOT,"data/checksums.csv"),DataFrame)))
    for ds in DATA
        inputs = Dict{String,DataFrame}()
        for kind in ["interactions","cooccurrences","occupancy"]
            name=ds*"_"*kind*".csv"; path=joinpath(ZERO_INPUT,name)
            digest=bytes2hex(sha256(read(path)))
            rawdigest=bytes2hex(sha256(read(joinpath(RAW,name))))
            @assert rawdigest==expected[name] "Raw inputs changed; review occurrence provenance before rerunning"
            df=CSV.read(path,DataFrame;types=String)
            @assert nrow(unique(df))==nrow(df) && !any(ismissing,Matrix(df))
            inputs[kind]=df
            push!(versions,(dataset=ds,kind=kind,sha256=digest,raw_sha256=rawdigest,rows=nrow(df)))
        end
        ints,co,occ=inputs["interactions"],inputs["cooccurrences"],inputs["occupancy"]
        @assert Set(Tuple.(eachrow(ints))) ⊆ Set(Tuple.(eachrow(co)))
        N=length(unique(vcat(ints.site,co.site,occ.site)))
        n=combine(groupby(co,[:consumer,:resource]),nrow=>:n)
        k=combine(groupby(ints,[:consumer,:resource]),nrow=>:K)
        p=leftjoin(n,k;on=[:consumer,:resource]);p.K=Int.(coalesce.(p.K,0))
        p.dataset=fill(ds,nrow(p));p.N=fill(N,nrow(p));p.realised=p.K.>0
        @assert all(0 .<= p.K .<= p.n .<= N)
        append!(allpairs,select(p,:dataset,:consumer,:resource,:N,:n,:K,:realised))
        linked=Set(vcat([(r.site,"consumer",r.consumer) for r in eachrow(ints)],
                       [(r.site,"resource",r.resource) for r in eachrow(ints)]))
        extra(g)=count(r->r.trophic_level==g && !((r.site,g,r.species) in linked),eachrow(occ))
        push!(audits,(dataset=ds,N=N,pairs=nrow(p),links=count(p.realised),total_support=sum(p.K),
            p_emp=sum(p.K)/sum(p.n),r_regional=mean(p.realised),
            consumer_presence=SCENARIO=="inclusive" && ds=="Gottin_PP" ? "includes uncertain missing-plant records" : "baseline occupancy",resource_presence="baseline plus documented additions; see presence_record_provenance.csv",
            unlinked_consumer_sites=extra("consumer"),unlinked_resource_sites=extra("resource"),
            occurrence_version="zero-presence sensitivity: "*SCENARIO,holdout_block=BLOCK[ds]))
    end
    savecsv("pairs",allpairs);savecsv("input_versions",versions);savecsv("occurrence_audit",audits)
    allpairs
end

function objective_rows(p,train)
    rows=NamedTuple[]
    for ds in train
        d=p[p.dataset.==ds,:]
        for r in eachrow(combine(groupby(d,[:n,:K]),nrow=>:count))
            push!(rows,(n=r.n,k=r.K,weight=r.count/(nrow(d)*length(train))))
        end
    end
    rows
end

function simplex(f,start;maxiter=1500)
    x=hcat(clamp.(start,LOWER,UPPER),clamp.(start+[.35,0],LOWER,UPPER),clamp.(start+[0,.7],LOWER,UPPER))
    values=[f(x[:,i]) for i in 1:3]
    for iter in 1:maxiter
        order=sortperm(values);x=x[:,order];values=values[order]
        if maximum(abs.(x[:,2:3].-x[:,1]))<1e-6 && maximum(abs.(values.-values[1]))<1e-9
            return (x=x[:,1],value=values[1],converged=true,iterations=iter)
        end
        center=(x[:,1]+x[:,2])/2
        xr=clamp.(2center-x[:,3],LOWER,UPPER);fr=f(xr)
        if fr<values[1]
            xe=clamp.(center+2(xr-center),LOWER,UPPER);fe=f(xe)
            x[:,3],values[3]=fe<fr ? (xe,fe) : (xr,fr)
        elseif fr<values[2]
            x[:,3]=xr;values[3]=fr
        else
            xc=clamp.(center+.5*(fr<values[3] ? xr-center : x[:,3]-center),LOWER,UPPER);fc=f(xc)
            if fc<min(fr,values[3])
                x[:,3]=xc;values[3]=fc
            else
                for i in 2:3;x[:,i]=(x[:,1]+x[:,i])/2;values[i]=f(x[:,i]);end
            end
        end
    end
    i=argmin(values)
    (x=x[:,i],value=values[i],converged=false,iterations=maxiter)
end

function fitmodels(p,train,fitid,starts)
    rows=objective_rows(p,train)
    mu0=sum(r.weight*r.k for r in rows)/sum(r.weight*r.n for r in rows)
    nll(mu,phi)=-sum(r.weight*logprob(r.n,r.k,mu,phi) for r in rows)
    best=nothing
    for (i,(mu,phi)) in enumerate(vcat([(mu0,v) for v in [.03,.3,3.,30.,300.,1e5]],[(.05,.3),(.5,.3)]))
        fit=simplex(x->nll(logistic(x[1]),exp(x[2])),[logit(mu),log(phi)])
        push!(starts,(fit_id=fitid,start=i,mu=logistic(fit.x[1]),phi=exp(fit.x[2]),nll=fit.value,
                      converged=fit.converged,iterations=fit.iterations))
        if isnothing(best) || fit.value<best.value;best=fit;end
    end
    @assert best.value<=nll(mu0,Inf)+1e-6 "M1 failed to recover binomial limit"
    mu1,phi1=logistic(best.x[1]),exp(best.x[2])
    [ (model="M0",mu=mu0,phi=Inf,nll=nll(mu0,Inf),converged=true,iterations=0,boundary=false),
      (model="M1",mu=mu1,phi=phi1,nll=best.value,converged=best.converged,iterations=best.iterations,
       boundary=any(abs.(best.x.-LOWER).<1e-4)||any(abs.(best.x.-UPPER).<1e-4)) ]
end

function survival(k,N,h)
    m=N*h;a=floor(Int,m);b=ceil(Int,m);w=m-a
    1-((1-w)*loss(k,N,a)+w*loss(k,N,b))
end

function evaluate(d,fit)
    N=first(d.N);hist=zeros(N+1);obs=zeros(N+1)
    for k in d.K;obs[k+1]+=1;end
    for r in eachrow(combine(groupby(d,:n),nrow=>:count))
        probs=[exp(logprob(r.n,k,fit.mu,fit.phi)) for k in 0:r.n]
        @assert isapprox(sum(probs),1;atol=1e-9)
        hist[1:r.n+1].+=r.count.*probs
    end
    L=sum(hist[2:end]);T=sum((0:N).*hist);Lo=sum(obs[2:end]);To=sum(d.K)
    @assert isapprox(T,fit.mu*sum(d.n);rtol=1e-9)
    R=[sum(hist[k+1]*survival(k,N,h) for k in 1:N)/L for h in H]
    Ro=[sum(obs[k+1]*survival(k,N,h) for k in 1:N)/Lo for h in H]
    G=R.-(1 .-H);Go=Ro.-(1 .-H)
    @assert isapprox(R[1],1) && all(diff(R).<=1e-10) && all(G.>=-1e-10)
    scores=(nlpd=-mean(logprob(r.n,r.K,fit.mu,fit.phi) for r in eachrow(d)),
        observed_total=To,predicted_total=T,total_relative_error=T/To-1,
        observed_links=Lo,predicted_links=L,links_relative_error=L/Lo-1,
        observed_singleton=obs[2]/Lo,predicted_singleton=hist[2]/L,
        observed_mean_support=To/Lo,predicted_mean_support=T/L,
        retention_mae=mean(abs.(R.-Ro)),observed_gap=area(Go),predicted_gap=area(G),gap_error=area(G)-area(Go))
    scores,hist,obs,R,Ro
end

function analyse(p)
    parameters,starts,metrics,hists,curves,conditional=NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[],NamedTuple[]
    jobs=Tuple{String,String,Vector{String},Vector{String}}[]
    for ds in DATA
        push!(jobs,("LODO",ds,filter(!=(ds),DATA),[ds]))
        push!(jobs,("Within",ds,[ds],[ds]))
    end
    for block in unique([BLOCK[d] for d in DATA])
        push!(jobs,("Grouped",block,filter(d->BLOCK[d]!=block,DATA),filter(d->BLOCK[d]==block,DATA)))
    end
    for (scheme,label,train,targets) in jobs
        id=scheme*":"*label
        fits=fitmodels(p,train,id,starts)
        for fit in fits
            push!(parameters,(scheme=scheme,fit_id=id,training_datasets=join(train,";"),test_datasets=join(targets,";"),fit...))
            for ds in targets
                d=p[p.dataset.==ds,:];N=first(d.N)
                scores,hist,obs,R,Ro=evaluate(d,fit)
                push!(metrics,(scheme=scheme,dataset=ds,model=fit.model,fit_id=id,scores...))
                if scheme=="Within" && fit.model=="M0"
                    @assert isapprox(scores.predicted_total,scores.observed_total;rtol=1e-12)
                end
                for k in 0:N
                    push!(hists,(scheme=scheme,dataset=ds,model=fit.model,k=k,observed_count=obs[k+1],predicted_count=hist[k+1],
                        observed_realised_frequency=k==0 ? NaN : obs[k+1]/scores.observed_links,
                        predicted_realised_frequency=k==0 ? NaN : hist[k+1]/scores.predicted_links))
                end
                for (i,h) in enumerate(H)
                    push!(curves,(scheme=scheme,dataset=ds,model=fit.model,removal=h,local_retention=1-h,
                        observed_regional=Ro[i],predicted_regional=R[i],observed_gap=Ro[i]-(1-h),predicted_gap=R[i]-(1-h)))
                end
                for n in sort(unique(d.n))
                    v=d[d.n.==n,:];real=v.K[v.K.>0]
                    realised_probability=-expm1(logprob(n,0,fit.mu,fit.phi))
                    push!(conditional,(scheme=scheme,dataset=ds,model=fit.model,n=n,pairs=nrow(v),realised_pairs=length(real),
                        observed_realisation=isempty(real) ? NaN : mean(real)/n,
                        predicted_realisation=fit.mu/realised_probability,
                        observed_regional_fraction=length(real)/nrow(v),predicted_regional_fraction=realised_probability))
                end
            end
        end
        println(id," M1 mu=",round(fits[2].mu;digits=4)," phi=",round(fits[2].phi;digits=4)," converged=",fits[2].converged)
    end
    for (name,rows) in [("parameters",parameters),("optimisation_starts",starts),("metrics",metrics),
                        ("support_histograms",hists),("retention_curves",curves),("conditional_realisation",conditional)]
        savecsv(name,rows)
    end
    m=DataFrame(metrics)
    summary=combine(groupby(m,[:scheme,:model]),:nlpd=>mean=>:mean_nlpd,
        :total_relative_error=>(x->mean(abs.(x)))=>:mean_absolute_total_error,
        :links_relative_error=>(x->mean(abs.(x)))=>:mean_absolute_links_error,
        :retention_mae=>mean=>:mean_retention_mae,:gap_error=>(x->mean(abs.(x)))=>:mean_absolute_gap_error)
    savecsv("equal_dataset_summary",summary)
end

function checks()
    @assert isapprox(exp(logprob(4,2,.3,Inf)),6*.3^2*.7^2;atol=1e-12)
    @assert all(isapprox.(exp.([logprob(4,k,.5,2.) for k in 0:4]),.2;atol=1e-12))
    @assert isapprox(logprob(30,7,.3,1e8),logprob(30,7,.3,Inf);atol=1e-5)
    @assert isapprox(survival(2,10,.5),1-2/9;atol=1e-12)
    @assert isapprox(survival(1,8,.8),.2;atol=1e-12)
    @assert isapprox(.3/(-expm1(logprob(1,0,.3,Inf))),1;atol=1e-12)
end

function exportfig(f,name)
    save(joinpath(DEST,name*".png"),f;px_per_unit=2)
    save(joinpath(DEST,name*".pdf"),f)
end
function legend(slot;local_line=false)
    elements=[LineElement(color=:black,linewidth=2),LineElement(color=COLORS[1],linewidth=2),LineElement(color=COLORS[2],linewidth=2)]
    labels=["Observed","M0: constant probability","M1: heterogeneous probabilities"]
    if local_line;push!(elements,LineElement(color=:gray60,linestyle=:dash));push!(labels,"Local support");end
    Legend(slot,elements,labels;orientation=:horizontal,framevisible=false)
end

function figures()
    readtable(n)=CSV.read(joinpath(DEST,n*".csv"),DataFrame)
    m=readtable("metrics");hs=readtable("support_histograms");cs=readtable("retention_curves");qs=readtable("conditional_realisation")
    set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,
        xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16),Legend=(labelsize=17,)))
    f=Figure(size=(1300,1100))
    for (j,ob,pr,label,islog) in [(1,:observed_total,:predicted_total,"Total local support",true),
        (2,:observed_links,:predicted_links,"Regional link richness",true),
        (3,:observed_gap,:predicted_gap,"Mean retention gap",false),
        (4,:observed_singleton,:predicted_singleton,"Fraction of realised links with K = 1",false)]
        ax=Axis(f[cld(j,2),mod1(j,2)],xlabel="Observed "*lowercase(label),ylabel="Predicted "*lowercase(label),
            xscale=islog ? log10 : identity,yscale=islog ? log10 : identity)
        d=m[m.scheme.=="LODO",:];lo=min(minimum(d[!,ob]),minimum(d[!,pr]));hi=max(maximum(d[!,ob]),maximum(d[!,pr]))
        limits=islog ? (lo*.75,hi*1.45) : (0.,hi*1.15)
        lines!(ax,collect(limits),collect(limits),color=:gray65,linestyle=:dash)
        for (model,color) in zip(MODELS,COLORS)
            v=d[d.model.==model,:]
            scatter!(ax,v[!,ob],v[!,pr],color=color,markersize=12)
        end
        xlims!(ax,limits...);ylims!(ax,limits...)
    end
    Legend(f[3,1:2],[MarkerElement(color=c,marker=:circle) for c in COLORS],MODELS,orientation=:horizontal,framevisible=false)
    exportfig(f,"01_transfer_accuracy")
    for (kind,data,name) in [("support",hs,"02_support_distributions"),("retention",cs,"03_retention_curves"),
                             ("conditional",qs,"06_conditional_realisation"),("regional",qs,"07_regional_realisation")]
        f=Figure(size=(1900,850))
        for (j,ds) in enumerate(DATA)
            v=data[(data.dataset.==ds).&(data.scheme.=="LODO"),:]
            ax=Axis(f[cld(j,5),mod1(j,5)],title=DISPLAY[ds],titlesize=18,
                xlabel=kind=="retention" ? "Sites removed (%)" : kind=="support" ? "Supporting sites, K" : "Co-occurring sites, n",
                ylabel=mod1(j,5)==1 ? (kind=="support" ? "Fraction with support ≥ K" : kind=="retention" ? "Regional links retained" :
                    kind=="conditional" ? "Mean K/n among realised links" : "Fraction of pairs realised") : "",
                xscale=kind=="retention" ? identity : log10)
            for (midx,model) in enumerate(MODELS)
                d=v[v.model.==model,:];color=COLORS[midx]
                if kind=="support"
                    d=sort(d[d.k.>0,:],:k);yp=reverse(cumsum(reverse(d.predicted_realised_frequency)))
                    yo=reverse(cumsum(reverse(d.observed_realised_frequency)))
                    lines!(ax,d.k,yp,color=color,linewidth=2.5)
                    midx==1 && lines!(ax,d.k,yo,color=:black,linewidth=2)
                elseif kind=="retention"
                    lines!(ax,100d.removal,d.predicted_regional,color=color,linewidth=2.5)
                    if midx==1
                        lines!(ax,100d.removal,d.observed_regional,color=:black,linewidth=2)
                        lines!(ax,100d.removal,d.local_retention,color=:gray60,linestyle=:dash)
                    end
                else
                    y=kind=="conditional" ? :predicted_realisation : :predicted_regional_fraction
                    o=kind=="conditional" ? :observed_realisation : :observed_regional_fraction
                    lines!(ax,d.n,d[!,y],color=color,linewidth=2.5)
                    if midx==1
                        keep=isfinite.(d[!,o])
                        counts=kind=="conditional" ? d.realised_pairs[keep] : d.pairs[keep]
                        scatter!(ax,d.n[keep],d[keep,o],color=(:black,.65),markersize=6 .+ 2log10.(counts))
                    end
                end
            end
            ylims!(ax,-.02,1.04)
            if kind!="retention"
                xmax=kind=="support" ? maximum(v.k) : maximum(v.n)
                ticks=filter(x->x<=xmax,[1,2,5,10,20,50,100,200])
                ax.xticks=(ticks,string.(ticks))
            end
        end
        legend(f[3,1:5];local_line=kind=="retention");colgap!(f.layout,25)
        exportfig(f,name)
    end
    f=Figure(size=(1450,650))
    for (j,model) in enumerate(MODELS)
        ax=Axis(f[1,j],xlabel="Dataset",ylabel="Mean negative log predictive probability",
            xticks=(1:10,[DISPLAY[x] for x in DATA]),xticklabelrotation=pi/4,title=model,titlesize=20)
        ylims!(ax,0,1.08maximum(m.nlpd))
        for (scheme,color,marker,offset) in [("LODO",COLORS[j],:circle,-.17),("Grouped","#AA4499",:utriangle,0.),("Within",:gray35,:rect,.17)]
            d=m[(m.model.==model).&(m.scheme.==scheme),:];v=[only(d.nlpd[d.dataset.==ds]) for ds in DATA]
            scatter!(ax,(1:10).+offset,v,color=color,marker=marker,markersize=13)
        end
    end
    Legend(f[2,1:2],[MarkerElement(color=:gray40,marker=x) for x in [:circle,:utriangle,:rect]],
        ["Leave one dataset out","Grouped holdout","Within-dataset fit (descriptive)"],orientation=:horizontal,framevisible=false)
    exportfig(f,"04_transfer_vs_local_fit")
    f=Figure(size=(1700,900));schemes=[("LODO","M0"),("LODO","M1"),("Grouped","M0"),("Grouped","M1"),("Within","M0"),("Within","M1")]
    for (j,metric,label) in [(1,:total_relative_error,"Total support error (%)"),(2,:links_relative_error,"Regional richness error (%)"),(3,:gap_error,"Mean gap error (percentage points)")]
        z=[100only(m[(m.dataset.==ds).&(m.scheme.==scheme).&(m.model.==model),metric]) for (scheme,model) in schemes, ds in DATA]
        bound=maximum(abs.(z));bound=max(bound,1.)
        ax=Axis(f[1,j],xticks=(1:6,[s*" "*v for (s,v) in schemes]),xticklabelrotation=pi/3,
            yticks=(1:10,[DISPLAY[x] for x in DATA]),yreversed=true)
        hm=heatmap!(ax,1:6,1:10,z,colormap=[:steelblue4,:white,:darkorange3],colorrange=(-bound,bound))
        for a in 1:6,b in 1:10
            text!(ax,a,b,text=@sprintf("%.1f",z[a,b]),align=(:center,:center),fontsize=12,color=abs(z[a,b])>.65bound ? :white : :black)
        end
        Colorbar(f[2,j],hm;vertical=false,label=label)
    end
    exportfig(f,"05_error_diagnostics")
end

mkpath(DEST)
checks()
if !("--plot-only" in ARGS)
    analyse(read_pairs())
end
figures()
println("Technical pilot saved to ",DEST)
end

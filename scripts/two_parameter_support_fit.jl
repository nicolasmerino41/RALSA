module TwoParameterSupportFit
using CSV, DataFrames, Statistics, Random, SHA, Printf, CairoMakie
include(joinpath(@__DIR__,"..","src","config.jl"))
const SOURCE=joinpath(OUT,"statistical description of local support")
const DEST=joinpath(SOURCE,"two_parameter_fit")
const B=999
const RNG=MersenneTwister(20261012)
const PALETTE=["#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499","#661100"]
const MODELS=["one_parameter","two_parameter"]
mkpath(DEST)

function geometric(r,M)
    M==1 && return [1.]
    r==0 && return vcat(1.,zeros(M-1))
    r==1 && return fill(1/M,M)
    v=exp.(log(r).*(0:M-1));v./sum(v)
end
function meanoffset(r,M)
    r==0 && return 0.
    r==1 && return (M-1)/2
    t=log(r)
    abs(M*t)<.01 && return (M-1)/2+(M^2-1)*t/12-(Float64(M)^4-1)*t^3/720
    r/(1-r)-M*exp(M*t)/(-expm1(M*t))
end
function estimate_r(target,M)
    M==1 && return NaN
    target<=0 && return 0.
    target>=(M-1)/2 && return 1.
    lo,hi=0.,1.
    for _ in 1:55
        mid=(lo+hi)/2
        if meanoffset(mid,M)<target;lo=mid;else;hi=mid;end
    end
    (lo+hi)/2
end
function fit(counts,model)
    N=length(counts);L=sum(counts);rep=L-counts[1]
    @assert L>0
    if model=="one_parameter"
        r=estimate_r(sum((0:N-1).*counts)/L,N);p=geometric(r,N)
        s=p[1];tailmean=N>1 && s<1 ? sum((2:N).*p[2:end])/(1-s) : NaN
        status=N==1 ? "N=1; decay not identifiable" : r in (0.,1.) ? "decay boundary" : "interior"
        pars=N==1 ? 0 : 1
    else
        s=counts[1]/L
        r=N<=2 || rep==0 ? NaN : estimate_r(sum((0:N-2).*counts[2:end])/rep,N-1)
        tail=N==1 ? Float64[] : geometric(isnan(r) ? 1. : r,N-1)
        p=vcat(s,(1-s).*tail)
        tailmean=N==1 || rep==0 ? NaN : sum((2:N).*tail)
        status=N==1 ? "N=1; tail absent" : N==2 ? "N=2; tail decay not identifiable" : rep==0 ? "no non-singletons; decay not identifiable" : rep<5 ? "fewer than 5 non-singletons; weakly informed decay" : r in (0.,1.) || s in (0.,1.) ? "boundary estimate" : "interior"
        pars=N==1 ? 0 : N==2 ? 1 : 2
    end
    @assert isapprox(sum(p),1;atol=1e-12)
    ll=sum(counts[k]*log(p[k]) for k in 1:N if counts[k]>0)
    aic=model=="two_parameter" && N>2 && rep==0 ? NaN : 2pars-2ll
    (s=s,r=r,p=p,tailmean=tailmean,loglik=ll,aic=aic,parameters=pars,status=status)
end
function discrepancy(counts,p)
    L=sum(counts);rep=L-counts[1]
    D=maximum(abs.(cumsum(counts)./L.-cumsum(p)))
    Dt=rep<5 || length(p)<2 || p[1]>=1 ? NaN : maximum(abs.(cumsum(counts[2:end])./rep.-cumsum(p[2:end])./(1-p[1])))
    D,Dt
end
function bootstrap(counts,model,observed)
    f=fit(counts,model);cdf=cumsum(f.p);cdf[end]=1.;L=sum(counts);N=length(counts)
    values=Float64[];tails=Float64[]
    for _ in 1:B
        sim=zeros(Int,N)
        for _ in 1:L;sim[searchsortedfirst(cdf,rand(RNG))]+=1;end
        sf=fit(sim,model);D,Dt=discrepancy(sim,sf.p)
        push!(values,D);isfinite(Dt) && push!(tails,Dt)
    end
    (p=(1+count(>=(observed[1]-1e-12),values))/(B+1),
     tail_p=isfinite(observed[2]) && !isempty(tails) ? (1+count(>=(observed[2]-1e-12),tails))/(length(tails)+1) : NaN,
     tail_replicates=length(tails))
end

@assert isapprox(meanoffset(.7,12),sum((0:11).*geometric(.7,12));atol=1e-10)
@assert isapprox(meanoffset(.999999,341),sum((0:340).*geometric(.999999,341));atol=1e-8)
@assert isnan(fit([3,0,0],"two_parameter").r)
@assert fit([2,3,0,0],"two_parameter").r==0
@assert fit([0,0,0,4],"two_parameter").r==1
@assert fit([2,3],"two_parameter").tailmean==2

linkpath=joinpath(SOURCE,"link_support.csv");sitepath=joinpath(SOURCE,"site_universe.csv")
links=CSV.read(linkpath,DataFrame);sites=CSV.read(sitepath,DataFrame;types=Dict(:site=>String))
versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
for row in eachrow(versions)
    @assert bytes2hex(sha256(read(joinpath(ROOT,row.path))))==row.sha256 "Input changed since descriptive analysis"
end
CSV.write(joinpath(DEST,"input_versions.csv"),versions)
CSV.write(joinpath(DEST,"descriptive_input_versions.csv"),DataFrame(path=[relpath(linkpath,ROOT),relpath(sitepath,ROOT)],sha256=[bytes2hex(sha256(read(linkpath))),bytes2hex(sha256(read(sitepath)))]))
summary=NamedTuple[];curves=NamedTuple[]
for ds in DATA
    v=links[links.dataset.==ds,:];N=only(unique(v.N));K=Int.(v.K);L=length(K)
    ints=CSV.read(joinpath(INPUT,ds*"_interactions.csv"),DataFrame;types=String)
    raw=combine(groupby(ints,[:consumer,:resource]),:site=>(x->length(unique(x)))=>:K)
    check=innerjoin(raw,v;on=[:consumer,:resource],makeunique=true)
    @assert nrow(check)==nrow(raw)==L && all(check.K.==check.K_1)
    @assert N==nrow(sites[sites.dataset.==ds,:]) && all(1 .<= K .<= N)
    counts=[count(==(k),K) for k in 1:N];of=counts./L;obs90=sort(K)[ceil(Int,.9L)]
    for model in MODELS
        f=fit(counts,model);D,Dt=discrepancy(counts,f.p);boot=bootstrap(counts,model,(D,Dt))
        fitted90=findfirst(>=(.9),cumsum(f.p));cut=obs90
        push!(summary,(dataset=ds,N=N,regional_links=L,non_singleton_links=L-counts[1],model=model,
            s=f.s,r=f.r,fitted_non_singleton_mean=f.tailmean,observed_non_singleton_mean=mean(K[K.>=2]),
            log_likelihood=f.loglik,parameters=f.parameters,AIC=f.aic,fit_status=f.status,
            observed_singleton_fraction=of[1],fitted_singleton_fraction=f.p[1],observed_mean=mean(K),fitted_mean=sum((1:N).*f.p),
            observed_p90=obs90,fitted_p90=fitted90,observed_fraction_above_p90=count(>(cut),K)/L,fitted_fraction_above_p90=sum(f.p[cut+1:end]),
            cdf_max_error=D,non_singleton_cdf_max_error=Dt,bootstrap_p=boot.p,tail_bootstrap_p=boot.tail_p,
            bootstrap_replicates=B,tail_bootstrap_replicates=boot.tail_replicates))
        obscc=reverse(cumsum(reverse(of)));fitcc=reverse(cumsum(reverse(f.p)))
        for k in 1:N
            push!(curves,(dataset=ds,model=model,K=k,observed_frequency=of[k],fitted_frequency=f.p[k],observed_ccdf=obscc[k],fitted_ccdf=fitcc[k]))
        end
        println(ds," ",model," r=",round(f.r;digits=4)," D=",round(D;digits=3)," p=",boot.p," tail p=",boot.tail_p)
    end
end
s=DataFrame(summary)
for ds in DATA
    a=s[(s.dataset.==ds).&(s.model.=="one_parameter"),:];b=s[(s.dataset.==ds).&(s.model.=="two_parameter"),:]
    @assert only(b.log_likelihood)>=only(a.log_likelihood)-1e-8
end
s.delta_AIC=s.AIC.-[minimum(s.AIC[s.dataset.==ds]) for ds in s.dataset]
CSV.write(joinpath(DEST,"fit_summary.csv"),s);d=DataFrame(curves);CSV.write(joinpath(DEST,"fitted_distributions.csv"),d)
set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=15,yticklabelsize=15)))
for (kind,name) in [("frequency","01_frequency_fits"),("ccdf","02_ccdf_fits")]
    f=Figure(size=(1850,900))
    for (j,ds) in enumerate(DATA)
        v=d[(d.dataset.==ds).&(d.model.=="two_parameter"),:];N=maximum(v.K);maxk=maximum(links.K[links.dataset.==ds])
        ax=Axis(f[cld(j,5),mod1(j,5)],title=DISPLAY[ds],titlesize=19,xlabel="Supporting sites, K",ylabel=mod1(j,5)==1 ? (kind=="frequency" ? "Proportion of regional links" : "Proportion of links with support ≥ K") : "",yscale=kind=="ccdf" ? log10 : identity)
        obs=kind=="frequency" ? v.observed_frequency : v.observed_ccdf
        mask=obs.>0
        if kind=="frequency"
            for r in eachrow(v[v.observed_frequency.>0,:]);lines!(ax,[r.K,r.K],[0,r.observed_frequency],color=(PALETTE[j],.45),linewidth=2);end
        end
        scatter!(ax,v.K[mask],obs[mask],color=PALETTE[j],markersize=7)
        for (model,col,style) in [("one_parameter",:gray40,:dash),("two_parameter",PALETTE[j],:solid)]
            z=d[(d.dataset.==ds).&(d.model.==model),:];y=kind=="frequency" ? z.fitted_frequency : z.fitted_ccdf
            keep=y.>0;lines!(ax,z.K[keep],y[keep],color=col,linestyle=style,linewidth=2.5)
        end
        ticks=filter(k->k<=min(N,maxk+2),[1,5,10,15,20,25,30,35,40]);ax.xticks=ticks
        xlims!(ax,.5,min(N,maxk+2)+.5)
        ylims!(ax,kind=="frequency" ? 0. : .5/only(s.regional_links[(s.dataset.==ds).&(s.model.=="two_parameter")]),kind=="frequency" ? .85 : 1.1)
    end
    Legend(f[3,1:5],[MarkerElement(color=:black,marker=:circle),LineElement(color=:gray40,linestyle=:dash),LineElement(color=:black)],
        ["Observed","One-parameter geometric","Singleton + geometric tail"];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,name*".png"),f;px_per_unit=2)
end
println("Saved to ",DEST)
end

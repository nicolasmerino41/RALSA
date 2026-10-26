module BetaPair
using CSV, DataFrames, Random, Statistics, LinearAlgebra, SHA, CairoMakie, Printf
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"../..")) : abspath(ARGS[1])
const SOURCE=joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory")
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(SOURCE,"beta_pair_conditioned_M")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const PHI=[.05,.1,.25,.5,1.,2.,4.,8.,16.,64.,256.]
const GRID=[(a=a,mu=1/(1+exp(-a)),phi=phi) for a in -4.:.5:4. for phi in PHI]
const REPS=64
const QUAD=128
const B=2500
mkpath(DEST)
logz(x)=log(-expm1(x))
lse(x)=maximum(x)+log(sum(exp.(x.-maximum(x))))
logistic(x)=1/(1+exp(-x))
include("conditional_support.jl")

# Gauss–Jacobi weights integrate over the beta density, including its endpoints.
function betanodes(mu,phi,n)
    alpha=mu*phi;beta=(1-mu)*phi;s=phi-2
    diag=zeros(n);off=zeros(n-1);diag[1]=(alpha-beta)/phi
    for j in 1:n-1;diag[j+1]=((alpha-1)^2-(beta-1)^2)/((2j+s)*(2j+s+2));end
    off[1]=sqrt(4alpha*beta/(phi^2*(phi+1)))
    for j in 2:n-1
        off[j]=sqrt(4j*(j+beta-1)*(j+alpha-1)*(j+s)/((2j+s)^2*((2j+s)^2-1)))
    end
    eig=eigen(SymTridiagonal(diag,off));p=(eig.values.+1)./2;w=vec(eig.vectors[1,:].^2)
    @assert all(0 .<p .<1) && abs(sum(w)-1)<1e-10 && abs(dot(w,p)-mu)<1e-10
    @assert abs(dot(w,(p.-mu).^2)-mu*(1-mu)/(phi+1))<1e-10
    p,w
end
function betapmfs(caps,mu,phi,n=QUAD)
    p,w=betanodes(mu,phi,n);lp=log.(p);lq=log1p.(-p);out=Dict{Int,Vector{Float64}}()
    for c in caps
        v=zeros(c);lc=0.
        for k in 1:c
            lc+=log((c-k+1)/k)
            v[k]=sum(w[t]*exp(lc+k*lp[t]+(c-k)*lq[t]-logz(c*lq[t])) for t in eachindex(p))
        end
        @assert abs(sum(v)-1)<1e-9
        out[c]=v/sum(v)
    end
    out
end
function summaries(O,N,K,T)
    pred=zeros(N);v=zeros(N)
    for (i,js) in enumerate(O)
        p=K[i]/length(js)
        for j in js;pred[j]+=p;v[j]+=p*(1-p);end
    end
    cluster=mean(((T.-pred).^2)./(v.+1))
    [log(mean(K)+.1),mean(K.==1),log(mean(K.^2)+.1),log1p(cluster)]
end
function simulate(rng,O,N,P)
    K=zeros(Int,length(O));T=zeros(Int,N)
    for (i,js) in enumerate(O)
        c=length(js);u=rand(rng);v=0.;k=c
        for j in 1:c;v+=P[c][j];if u<=v;k=j;break;end;end
        K[i]=k
        positions=randperm(rng,c)
        for j in @view positions[1:k];T[js[j]]+=1;end
    end
    summaries(O,N,K,T)
end
function loadinputs()
    ver=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
    for r in eachrow(ver);@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256;end
    ver.path=String.(ver.path);ver.sha256=String.(ver.sha256)
    for file in ["quality_control.csv","fitted_parameters.csv","prediction_metrics.csv","conditioned_M/summary.csv","conditioned_M/curves.csv","conditioned_M/support_frequencies.csv","conditioned_M/input_versions.csv"]
        path=joinpath(SOURCE,file);push!(ver,(relpath(path,ROOT),bytes2hex(sha256(read(path)))))
    end
    for r in eachrow(CSV.read(joinpath(SOURCE,"conditioned_M/input_versions.csv"),DataFrame));@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256;end
    qc=CSV.read(joinpath(SOURCE,"quality_control.csv"),DataFrame)
    sites=CSV.read(joinpath(ROOT,"outputs/statistical description of local support/site_universe.csv"),DataFrame;types=Dict(:site=>String))
    systems=NamedTuple[]
    for ds in DATA
        g=qc[qc.dataset.==ds,:];N=only(unique(g.N));ss=sites.site[sites.dataset.==ds];siteid=Dict(s=>i for (i,s) in enumerate(ss));@assert length(ss)==N
        occ=CSV.read(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv"),DataFrame;types=String)
        pres=Dict((v.trophic_level[1],v.species[1])=>Set(v.site) for v in groupby(occ,[:trophic_level,:species]))
        O=[sort([siteid[s] for s in intersect(pres[("consumer",r.consumer)],pres[("resource",r.resource)])]) for r in eachrow(g)]
        K=Int.(g.K);C=length.(O);@assert C==g.C && all(1 .<=K .<=C .<=N)
        ints=CSV.read(joinpath(ROOT,"data/processed",ds*"_interactions.csv"),DataFrame;types=String)
        pairid=Dict((r.consumer,r.resource)=>i for (i,r) in enumerate(eachrow(g)))
        A=falses(length(K),N)
        for r in eachrow(ints);A[pairid[(r.consumer,r.resource)],siteid[r.site]]=true;end
        @assert vec(sum(A;dims=2))==K
        obs=summaries(O,N,K,vec(sum(A;dims=1)))
        push!(systems,(name=ds,N=N,C=C,K=K,O=O,obs=obs))
    end
    CSV.write(joinpath(DEST,"input_versions.csv"),unique(ver));systems
end
function library(systems)
    file=joinpath(DEST,"simulation_library.csv")
    signature=bytes2hex(sha256(string(GRID,REPS,QUAD,[(s.name,s.N,s.O) for s in systems])))
    if isfile(file)
        @assert strip(read(joinpath(DEST,"library_signature.txt"),String))==signature
        return CSV.read(file,DataFrame)
    end
    rows=NamedTuple[];caps=sort(unique(vcat([s.C for s in systems]...)))
    rngs=[MersenneTwister(224000+i) for i in eachindex(systems)]
    for (gi,p) in enumerate(GRID)
        P=betapmfs(caps,p.mu,p.phi)
        for (di,s) in enumerate(systems),rep in 1:REPS
            x=simulate(rngs[di],s.O,s.N,P)
            push!(rows,(dataset=s.name,grid=gi,replicate=rep,x1=x[1],x2=x[2],x3=x[3],x4=x[4]))
        end
        gi%17==0 && (println("Beta library ",gi,"/",length(GRID));flush(stdout))
    end
    tab=DataFrame(rows);CSV.write(file,tab);write(joinpath(DEST,"library_signature.txt"),signature*"\n");tab
end
function fitmodels(systems,lib)
    refs=Dict()
    for g in groupby(lib,[:dataset,:grid])
        X=Matrix(g[:,[:x1,:x2,:x3,:x4]]);mu=vec(mean(X;dims=1));V=cov(X);V=.8V+.2Diagonal(diag(V))+1e-5I
        refs[(g.dataset[1],g.grid[1])]=(mu,inv(V),logdet(V))
    end
    fits=NamedTuple[];scores=NamedTuple[]
    for target in systems
        train=filter(s->s.name!=target.name,systems);values=Float64[]
        for (gi,p) in enumerate(GRID)
            score=mean(begin mu,Q,ld=refs[(s.name,gi)];z=s.obs.-mu;(dot(z,Q*z)+ld)/4 end for s in train)
            push!(values,score);push!(scores,(held_out_dataset=target.name,grid=gi,a=p.a,mu=p.mu,phi=p.phi,score=score))
        end
        gi=argmin(values);p=GRID[gi]
        push!(fits,(dataset=target.name,grid=gi,a=p.a,mu=p.mu,phi=p.phi,score=values[gi],boundary=abs(p.a)==4 || p.phi in extrema(PHI),training_datasets=join([s.name for s in train],";")))
    end
    CSV.write(joinpath(DEST,"fitted_parameters.csv"),DataFrame(fits));CSV.write(joinpath(DEST,"training_grid_scores.csv"),DataFrame(scores));DataFrame(fits)
end
function checks()
    rows=NamedTuple[]
    # For C=2 and a uniform beta, P(K=1)=2-2log(2).
    p=betapmfs([1,2],.5,2.,QUAD)
    err=abs(p[2][1]-(2-2log(2)));@assert err<1e-10 && p[1]==[1.]
    push!(rows,(check="Uniform beta analytic C=2 probability",error=err,tolerance=1e-10))
    C=[2,3,4];M=6;P=betapmfs(C,.4,1.)
    D,W=dp([log.(P[c]) for c in C],M);f=exactfreq(M,D,W,4)
    kk=[collect(k) for k in Iterators.product((1:c for c in C)...) if sum(k)==M]
    w=[prod(P[C[i]][k[i]] for i in 1:3) for k in kk];w/=sum(w);truth=zeros(4)
    for (k,v) in zip(kk,w),x in k;truth[x]+=v/3;end
    err=maximum(abs.(f.-truth));@assert err<1e-10
    push!(rows,(check="Fixed-M beta DP versus enumeration",error=err,tolerance=1e-10))
    CSV.write(joinpath(DEST,"algorithm_checks.csv"),DataFrame(rows))
end
function evaluate(systems,fits)
    prior=CSV.read(joinpath(SOURCE,"conditioned_M/summary.csv"),DataFrame)
    pcurves=CSV.read(joinpath(SOURCE,"conditioned_M/curves.csv"),DataFrame)
    pfreq=CSV.read(joinpath(SOURCE,"conditioned_M/support_frequencies.csv"),DataFrame)
    rows=NamedTuple[];curves=NamedTuple[];freqs=NamedTuple[];draws=NamedTuple[];quadrature=NamedTuple[]
    for (di,s) in enumerate(systems)
        fit=only(eachrow(fits[fits.dataset.==s.name,:]));@assert !(s.name in split(fit.training_datasets,";"))
        caps=unique(s.C);P=betapmfs(caps,fit.mu,fit.phi,QUAD);fine=betapmfs(caps,fit.mu,fit.phi,2QUAD)
        err=maximum(maximum(abs.(P[c].-fine[c])) for c in caps);@assert err<1e-8
        push!(quadrature,(dataset=s.name,nodes=QUAD,comparison_nodes=2QUAD,maximum_probability_difference=err,passed=true))
        L=length(s.K);N=s.N;M=sum(s.K);H=gapmatrix(N);idx=1:fld(4N,5)+1;coef=vec(mean(H[idx,:];dims=1));emp=vec(mean(H[:,s.K];dims=2))
        uncond=zeros(N)
        for c in s.C;uncond[1:c].+=P[c]/L;end
        D,W=dp([log.(P[c]) for c in s.C],M);f=exactfreq(M,D,W,N);pred=H*f;errors=100(pred[idx].-emp[idx]);rng=MersenneTwister(993000+di)
        gapdraws=zeros(B);singledraws=zeros(B)
        for b in 1:B
            K=drawK(rng,M,D,W);@assert sum(K)==M && all(1 .<=K .<=s.C)
            gapdraws[b]=100mean(coef[K]);singledraws[b]=mean(K.==1)
            push!(draws,(dataset=s.name,replicate=b,gap_pp=gapdraws[b],singleton_fraction=singledraws[b]))
        end
        @assert abs(mean(gapdraws)-100mean(pred[idx]))<6std(gapdraws)/sqrt(B)
        lo,hi=quantile(gapdraws,[.025,.975]);sl,sh=quantile(singledraws,[.025,.975])
        common=only(eachrow(prior[(prior.dataset.==s.name).&(prior.model.=="Common"),:]));pair=only(eachrow(prior[(prior.dataset.==s.name).&(prior.model.=="Pair"),:]))
        ref=pcurves[(pcurves.dataset.==s.name).&(pcurves.model.=="Pair"),:];@assert ref.m==collect(0:N) && maximum(abs.(100emp.-ref.empirical_gap_pp))<1e-8
        push!(rows,(dataset=s.name,N=N,L=L,M=M,mu=fit.mu,phi=fit.phi,boundary=fit.boundary,empirical_mean_gap_pp=100mean(emp[idx]),predicted_mean_gap_pp=100mean(pred[idx]),lower_gap_pp=lo,upper_gap_pp=hi,mae_pp=mean(abs.(errors)),maximum_error_pp=maximum(abs.(errors)),signed_error_pp=mean(errors),common_mae_pp=common.mae_pp,two_class_mae_pp=pair.mae_pp,mae_change_vs_two_class_pp=mean(abs.(errors))-pair.mae_pp,empirical_singleton_percent=100mean(s.K.==1),predicted_singleton_percent=100f[1],lower_singleton_percent=100sl,upper_singleton_percent=100sh,singleton_error_pp=100(f[1]-mean(s.K.==1)),two_class_singleton_error_pp=pair.singleton_error_pp,unconditioned_beta_mae_pp=100mean(abs.((H*uncond)[idx].-emp[idx])),unconditioned_mean_support=dot(1:N,uncond)))
        for m in 0:N;push!(curves,(dataset=s.name,m=m,removal_fraction=m/N,empirical_gap_pp=100emp[m+1],beta_gap_pp=100pred[m+1],two_class_gap_pp=ref.predicted_gap_pp[m+1],common_gap_pp=only(pcurves.predicted_gap_pp[(pcurves.dataset.==s.name).&(pcurves.model.=="Common").&(pcurves.m.==m)]),unconditioned_beta_gap_pp=100dot(H[m+1,:],uncond)));end
        for k in 1:N;push!(freqs,(dataset=s.name,K=k,empirical=mean(s.K.==k),beta=f[k],two_class=only(pfreq.predicted_frequency[(pfreq.dataset.==s.name).&(pfreq.model.=="Pair").&(pfreq.K.==k)])));end
        println(s.name,": beta fixed-M MAE ",round(mean(abs.(errors));digits=3),"; two-class ",round(pair.mae_pp;digits=3),"; mu ",round(fit.mu;digits=3)," phi ",fit.phi);flush(stdout)
    end
    for (file,v) in [("summary",rows),("curves",curves),("support_frequencies",freqs),("draws",draws),("quadrature_checks",quadrature)];CSV.write(joinpath(DEST,file*".csv"),DataFrame(v));end
end
function plots()
    t=CSV.read(joinpath(DEST,"summary.csv"),DataFrame);curves=CSV.read(joinpath(DEST,"curves.csv"),DataFrame)
    p=CSV.read(joinpath(SOURCE,"conditioned_M/summary.csv"),DataFrame)
    beta="#228877";pair="#CC6677";common="#888888"
    set_theme!(Theme(font="Arial",fontsize=20,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=17,yticklabelsize=17)))
    f=Figure(size=(1500,800))
    for (j,kind,label) in [(1,"gap","Mean gap (percentage points)"),(2,"singleton","Singleton links (%)")]
        ax=Axis(f[1,j];xlabel=label,yticks=(1:10,replace.(DATA,"_"=>" ")))
        for (i,ds) in enumerate(DATA)
            r=only(eachrow(t[t.dataset.==ds,:]))
            for (model,col,offset) in [("Common",common,-.22),("Pair",pair,0.),("Beta",beta,.22)]
                z=model=="Beta" ? r : only(eachrow(p[(p.dataset.==ds).&(p.model.==model),:]))
                lo=kind=="gap" ? z.lower_gap_pp : z.lower_singleton_percent;hi=kind=="gap" ? z.upper_gap_pp : z.upper_singleton_percent
                mu=kind=="gap" ? z.predicted_mean_gap_pp : z.predicted_singleton_percent
                lines!(ax,[lo,hi],[i+offset,i+offset];color=col,linewidth=3);scatter!(ax,[mu],[i+offset];color=col,markersize=9)
            end
            obs=kind=="gap" ? r.empirical_mean_gap_pp : r.empirical_singleton_percent
            scatter!(ax,[obs],[i];color=:black,marker=:diamond,markersize=13)
        end
        ylims!(ax,.4,10.6);ax.yreversed=true
    end
    Legend(f[2,1:2],[LineElement(color=common),LineElement(color=pair),LineElement(color=beta),MarkerElement(color=:black,marker=:diamond)],["Common / allocation null","Two-class pair","Beta pair","Empirical"];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"01_allocations.png"),f;px_per_unit=2)
    f=Figure(size=(1450,790))
    for (j,kind,label) in [(1,"gap","Gap MAE (percentage points)"),(2,"singleton","Absolute singleton error (pp)")]
        ax=Axis(f[1,j];xlabel=label,yticks=(1:11,vcat(replace.(DATA,"_"=>" "),["Equal-dataset mean"])))
        a=kind=="gap" ? t.two_class_mae_pp : abs.(t.two_class_singleton_error_pp);b=kind=="gap" ? t.mae_pp : abs.(t.singleton_error_pp)
        a=vcat(a,mean(a));b=vcat(b,mean(b));hlines!(ax,[10.5];color="#DDDDDD")
        for i in 1:11
            lines!(ax,[a[i],b[i]],[i,i];color="#BBBBBB",linewidth=2)
        end
        scatter!(ax,a,1:11;color=pair,markersize=vcat(fill(11,10),17),marker=:diamond)
        scatter!(ax,b,1:11;color=beta,markersize=vcat(fill(10,10),16))
        ylims!(ax,.4,11.6);ax.yreversed=true
    end
    Legend(f[2,1:2],[MarkerElement(color=pair,marker=:diamond),MarkerElement(color=beta,marker=:circle)],["Two-class pair","Beta pair"];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"02_two_class_vs_beta.png"),f;px_per_unit=2)
    f=Figure(size=(1700,820))
    for (i,ds) in enumerate(DATA)
        ax=Axis(f[cld(i,5),mod1(i,5)];title=replace(ds,"_"=>" "),xlabel="Sites removed (%)",ylabel=mod1(i,5)==1 ? "Gap (pp)" : "")
        v=curves[(curves.dataset.==ds).&(curves.removal_fraction.<=.8),:]
        for (field,col,style) in [(:common_gap_pp,common,:solid),(:two_class_gap_pp,pair,:solid),(:beta_gap_pp,beta,:solid),(:empirical_gap_pp,:black,:dash)]
            lines!(ax,100v.removal_fraction,v[!,field];color=col,linestyle=style,linewidth=2)
        end
        xlims!(ax,0,80);ylims!(ax,0,40)
    end
    Legend(f[3,1:5],[LineElement(color=common),LineElement(color=pair),LineElement(color=beta),LineElement(color=:black,linestyle=:dash)],["Common / allocation null","Two-class pair","Beta pair","Empirical"];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"03_gap_curves.png"),f;px_per_unit=2)
end
function report()
    t=CSV.read(joinpath(DEST,"summary.csv"),DataFrame);curves=CSV.read(joinpath(DEST,"curves.csv"),DataFrame);freq=CSV.read(joinpath(DEST,"support_frequencies.csv"),DataFrame)
    base=CSV.read(joinpath(SOURCE,"conditioned_M/summary.csv"),DataFrame);rows=NamedTuple[]
    for (model,mae,se) in [("Common",t.common_mae_pp,base.singleton_error_pp[base.model.=="Common"]),("Two-class pair",t.two_class_mae_pp,t.two_class_singleton_error_pp),("Beta pair",t.mae_pp,t.singleton_error_pp)]
        push!(rows,(model=model,mean_mae_pp=mean(mae),median_mae_pp=median(mae),mean_absolute_singleton_error_pp=mean(abs.(se)),datasets_improved_vs_common=count(mae.<t.common_mae_pp),datasets_improved_vs_two_class=count(mae.<t.two_class_mae_pp)))
    end
    CSV.write(joinpath(DEST,"across_dataset_summary.csv"),DataFrame(rows));checks=NamedTuple[]
    for r in eachrow(t)
        v=curves[curves.dataset.==r.dataset,:];f=freq[freq.dataset.==r.dataset,:];mask=v.removal_fraction.<=.8
        @assert v.m==collect(0:r.N) && abs(sum(f.beta)-1)<1e-8 && abs(dot(f.K,f.beta)-r.M/r.L)<1e-8
        @assert abs(mean(abs.(v.beta_gap_pp[mask].-v.empirical_gap_pp[mask]))-r.mae_pp)<1e-8
        @assert maximum(abs.(100gapmatrix(r.N)*f.beta.-v.beta_gap_pp))<1e-8
        push!(checks,(dataset=r.dataset,fixed_M=true,curve_from_support=true,metric_averaging=true))
    end
    CSV.write(joinpath(DEST,"result_checks.csv"),DataFrame(checks))
    fmt(x)=@sprintf("%.3f",x)
    open(joinpath(DEST,"FINDINGS.md"),"w") do io
        println(io,"# Continuous beta versus two-class pair heterogeneity\n")
        println(io,"Both approaches use known regional links, positive support, the same confirmed opportunities and exact observed M. Beta parameters are trained on the other nine datasets using the same four-summary synthetic-likelihood criterion as the saved two-class fit. No target-specific parameter tuning was performed.\n")
        println(io,"| Model | Equal-dataset mean gap MAE (pp) | Median MAE (pp) | Mean absolute singleton error (pp) |\n|---|---:|---:|---:|")
        for r in rows;println(io,"| $(r.model) | $(fmt(r.mean_mae_pp)) | $(fmt(r.median_mae_pp)) | $(fmt(r.mean_absolute_singleton_error_pp)) |");end
        println(io,"\nBeta improves on Common in $(count(t.mae_pp.<t.common_mae_pp))/10 datasets and on the two-class model in $(count(t.mae_pp.<t.two_class_mae_pp))/10. These comparisons give each dataset equal weight.\n")
        println(io,"The average gap error is slightly lower for Beta, but the dataset wins are evenly split and its mean absolute singleton error is higher. Both heterogeneity descriptions improve substantially on Common. The defensible conclusion is that the gap improvement survives replacing two classes with continuous pair heterogeneity; this does not establish actual ecological classes or a uniquely correct distribution. Beta's largest gap improvements occur in the two Gottin datasets, which share sites, so they are not independent replications of its advantage.\n")
        println(io,"| Dataset | Two-class gap MAE | Beta gap MAE | Beta mu | Beta concentration |\n|---|---:|---:|---:|---:|")
        for r in eachrow(t);println(io,"| $(replace(r.dataset,"_"=>" ")) | $(fmt(r.two_class_mae_pp)) | $(fmt(r.mae_pp)) | $(fmt(r.mu)) | $(r.phi) |");end
        println(io,"\nBoundary selections: $(count(t.boundary))/10. Fits use a prespecified finite grid, not unconstrained maximum likelihood. Parameters minimise training summary discrepancy, not held-out gap errors. The library has 64 replicates per dataset/grid point; simulation and grid resolution can influence close rankings.\n")
        println(io,"Nine folds select mu=0.5 and phi=0.5, corresponding to Beta(0.25,0.25), a continuous U-shaped distribution concentrated toward low and high probabilities. Garraf PP2 selects mu≈0.378 and phi=2. Thus the continuous fit still favours substantial probability heterogeneity; it does not imply that probabilities are concentrated around one common value. These are latent probabilities before the positivity and total-support conditioning.\n")
        println(io,"The beta distribution is placed among known pairs, matching the existing two-class convention. The support distribution integrates a positive-truncated binomial over that beta density. This is a beta mixture of positive binomials, not the conventional beta-binomial truncated after mixing. The distinction is deliberate: changing it would change selection as well as the shape of pair heterogeneity. See README for the equation.\n")
        println(io,"The comparison can show whether the gap result survives replacing two classes with a continuous distribution, and whether one better describes these outputs. It cannot establish two genuine ecological classes, a true beta distribution, or a unique cause of heterogeneity. There are only ten related datasets and the test conditions on observed M.\n")
        println(io,"All input fingerprints, opportunity/support constraints, fixed-M marginals and curve/metric identities were checked. Gaussian quadrature was checked against an analytic example and doubled in resolution for each selected fit. Exact finite enumeration validates the conditional dynamic program. All conditional means use dynamic programming; 2,500 independent draws per dataset supply central 95% allocation intervals, not confidence intervals.\n")
        println(io,"Confirmed-restoration flags, interaction-derived occupancy, excluded uncertain records, shared Gottin sites and unresolved shared-collection overlap remain unchanged. This is exploratory conditional prediction for known links, not independent landscape validation or discovery of unobserved interactions.")
    end
    CSV.write(joinpath(DEST,"code_versions.csv"),DataFrame(file=["beta_pair.jl","conditional_support.jl"],sha256=[bytes2hex(sha256(read(joinpath(@__DIR__,file)))) for file in ["beta_pair.jl","conditional_support.jl"]]))
end
if !(length(ARGS)>=3 && ARGS[3] in ["plot","report"])
    checks();systems=loadinputs();lib=library(systems);fits=fitmodels(systems,lib);evaluate(systems,fits)
    for file in ["metadata_relationships.csv","metadata_source_versions.csv"];cp(joinpath(SOURCE,file),joinpath(DEST,file);force=true);end
end
if !(length(ARGS)>=3 && ARGS[3]=="report");plots();end
report()
end

module ConditionedM
using CSV, DataFrames, Statistics, Random, SHA, CairoMakie, LinearAlgebra, Printf
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"../..")) : abspath(ARGS[1])
const SOURCE=joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory")
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(SOURCE,"conditioned_M")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const B=2500
const CHAINS=4
const BURN=1000
const STEPS=4000
const THIN=5
mkpath(DEST)
logistic(x)=1/(1+exp(-x))
logz(x)=log(-expm1(x))
lse(x)=maximum(x)+log(sum(exp.(x.-maximum(x))))
function gapmatrix(N)
    H=zeros(N+1,N)
    for m in 0:N,k in 1:N
        H[m+1,k]=m/N-(m<k ? 0. : prod((m-j)/(N-j) for j in 0:k-1))
    end
    H
end
function pairweights(C,a,h)
    [begin
        qs=logistic.([a-h,a+h])
        [sum(log((c-j+1)/j) for j in 1:k)+lse([k*log(q)+(c-k)*log1p(-q)-logz(c*log1p(-q)) for q in qs])-log(2) for k in 1:c]
    end for c in C]
end
function dp(logw,M)
    C=length.(logw);L=length(C)
    low=-60.;high=60.
    for _ in 1:65
        t=(low+high)/2
        mu=sum(begin v=w.+t.*collect(1:length(w));p=exp.(v.-lse(v));sum(p.*collect(1:length(w))) end for w in logw)
        if mu>M;high=t;else;low=t;end
    end
    tilt=(low+high)/2
    P=[exp.(w.+tilt.*collect(1:length(w)).-lse(w.+tilt.*collect(1:length(w)))) for w in logw]
    D=zeros(L+1,M+1);D[1,1]=1.;cap=0
    for i in 1:L
        cap+=C[i]
        for total in i:min(cap,M),k in 1:min(C[i],total-i+1)
            D[i+1,total+1]+=P[i][k]*D[i,total-k+1]
        end
        D[i+1,:]./=maximum(@view D[i+1,:])
    end
    @assert D[end,M+1]>0
    D,P
end
function drawK(rng,M,D,P)
    K=zeros(Int,length(P));left=M
    for i in length(P):-1:1
        w=[P[i][k]*D[i,left-k+1] for k in 1:min(length(P[i]),left-i+1)]
        u=rand(rng)*sum(w);v=0.;chosen=length(w)
        for k in eachindex(w);v+=w[k];if u<=v;chosen=k;break;end;end
        K[i]=chosen;left-=chosen
    end
    @assert left==0
    K
end
function exactfreq(M,D,P,N)
    state=zeros(M+1);state[M+1]=1.;freq=zeros(N)
    for i in length(P):-1:1
        prior=zeros(M+1)
        for total in i:M
            state[total+1]<1e-300 && continue
            upto=min(length(P[i]),total-i+1)
            den=sum(P[i][k]*D[i,total-k+1] for k in 1:upto)
            den==0 && continue
            for k in 1:upto
                mass=state[total+1]*P[i][k]*D[i,total-k+1]/den
                freq[k]+=mass;prior[total-k+1]+=mass
            end
        end
        state=prior
    end
    @assert isapprox(sum(freq),length(P);atol=1e-7)
    @assert isapprox(sum(freq.*collect(1:N)),M;atol=1e-7)
    freq/length(P)
end

# Joint MCMC targets the original positive-row model conditional on exact M.
function sitechain(O,N,M,a,h,H,seed,chain;burn=BURN,steps=STEPS,thin=THIN)
    rng=MersenneTwister(seed);L=length(O);C=length.(O)
    ri=Int[];sj=Int[];rowcells=[Int[] for _ in 1:L];neighbors=[Int[] for _ in 1:N]
    for i in 1:L,j in O[i]
        push!(ri,i);push!(sj,j);push!(rowcells[i],length(ri));push!(neighbors[j],i)
    end
    E=length(ri);x=falses(E);K=ones(Int,L)
    for cells in rowcells;x[rand(rng,cells)]=true;end
    free=findall(.!x);shuffle!(rng,free);x[free[1:M-L]].=true
    K=[count(x[cells]) for cells in rowcells]
    onesidx=findall(x);zerosidx=findall(.!x)
    z=chain==1 ? fill(-1,N) : chain==2 ? fill(1,N) : rand(rng,[-1,1],N)
    q=logistic.(a .+ h.*z);lp=log.(q);lq=log1p.(-q)
    rowzero=[sum(lq[js]) for js in O];T=[count(x[t] for t in eachindex(x) if sj[t]==j) for j in 1:N]
    occ=length.(neighbors);idx=1:fld(4N,5)+1;coef=vec(mean(H[idx,:];dims=1))
    records=NamedTuple[];freq=zeros(N);acceptswap=0;acceptsign=0;saved=0
    for sweep in 1:burn+steps
        if !isempty(zerosidx)
            for _ in 1:M
                u=rand(rng,eachindex(onesidx));v=rand(rng,eachindex(zerosidx));e=onesidx[u];f=zerosidx[v]
                i=ri[e];ii=ri[f];j=sj[e];jj=sj[f]
                i!=ii && K[i]==1 && continue
                if log(rand(rng))<h*(z[jj]-z[j])
                    x[e]=false;x[f]=true;onesidx[u]=f;zerosidx[v]=e
                    K[i]-=1;K[ii]+=1;T[j]-=1;T[jj]+=1;acceptswap+=1
                end
            end
        end
        for j in randperm(rng,N)
            qq=logistic(a-h*z[j]);newlp=log(qq);newlq=log1p(-qq);delta=newlq-lq[j]
            ratio=T[j]*(newlp-lp[j])+(occ[j]-T[j])*delta
            for i in neighbors[j];ratio-=logz(rowzero[i]+delta)-logz(rowzero[i]);end
            if log(rand(rng))<ratio
                z[j]*=-1;lp[j]=newlp;lq[j]=newlq
                for i in neighbors[j];rowzero[i]+=delta;end
                acceptsign+=1
            end
        end
        if sweep>burn && (sweep-burn)%thin==0
            @assert sum(K)==M && all(1 .<=K .<=C)
            saved+=1
            for k in K;freq[k]+=1;end
            push!(records,(chain=chain,draw=saved,gap_pp=100mean(coef[K]),singleton_fraction=mean(K.==1),high_site_fraction=mean(z.==1)))
        end
    end
    records,freq/(saved*L),(swap_acceptance=acceptswap/((burn+steps)*M),sign_acceptance=acceptsign/((burn+steps)*N))
end
function diagnostics(X)
    n,c=size(X);half=fld(n,2);Y=hcat(X[1:half,:],X[end-half+1:end,:]);n,c=size(Y)
    W=mean(var(Y;dims=1));between=n*var(vec(mean(Y;dims=1)));v=(n-1)/n*W+between/n
    W==0 && return (rhat=1.,ess=Float64(n*c),mcse=0.)
    rho=Float64[]
    for lag in 0:n-1
        ac=mean([sum((Y[1:n-lag,j].-mean(Y[:,j])).*(Y[1+lag:n,j].-mean(Y[:,j])))/n for j in 1:c])
        push!(rho,1-(W-ac)/v)
    end
    pairs=Float64[]
    for t in 1:2:length(rho)-1
        p=rho[t]+rho[t+1];p<0 && break
        push!(pairs,isempty(pairs) ? p : min(p,last(pairs)))
    end
    ess=min(n*c,n*c/max(1.,-1+2sum(pairs)))
    (rhat=sqrt(v/W),ess=ess,mcse=sqrt(v/ess))
end
function samplerchecks()
    O=[[1,2],[2,3],[1,2,3]];C=length.(O);N=3;M=5;a=-.5;h=.7;H=gapmatrix(N)
    configs=[collect(k) for k in Iterators.product((1:c for c in C)...) if sum(k)==M]
    errs=Float64[]
    for hh in [0.,h]
        w=pairweights(C,a,hh);D,P=dp(w,M)
        truth=[exp(sum(w[i][k[i]] for i in eachindex(C))) for k in configs];truth/=sum(truth)
        f=zeros(N)
        for (k,p) in zip(configs,truth),v in k;f[v]+=p/length(C);end
        push!(errs,maximum(abs.(f-exactfreq(M,D,P,N))))
        if hh==0
            neutral=[prod(binomial(c,k) for (c,k) in zip(C,v)) for v in configs];neutral=neutral/sum(neutral)
            push!(errs,maximum(abs.(neutral-truth)))
        end
    end
    cells=[(i,j) for i in eachindex(O) for j in O[i]];weights=Float64[];gaps=Float64[];single=Float64[]
    coef=vec(mean(H[1:fld(4N,5)+1,:];dims=1))
    for z in Iterators.product(ntuple(_->[-1,1],N)...),mask in 0:2^length(cells)-1
        count_ones(mask)==M || continue
        K=zeros(Int,length(O));q=logistic.(a .+h.*collect(z));ll=0.
        for (t,(i,j)) in enumerate(cells)
            hit=(mask>>(t-1))&1;K[i]+=hit;ll+=hit==1 ? log(q[j]) : log1p(-q[j])
        end
        all(K.>0) || continue
        ll-=sum(logz(sum(log1p(-q[j]) for j in js)) for js in O)
        push!(weights,exp(ll));push!(gaps,100mean(coef[K]));push!(single,mean(K.==1))
    end
    weights/=sum(weights);truthgap=dot(weights,gaps);truthsingle=dot(weights,single)
    rr=NamedTuple[]
    for ch in 1:4
        r,_,_=sitechain(O,N,M,a,h,H,9900+ch,ch;burn=1000,steps=20000,thin=5);append!(rr,r)
    end
    tab=DataFrame(rr);out=NamedTuple[]
    for (field,truth) in [(:gap_pp,truthgap),(:singleton_fraction,truthsingle)]
        X=hcat([tab[tab.chain.==ch,field] for ch in 1:4]...);d=diagnostics(X);error=abs(mean(X)-truth)
        @assert error<6d.mcse+1e-5 && d.rhat<1.02
        push!(out,(check="Site exact enumeration: $field",error=error,tolerance=6d.mcse+1e-5,passed=true))
    end
    @assert maximum(errs)<1e-10
    push!(out,(check="Pair exact marginals and Common-null equivalence",error=maximum(errs),tolerance=1e-10,passed=true))
    CSV.write(joinpath(DEST,"sampler_checks.csv"),DataFrame(out))
end
function main()
    samplerchecks();println("Sampler verification passed");flush(stdout)
    ver=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
    for r in eachrow(ver);@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256;end
    ver.path=String.(ver.path);ver.sha256=String.(ver.sha256)
    function readsource(file)
        path=joinpath(SOURCE,file);push!(ver,(relpath(path,ROOT),bytes2hex(sha256(read(path)))));CSV.read(path,DataFrame)
    end
    qc=readsource("quality_control.csv");fits=readsource("fitted_parameters.csv");old=readsource("prediction_metrics.csv")
    neutral=readsource("neutral_allocation/summary.csv");ndraw=readsource("neutral_allocation/draws.csv")
    sites=CSV.read(joinpath(ROOT,"outputs/statistical description of local support/site_universe.csv"),DataFrame;types=Dict(:site=>String))
    summaries=NamedTuple[];curves=NamedTuple[];drawrows=NamedTuple[];diagrows=NamedTuple[];freqrows=NamedTuple[]
    for (di,ds) in enumerate(DATA)
        g=qc[qc.dataset.==ds,:];N=only(unique(g.N));C=Int.(g.C);K=Int.(g.K);M=sum(K);L=length(K);H=gapmatrix(N);idx=1:fld(4N,5)+1
        coef=vec(mean(H[idx,:];dims=1));emp=vec(mean(H[:,K];dims=2));empmean=100mean(emp[idx]);empsingle=mean(K.==1)
        ss=sites.site[sites.dataset.==ds];siteid=Dict(s=>j for (j,s) in enumerate(ss));@assert length(ss)==N
        occ=CSV.read(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv"),DataFrame;types=String)
        pres=Dict((v.trophic_level[1],v.species[1])=>Set(v.site) for v in groupby(occ,[:trophic_level,:species]))
        O=[sort([siteid[s] for s in intersect(pres[("consumer",r.consumer)],pres[("resource",r.resource)])]) for r in eachrow(g)]
        @assert length.(O)==C && all(1 .<=K .<=C .<=N)
        for model in ["Common","Pair","Site"]
            p=only(eachrow(fits[(fits.dataset.==ds).&(fits.scope.=="LODO").&(fits.model.==model),:]))
            @assert Set(split(p.training_datasets,";"))==Set(filter(!=(ds),DATA))
            rng=MersenneTwister(814000+di);records=NamedTuple[];f=zeros(N);mcse=0.;rhat=1.;ess=Float64(B)
            if model!="Site"
                D,P=dp(pairweights(C,p.a,p.h),M);f=exactfreq(M,D,P,N)
                if model=="Pair"
                    for b in 1:B
                        kk=drawK(rng,M,D,P);@assert sum(kk)==M && all(1 .<=kk .<=C)
                        push!(records,(chain=0,draw=b,gap_pp=100mean(coef[kk]),singleton_fraction=mean(kk.==1),high_site_fraction=NaN))
                    end
                else
                    n=ndraw[ndraw.dataset.==ds,:]
                    for (b,r) in enumerate(eachrow(n));push!(records,(chain=0,draw=b,gap_pp=r.gap_pp,singleton_fraction=r.singleton_fraction,high_site_fraction=NaN));end
                    mu=100dot(coef,f);@assert abs(mu-mean(n.gap_pp))<6std(n.gap_pp)/sqrt(nrow(n))
                end
            else
                extended=ds in ["Gottin_PP","Montseny"]
                burn=extended ? 5000 : BURN;steps=extended ? 40000 : STEPS;thin=extended ? 20 : THIN
                for chain in 1:CHAINS
                    r,ff,accept=sitechain(O,N,M,p.a,p.h,H,920000+100di+chain,chain;burn,steps,thin)
                    append!(records,r);f.+=ff/CHAINS
                    push!(diagrows,(dataset=ds,quantity="acceptance_chain_$chain",rhat=NaN,ess=NaN,mcse=NaN,swap_acceptance=accept.swap_acceptance,sign_acceptance=accept.sign_acceptance,passed=true))
                end
                temp=DataFrame(records)
                for field in [:gap_pp,:singleton_fraction,:high_site_fraction]
                    X=hcat([temp[temp.chain.==ch,field] for ch in 1:CHAINS]...);d=diagnostics(X)
                    push!(diagrows,(dataset=ds,quantity=String(field),rhat=d.rhat,ess=d.ess,mcse=d.mcse,swap_acceptance=NaN,sign_acceptance=NaN,passed=d.rhat<1.01 && d.ess>=400))
                    if field==:gap_pp;mcse=d.mcse;rhat=d.rhat;ess=d.ess;end
                end
            end
            @assert isapprox(sum(f),1;atol=1e-8) && isapprox(dot(f,1:N),M/L;atol=1e-8)
            pred=H*f;err=100 .* (pred[idx].-emp[idx]);drawtab=DataFrame(records);lo,hi=quantile(drawtab.gap_pp,[.025,.975]);sl,sh=quantile(drawtab.singleton_fraction,[.025,.975])
            previous=only(eachrow(old[(old.dataset.==ds).&(old.scope.=="LODO").&(old.model.==model),:]))
            push!(summaries,(dataset=ds,model=model,N=N,L=L,M=M,a=p.a,h=p.h,training_datasets=p.training_datasets,draws=nrow(drawtab),empirical_mean_gap_pp=empmean,predicted_mean_gap_pp=100mean(pred[idx]),lower_gap_pp=lo,upper_gap_pp=hi,mae_pp=mean(abs.(err)),maximum_error_pp=maximum(abs.(err)),signed_error_pp=mean(err),unconditioned_mae_pp=previous.mae_pp,mae_change_vs_unconditioned_pp=mean(abs.(err))-previous.mae_pp,empirical_singleton_percent=100empsingle,predicted_singleton_percent=100f[1],lower_singleton_percent=100sl,upper_singleton_percent=100sh,singleton_error_pp=100(f[1]-empsingle),mean_gap_mcse_pp=mcse,split_rhat=rhat,effective_sample_size=ess))
            for m in 0:N;push!(curves,(dataset=ds,model=model,m=m,removal_fraction=m/N,empirical_gap_pp=100emp[m+1],predicted_gap_pp=100pred[m+1]));end
            for k in 1:N;push!(freqrows,(dataset=ds,model=model,K=k,empirical_frequency=mean(K.==k),predicted_frequency=f[k]));end
            for r in records;push!(drawrows,merge((dataset=ds,model=model),r));end
            println(ds," ",model," fixed M: MAE ",round(mean(abs.(err));digits=3)," pp; Rhat ",round(rhat;digits=3));flush(stdout)
        end
        CSV.write(joinpath(DEST,"summary.csv"),DataFrame(summaries));CSV.write(joinpath(DEST,"curves.csv"),DataFrame(curves));CSV.write(joinpath(DEST,"draws.csv"),DataFrame(drawrows));CSV.write(joinpath(DEST,"site_diagnostics.csv"),DataFrame(diagrows));CSV.write(joinpath(DEST,"support_frequencies.csv"),DataFrame(freqrows))
    end
    CSV.write(joinpath(DEST,"input_versions.csv"),unique(ver));CSV.write(joinpath(DEST,"parameter_estimates_used.csv"),fits[(fits.scope.=="LODO").&in.(fits.model,Ref(["Common","Pair","Site"])),:])
    for file in ["metadata_relationships.csv","metadata_source_versions.csv"];cp(joinpath(SOURCE,file),joinpath(DEST,file);force=true);end
    tab=DataFrame(summaries);agg=NamedTuple[]
    for model in ["Common","Pair","Site"]
        v=tab[tab.model.==model,:];base=tab[tab.model.=="Common",:]
        push!(agg,(model=model,mean_mae_pp=mean(v.mae_pp),median_mae_pp=median(v.mae_pp),mean_absolute_singleton_error_pp=mean(abs.(v.singleton_error_pp)),mean_signed_gap_error_pp=mean(v.signed_error_pp),datasets_improved_vs_common=count(v.mae_pp.<base.mae_pp),unconditioned_mean_mae_pp=mean(v.unconditioned_mae_pp),datasets_improved_vs_unconditioned=count(v.mae_pp.<v.unconditioned_mae_pp)))
    end
    CSV.write(joinpath(DEST,"across_dataset_summary.csv"),DataFrame(agg))
    plots(tab,DataFrame(curves),DataFrame(drawrows))
    println("Done: ",DEST)
end
function plots(tab,curves,draws)
    set_theme!(Theme(font="Arial",fontsize=20,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=17,yticklabelsize=17)))
    colors=Dict("Common"=>"#777777","Pair"=>"#CC6677","Site"=>"#4477AA");offsets=Dict("Common"=>-.22,"Pair"=>0.,"Site"=>.22)
    fig=Figure(size=(1550,850))
    for (j,quantity,label) in [(1,"gap","Mean gap (percentage points)"),(2,"singleton","Singleton links (%)")]
        ax=Axis(fig[1,j];xlabel=label,yticks=(1:10,replace.(DATA,"_"=>" ")))
        for (i,ds) in enumerate(DATA),model in ["Common","Pair","Site"]
            r=only(eachrow(tab[(tab.dataset.==ds).&(tab.model.==model),:]));y=i+offsets[model]
            l=quantity=="gap" ? r.lower_gap_pp : r.lower_singleton_percent;u=quantity=="gap" ? r.upper_gap_pp : r.upper_singleton_percent
            mu=quantity=="gap" ? r.predicted_mean_gap_pp : r.predicted_singleton_percent
            obs=quantity=="gap" ? r.empirical_mean_gap_pp : r.empirical_singleton_percent
            lines!(ax,[l,u],[y,y];color=colors[model],linewidth=3);scatter!(ax,[mu],[y];color=colors[model],markersize=9)
            model=="Pair" && scatter!(ax,[obs],[i];color=:black,marker=:diamond,markersize=14)
        end
        ylims!(ax,.4,10.6);ax.yreversed=true
    end
    Legend(fig[2,1:2],vcat([LineElement(color=colors[m],linewidth=3) for m in ["Common","Pair","Site"]],[MarkerElement(color=:black,marker=:diamond)]),["Common / allocation null","Persistent pair","Shared site","Empirical"];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"01_conditional_allocations.png"),fig;px_per_unit=2)
    fig=Figure(size=(1450,760))
    for (j,label) in [(1,"Gap MAE with M fixed (pp)"),(2,"Change in MAE after fixing M (pp)")]
        ax=Axis(fig[1,j];xlabel=label,yticks=(1:10,replace.(DATA,"_"=>" ")))
        j==2 && vlines!(ax,[0];color=:gray,linestyle=:dash)
        for (i,ds) in enumerate(DATA),model in ["Common","Pair","Site"]
            r=only(eachrow(tab[(tab.dataset.==ds).&(tab.model.==model),:]))
            scatter!(ax,[j==1 ? r.mae_pp : r.mae_change_vs_unconditioned_pp],[i+offsets[model]];color=colors[model],markersize=11)
        end
        ylims!(ax,.4,10.6);ax.yreversed=true
    end
    Legend(fig[2,1:2],[MarkerElement(color=colors[m],marker=:circle) for m in ["Common","Pair","Site"]],["Common / allocation null","Persistent pair","Shared site"];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"02_prediction_errors.png"),fig;px_per_unit=2)
    fig=Figure(size=(1700,820))
    for (i,ds) in enumerate(DATA)
        ax=Axis(fig[cld(i,5),mod1(i,5)];title=replace(ds,"_"=>" "),xlabel="Sites removed (%)",ylabel=mod1(i,5)==1 ? "Gap (pp)" : "")
        for model in ["Common","Pair","Site"]
            v=curves[(curves.dataset.==ds).&(curves.model.==model).&(curves.removal_fraction.<=.8),:]
            lines!(ax,100v.removal_fraction,v.predicted_gap_pp;color=colors[model],linewidth=model=="Common" ? 3 : 2,linestyle=model=="Site" ? :dashdot : :solid)
            model=="Common" && lines!(ax,100v.removal_fraction,v.empirical_gap_pp;color=:black,linewidth=2.5,linestyle=:dash)
        end
        xlims!(ax,0,80);ylims!(ax,0,40)
    end
    Legend(fig[3,1:5],vcat([LineElement(color=colors[m],linewidth=2,linestyle=m=="Site" ? :dashdot : :solid) for m in ["Common","Pair","Site"]],[LineElement(color=:black,linestyle=:dash)]),["Common / allocation null","Persistent pair","Shared site","Empirical"];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"03_gap_curves.png"),fig;px_per_unit=2)
end
function report()
    tab=CSV.read(joinpath(DEST,"summary.csv"),DataFrame);curves=CSV.read(joinpath(DEST,"curves.csv"),DataFrame)
    freq=CSV.read(joinpath(DEST,"support_frequencies.csv"),DataFrame);diag=CSV.read(joinpath(DEST,"site_diagnostics.csv"),DataFrame)
    agg=CSV.read(joinpath(DEST,"across_dataset_summary.csv"),DataFrame)
    checks=NamedTuple[]
    for r in eachrow(tab)
        v=curves[(curves.dataset.==r.dataset).&(curves.model.==r.model),:]
        f=freq[(freq.dataset.==r.dataset).&(freq.model.==r.model),:]
        @assert v.m==collect(0:r.N) && all(v.removal_fraction.==v.m./r.N)
        @assert isapprox(sum(f.predicted_frequency),1;atol=1e-9)
        @assert isapprox(dot(f.K,f.predicted_frequency),r.M/r.L;atol=1e-8)
        @assert maximum(abs.(100gapmatrix(r.N)*f.predicted_frequency.-v.predicted_gap_pp))<1e-8
        e=(v.predicted_gap_pp.-v.empirical_gap_pp)[v.removal_fraction.<=.8]
        @assert isapprox(mean(abs.(e)),r.mae_pp;atol=1e-9)
        @assert isapprox(mean(e),r.signed_error_pp;atol=1e-9)
        push!(checks,(dataset=r.dataset,model=r.model,fixed_M=true,curve_from_support=true,metric_averaging=true))
    end
    CSV.write(joinpath(DEST,"result_checks.csv"),DataFrame(checks))
    finite=diag[isfinite.(diag.rhat),:];failed=finite[.!finite.passed,:]
    fmt(x)=@sprintf("%.3f",x)
    open(joinpath(DEST,"FINDINGS.md"),"w") do io
        println(io,"# Fixed-total-support comparison\n")
        println(io,"Existing entire-dataset-holdout parameters were kept fixed. Each model receives the target's opportunity matrix, known regional links, site count, and observed total support M. This evaluates support allocation conditional on M; it is not the earlier prediction task without target support information.\n")
        println(io,"| Model | Mean gap MAE (pp) | Median MAE (pp) | Mean absolute singleton error (pp) | Earlier MAE without fixed M (pp) |\n|---|---:|---:|---:|---:|")
        for r in eachrow(agg);println(io,"| $(r.model) | $(fmt(r.mean_mae_pp)) | $(fmt(r.median_mae_pp)) | $(fmt(r.mean_absolute_singleton_error_pp)) | $(fmt(r.unconditioned_mean_mae_pp)) |");end
        pair=only(eachrow(agg[agg.model.=="Pair",:]));common=only(eachrow(agg[agg.model.=="Common",:]))
        println(io,"\nPair improves gap MAE over the fixed-M Common/allocation null in $(pair.datasets_improved_vs_common) of ten datasets, reducing equal-dataset mean error by $(round(100(1-pair.mean_mae_pp/common.mean_mae_pp);digits=1))%. It also better describes singleton abundance on average. The Site model produces gaps and singleton fractions close to the allocation null; its transferred site heterogeneity supplies little improvement in this test. This concerns the existing fitted Site specification, not all possible site effects.\n")
        println(io,"| Dataset | Common MAE (pp) | Pair MAE (pp) | Site MAE (pp) |\n|---|---:|---:|---:|")
        for ds in DATA
            vals=[only(tab.mae_pp[(tab.dataset.==ds).&(tab.model.==m)]) for m in ["Common","Pair","Site"]]
            println(io,"| $(replace(ds,"_"=>" ")) | $(join(fmt.(vals)," | ")) |")
        end
        println(io,"\nNahuel is the exception to Pair's improvement over Common. Pair's errors are not uniformly signed across systems: some gaps remain too large and others too small. Its aggregate advantage is therefore not a universal fit or evidence of a uniquely identified mechanism. Most systems use the same transferred Pair parameters (a=0, h=2); Garraf PP2 uses a=−1, h=2. No target-specific parameter adjustment was made.\n")
        println(io,"Fixing M improves Pair relative to its earlier unconditioned curve in $(pair.datasets_improved_vs_unconditioned) datasets. The remaining advantage over Common is evidence that persistent pair differences provide a useful statistical description of allocation beyond opportunity counts and total support. Supplying M is additional target information, so the reduction relative to the earlier exercise should not be described as a like-for-like validation gain.\n")
        println(io,"## Verification and limits\n")
        println(io,"Pair and Common means use exact finite dynamic-programming marginals. Pair intervals use 2,500 independent exact conditional draws; Common intervals reuse 5,000 existing draws. Site intervals and means use four MCMC chains. A small enumerated system checks the samplers, and all 30 dataset/model combinations pass fixed-total, curve and metric consistency checks.\n")
        println(io,"Site diagnostic maximum split R-hat: $(fmt(maximum(finite.rhat))); minimum effective sample size: $(round(minimum(finite.ess);digits=1)). Maximum mean-gap Monte Carlo SE: $(fmt(maximum(tab.mean_gap_mcse_pp))) pp. Gottin PP and Montseny were extended after the initial diagnostic screen. These diagnostics do not prove convergence.\n")
        if nrow(failed)==0
            println(io,"All monitored Site summaries pass the prespecified screening thresholds (R-hat <1.01; ESS ≥400).\n")
        else
            println(io,"Remaining sampling flags: ",join(["$(r.dataset): $(r.quantity), R-hat=$(fmt(r.rhat)), ESS=$(round(r.ess;digits=1))" for r in eachrow(failed)],"; "),". Treat those Site summaries as provisional.\n")
        end
        println(io,"Intervals describe conditional allocation variation, not confidence intervals. Gap errors average integer removal levels through 80%, including zero, with equal dataset weighting for aggregate performance. Shared-site/collection dependence, interaction-derived occurrences, excluded uncertain restorations and known-link conditioning remain unchanged; see README and the copied metadata files. No new ecological mechanism is established.")
    end
    CSV.write(joinpath(DEST,"code_versions.csv"),DataFrame(file=["conditioned_M.jl"],sha256=[bytes2hex(sha256(read(@__FILE__)))]))
end
if length(ARGS)>=3 && ARGS[3]=="plot"
    plots(CSV.read(joinpath(DEST,"summary.csv"),DataFrame),CSV.read(joinpath(DEST,"curves.csv"),DataFrame),CSV.read(joinpath(DEST,"draws.csv"),DataFrame))
else
    main()
end
report()
end

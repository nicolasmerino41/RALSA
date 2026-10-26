module OpportunitySupportOrganisation
using CSV, DataFrames, Statistics, LinearAlgebra, Random, SHA, CairoMakie, Printf
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"../..")) : abspath(ARGS[1])
const SOURCE=joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory")
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(SOURCE,"opportunity_support_organisation")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const B=499
const R=300
const COLORS=["#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499","#661100"]
logz(x)=log(-expm1(x));lse(x)=maximum(x)+log(sum(exp.(x.-maximum(x))));logistic(x)=1/(1+exp(-x))
include(joinpath(ROOT,"scripts/beta_pair_conditioned_M/conditional_support.jl"))
mkpath(DEST)
function lorenz(x,q)
    v=sort(x);s=sum(v);s==0 && return NaN
    u=q*length(v);a=floor(Int,u)
    (sum(@view v[1:a])+(a<length(v) ? (u-a)*v[a+1] : 0.))/s
end
function spread(K,C)
    groups=[findall(==(c),C) for c in sort(unique(C))];groups=filter(g->length(g)>=3,groups)
    isempty(groups) && return NaN
    sum(sum(abs2,K[g]./C[g].-mean(K[g]./C[g])) for g in groups)/sum(length,groups)
end
function stats(K,C)
    isempty(K) && return (singleton=NaN,gini=NaN,top10=NaN,within_C_variance=NaN,mean_support=NaN,cv=NaN,p90=NaN)
    v=sort(K);L=length(v)
    (singleton=mean(K.==1),gini=2dot(1:L,v)/(L*sum(v))-(L+1)/L,top10=1-lorenz(v,.9),within_C_variance=spread(K,C),mean_support=mean(K),cv=std(K;corrected=false)/mean(K),p90=Float64(v[ceil(Int,.9L)]))
end
function enumeration_check()
    C=[2,3,4];M=6;N=4
    out=NamedTuple[]
    for h in [0.,2.]
        w=pairweights(C,0.,h);D,P=dp(w,M);f=exactfreq(M,D,P,N)
        values=[collect(k) for k in Iterators.product((1:c for c in C)...) if sum(k)==M]
        weights=exp.([sum(w[i][k[i]] for i in eachindex(C)) for k in values]);weights/=sum(weights);truth=zeros(N)
        for (k,p) in zip(values,weights),v in k;truth[v]+=p/length(C);end
        error=maximum(abs.(truth.-f));@assert error<1e-10
        push!(out,(check="Exact enumeration h=$h",error=error,passed=true))
    end
    CSV.write(joinpath(DEST,"algorithm_checks.csv"),DataFrame(out))
end
function main()
    enumeration_check()
    versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
    for r in eachrow(versions);@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256;end
    used=NamedTuple[]
    function load(path;kwargs...)
        push!(used,(path=relpath(path,ROOT),sha256=bytes2hex(sha256(read(path)))))
        CSV.read(path,DataFrame;kwargs...)
    end
    qc=load(joinpath(SOURCE,"quality_control.csv"));fits=load(joinpath(SOURCE,"fitted_parameters.csv"));old=load(joinpath(SOURCE,"conditioned_M/support_frequencies.csv"))
    sites=load(joinpath(ROOT,"outputs/statistical description of local support/site_universe.csv");types=Dict(:site=>String))
    empirical=NamedTuple[];freqs=NamedTuple[];concentration=NamedTuple[];drawrows=NamedTuple[];modelrows=NamedTuple[];conditional=NamedTuple[];rarefaction=NamedTuple[];rarefreq=NamedTuple[];checks=NamedTuple[]
    for (di,ds) in enumerate(DATA)
        g=qc[qc.dataset.==ds,:];K=Int.(g.K);C=Int.(g.C);N=only(unique(g.N));L=length(K);M=sum(K)
        ss=String.(sites.site[sites.dataset.==ds]);@assert length(ss)==N
        si=Dict(s=>i for (i,s) in enumerate(ss));li=Dict((r.consumer,r.resource)=>i for (i,r) in enumerate(eachrow(g)))
        occ=load(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv");types=String)
        ints=load(joinpath(ROOT,"data/processed",ds*"_interactions.csv");types=String)
        pres=Dict((v.trophic_level[1],v.species[1])=>Set(v.site) for v in groupby(occ,[:trophic_level,:species]))
        A=falses(L,N);O=falses(L,N)
        for r in eachrow(ints);A[li[(r.consumer,r.resource)],si[r.site]]=true;end
        for (i,r) in enumerate(eachrow(g)),s in intersect(pres[("consumer",r.consumer)],pres[("resource",r.resource)]);O[i,si[s]]=true;end
        @assert vec(sum(A;dims=2))==K && vec(sum(O;dims=2))==C && all(A.<=O) && all(1 .<=K .<=C .<=N)
        push!(checks,(dataset=ds,N=N,L=L,M=M,interaction_and_occurrence_verified=true))
        for (kind,values) in [("Support",K),("Opportunity",C)]
            z=stats(values,C);push!(empirical,merge((dataset=ds,quantity=kind,N=N,L=L,total=sum(values)),z))
            for k in 1:N;push!(freqs,(dataset=ds,quantity=kind,K=k,frequency=mean(values.==k),ccdf=mean(values.>=k)));end
            for q in 0.:.01:1.;push!(concentration,(dataset=ds,quantity=kind,link_fraction=q,mean_share=lorenz(values,q),lower_share=NaN,upper_share=NaN));end
        end
        observed=stats(K,C);classes=sort(unique(C));groups=[findall(==(c),C) for c in classes]
        for model in ["Common","Pair"]
            par=only(eachrow(fits[(fits.dataset.==ds).&(fits.scope.=="LODO").&(fits.model.==model),:]))
            @assert Set(split(par.training_datasets,";"))==Set(filter(!=(ds),DATA))
            D,P=dp(pairweights(C,par.a,par.h),M);f=exactfreq(M,D,P,N)
            ref=old[(old.dataset.==ds).&(old.model.==model),:];@assert maximum(abs.(f.-ref.predicted_frequency))<1e-8
            for k in 1:N;push!(freqs,(dataset=ds,quantity=model,K=k,frequency=f[k],ccdf=sum(f[k:end])));end
            rng=MersenneTwister(175000+di*10+(model=="Pair"));sim=NamedTuple[];lc=zeros(B,101);means=zeros(B,length(classes));sds=similar(means);samplefreq=zeros(N)
            for b in 1:B
                kk=drawK(rng,M,D,P);@assert sum(kk)==M && all(1 .<=kk .<=C)
                z=stats(kk,C);push!(sim,z);push!(drawrows,merge((dataset=ds,model=model,draw=b),z))
                for k in kk;samplefreq[k]+=1/(B*L);end
                for (j,q) in enumerate(0.:.01:1.);lc[b,j]=lorenz(kk,q);end
                for (j,ids) in enumerate(groups);means[b,j]=mean(kk[ids]);sds[b,j]=std(kk[ids];corrected=false);end
            end
            for field in [:singleton,:gini,:top10,:within_C_variance]
                vals=[z[field] for z in sim];vals=filter(isfinite,vals);obs=observed[field]
                isempty(vals) && continue
                lo,hi=quantile(vals,[.025,.975]);push!(modelrows,(dataset=ds,model=model,metric=String(field),observed=obs,predicted_mean=mean(vals),lower=lo,upper=hi,mean_mcse=std(vals)/sqrt(length(vals)),empirical_position=(count(<(obs),vals)+.5count(==(obs),vals))/length(vals),outside_interval=obs<lo || obs>hi,a=par.a,h=par.h,draws=B))
            end
            for (j,q) in enumerate(0.:.01:1.);lo,hi=quantile(lc[:,j],[.025,.975]);push!(concentration,(dataset=ds,quantity=model,link_fraction=q,mean_share=mean(lc[:,j]),lower_share=lo,upper_share=hi));end
            for (j,ids) in enumerate(groups)
                ml,mh=quantile(means[:,j],[.025,.975]);sl,sh=quantile(sds[:,j],[.025,.975])
                push!(conditional,(dataset=ds,model=model,C=classes[j],links=length(ids),display=length(ids)>=3,observed_mean_K=mean(K[ids]),observed_sd_K=std(K[ids];corrected=false),predicted_mean_K=mean(means[:,j]),lower_mean_K=ml,upper_mean_K=mh,predicted_sd_K=mean(sds[:,j]),lower_sd_K=sl,upper_sd_K=sh))
            end
            println(ds," ",model,": conditional draws complete");flush(stdout)
        end
        sizes=sort(unique(vcat(filter(n->n<=N,[4,8,14,20,32,64,128,256]),[N])))
        rng=MersenneTwister(275000+di);f8=zeros(8);c8=zeros(8);valid8=0
        for rep in 1:R
            order=randperm(rng,N);kk=zeros(Int,L);cc=zeros(Int,L);prev=0
            for n in sizes
                for j in order[prev+1:n];kk .+=A[:,j];cc .+=O[:,j];end
                keep=kk.>0;z=stats(kk[keep],cc[keep]);cz=stats(cc[keep],cc[keep]);retained=count(keep)
                push!(rarefaction,merge((dataset=ds,replicate=rep,sampled_sites=n,original_sites=N,original_links=L,retained_links=retained,lost_link_fraction=1-retained/L,zero_K_with_opportunity_fraction=count((kk.==0).&(cc.>0))/L,opportunity_gini=cz.gini,opportunity_singleton=cz.singleton,opportunity_top10=cz.top10),z))
                if n==8 && retained>0
                    valid8+=1
                    for k in 1:8;f8[k]+=mean(kk[keep].==k);c8[k]+=mean(cc[keep].==k);end
                end
                n==N && @assert kk==K && cc==C
                prev=n
            end
        end
        for k in 1:8;push!(rarefreq,(dataset=ds,K=k,support_frequency=f8[k]/valid8,opportunity_frequency=c8[k]/valid8,valid_replicates=valid8));end
        for (file,rows) in [("empirical_summary",empirical),("distributions",freqs),("concentration_curves",concentration),("allocation_draws",drawrows),("allocation_summary",modelrows),("conditional_C",conditional),("site_subsampling_draws",rarefaction),("eight_site_distributions",rarefreq),("input_checks",checks)]
            CSV.write(joinpath(DEST,file*".csv"),DataFrame(rows))
        end
        println(ds," site subsampling complete");flush(stdout)
    end
    rt=DataFrame(rarefaction);summary=NamedTuple[]
    for v in groupby(rt,[:dataset,:sampled_sites]),field in [:singleton,:gini,:top10,:mean_support,:opportunity_gini,:opportunity_singleton,:opportunity_top10,:lost_link_fraction,:retained_links,:zero_K_with_opportunity_fraction]
        values=filter(isfinite,v[!,field]);lo,hi=quantile(values,[.025,.975]);push!(summary,(dataset=v.dataset[1],sampled_sites=v.sampled_sites[1],original_sites=v.original_sites[1],metric=String(field),mean=mean(values),lower=lo,upper=hi,replicates=length(values)))
    end
    CSV.write(joinpath(DEST,"site_subsampling_summary.csv"),DataFrame(summary));CSV.write(joinpath(DEST,"input_versions.csv"),versions)
    helper=joinpath(ROOT,"scripts/beta_pair_conditioned_M/conditional_support.jl");push!(used,(path=relpath(helper,ROOT),sha256=bytes2hex(sha256(read(helper)))))
    CSV.write(joinpath(DEST,"source_versions.csv"),DataFrame(used));CSV.write(joinpath(DEST,"parameters_used.csv"),fits[(fits.scope.=="LODO").&in.(fits.model,Ref(["Common","Pair"])),:])
    for file in ["metadata_relationships.csv","metadata_source_versions.csv"];cp(joinpath(SOURCE,file),joinpath(DEST,file);force=true);end
end
include("figures_and_findings.jl")
if !(length(ARGS)>=3 && ARGS[3]=="plot");main();end
plots();report()
end

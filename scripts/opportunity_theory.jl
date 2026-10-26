module OpportunityTheory
using CSV, DataFrames, Statistics, Random, LinearAlgebra, SHA, CairoMakie
const ROOT = isempty(ARGS) ? normpath(joinpath(@__DIR__,"..")) : abspath(ARGS[1])
const DEST = length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(ROOT,"outputs","statistical description of local support","opportunity_theory")
const SRC=joinpath(ROOT,"outputs","statistical description of local support")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const COLORS=["#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499","#661100"]
const MODELS=["Common","Pair","Site"]
const RNG=MersenneTwister(731029)
const REPS=64
const PREDREPS=400
mkpath(DEST)
logistic(x)=1/(1+exp(-x))
lse(x)=maximum(x)+log(sum(exp.(x.-maximum(x))))
struct SystemData
    name::String
    N::Int
    O::Vector{Vector{Int}}
    A::BitMatrix
    C::Vector{Int}
    K::Vector{Int}
    keep::BitVector
end
function survival(N)
    H=ones(N+1,N)
    for m in 0:N
        e=1.
        for k in 1:m
            e*=(m-k+1)/(N-k+1);H[m+1,k]=1-e
        end
    end
    H
end
function gap(K,N,H=survival(N))
    f=[count(==(k),K)/length(K) for k in 1:N]
    H*f .- (1 .- collect(0:N)./N)
end
mainidx(N)=1:fld(4N,5)+1
function extremes(C,M)
    L=length(C);@assert L<=M<=sum(C)
    concentrated=ones(Int,L);left=M-L
    for i in sortperm(C;rev=true)
        extra=min(left,C[i]-1);concentrated[i]+=extra;left-=extra
    end
    balanced=ones(Int,L);left=M-L
    while left>0
        for i in sortperm(balanced)
            if balanced[i]<C[i];balanced[i]+=1;left-=1;end
            left==0 && break
        end
    end
    @assert sum(concentrated)==sum(balanced)==M
    concentrated,balanced
end
function checkbounds()
    cases=0
    for N in 2:5,L in 1:4,caps in Iterators.product(ntuple(_->1:N,L)...)
        C=collect(caps);H=survival(N)
        vectors=[collect(k) for k in Iterators.product((1:c for c in C)...)]
        for M in L:sum(C)
            lo,hi=extremes(C,M);glo=gap(lo,N,H);ghi=gap(hi,N,H)
            for K in vectors
                sum(K)==M || continue
                g=gap(K,N,H)
                @assert all(glo.-1e-12 .<= g) && all(g .<= ghi.+1e-12)
                for r in (.1,.37,.8)
                    @assert r-mean(r.^lo)-1e-12<=r-mean(r.^K)<=r-mean(r.^hi)+1e-12
                end
                cases+=1
            end
        end
    end
    CSV.write(joinpath(DEST,"algorithm_checks.csv"),DataFrame(check=["Exhaustive fixed-count and independent-removal bounds"],configurations=[cases],passed=[true]))
end
function loadinputs()
    ver=CSV.read(joinpath(SRC,"cooccurrence_prediction/input_versions.csv"),DataFrame)
    for z in eachrow(ver);@assert bytes2hex(sha256(read(joinpath(ROOT,z.path))))==z.sha256;end
    links=CSV.read(joinpath(SRC,"link_support.csv"),DataFrame)
    sites=CSV.read(joinpath(SRC,"site_universe.csv"),DataFrame;types=Dict(:site=>String))
    systems=SystemData[];qc=NamedTuple[]
    for ds in DATA
        l=links[links.dataset.==ds,:];ss=sites.site[sites.dataset.==ds];N=length(ss);siteid=Dict(s=>j for (j,s) in enumerate(ss))
        occ=CSV.read(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv"),DataFrame;types=String)
        pres=Dict((g.trophic_level[1],g.species[1])=>Set(g.site) for g in groupby(occ,[:trophic_level,:species]))
        O=[sort([siteid[s] for s in intersect(pres[("consumer",z.consumer)],pres[("resource",z.resource)])]) for z in eachrow(l)]
        A=falses(nrow(l),N);pairid=Dict((z.consumer,z.resource)=>i for (i,z) in enumerate(eachrow(l)))
        path=joinpath(ROOT,"data/processed",ds*"_interactions.csv")
        records=CSV.read(path,DataFrame;types=String)
        push!(ver,(relpath(path,ROOT),bytes2hex(sha256(read(path)))))
        for z in eachrow(records);A[pairid[(z.consumer,z.resource)],siteid[z.site]]=true;end
        K=vec(sum(A;dims=2));C=length.(O);@assert K==l.K
        for i in eachindex(K)
            valid=1<=K[i]<=C[i]<=N && all(j in O[i] for j in findall(A[i,:]))
            push!(qc,(dataset=ds,consumer=l.consumer[i],resource=l.resource[i],K=K[i],C=C[i],N=N,valid=valid))
        end
        keep=trues(N);keep[randperm(RNG,N)[1:max(1,round(Int,.2N))]].=false
        push!(systems,SystemData(ds,N,O,A,C,K,keep))
    end
    q=DataFrame(qc);CSV.write(joinpath(DEST,"quality_control.csv"),q);@assert all(q.valid)
    CSV.write(joinpath(DEST,"input_versions.csv"),ver)
    for f in ["metadata_relationships.csv","metadata_source_versions.csv"]
        cp(joinpath(SRC,"parameter_transfer",f),joinpath(DEST,f);force=true)
    end
    CSV.write(joinpath(DEST,"site_split.csv"),DataFrame([(dataset=s.name,site_index=j,training=s.keep[j]) for s in systems for j in 1:s.N]))
    systems
end
function bounds(systems)
    rows=NamedTuple[];summ=NamedTuple[];vectors=NamedTuple[]
    for s in systems
        N=s.N;L=length(s.K);M=sum(s.K);H=survival(N);emp=gap(s.K,N,H);idx=mainidx(N)
        for (label,C) in [("A",fill(N,L)),("B",s.C)]
            lo,hi=extremes(C,M);gl=gap(lo,N,H);gh=gap(hi,N,H);width=gh.-gl
            for i in 1:L;push!(vectors,(dataset=s.name,envelope=label,link_index=i,min_gap_support=lo[i],max_gap_support=hi[i]));end
            for m in 0:N
                w=width[m+1];z=w>1e-10 ? (emp[m+1]-gl[m+1])/w : NaN
                @assert gl[m+1]-1e-10<=emp[m+1]<=gh[m+1]+1e-10
                push!(rows,(dataset=s.name,envelope=label,m=m,removal_fraction=m/N,empirical=emp[m+1],minimum=gl[m+1],maximum=gh[m+1],width=w,Z=z))
            end
            il=mean(gl[idx]);ih=mean(gh[idx]);ie=mean(emp[idx])
            push!(summ,(dataset=s.name,envelope=label,N=N,L=L,M=M,integrated_empirical_pp=100ie,integrated_minimum_pp=100il,integrated_maximum_pp=100ih,mean_width_pp=100(ih-il),integrated_Z=ih-il>1e-10 ? (ie-il)/(ih-il) : NaN))
        end
    end
    CSV.write(joinpath(DEST,"envelope_curves.csv"),DataFrame(rows));CSV.write(joinpath(DEST,"envelope_summary.csv"),DataFrame(summ));CSV.write(joinpath(DEST,"extremal_support_vectors.csv"),DataFrame(vectors))
end
# Conditional Bernoulli sampling: force positivity through sequential conditional probabilities.
function drawrow!(rng,A,i,js,p;positive=true)
    isempty(js) && return
    suffix=zeros(length(js)+1)
    for t in length(js):-1:1;suffix[t]=suffix[t+1]+log1p(-p[t]);end
    found=!positive
    for t in eachindex(js)
        prob=found ? p[t] : p[t]/(-expm1(suffix[t]))
        hit=rand(rng)<clamp(prob,0,1);A[i,js[t]]=hit;found|=hit
    end
    @assert found
end
function simulate(rng,s,model,a,h)
    A=falses(length(s.O),s.N);sitep=logistic.(a .+ h.*rand(rng,[-1.,1.],s.N))
    for (i,js) in enumerate(s.O)
        p=model=="Site" ? sitep[js] : fill(logistic(a+(model=="Pair" ? h*rand(rng,[-1.,1.]) : 0)),length(js))
        drawrow!(rng,A,i,js,p)
    end
    A
end
function summaries(s,A,keep)
    K=vec(sum(A[:,keep];dims=2));T=vec(sum(A;dims=1));pred=zeros(s.N);v=zeros(s.N)
    for (i,js) in enumerate(s.O)
        c=count(j->keep[j],js);c==0 && continue
        p=K[i]/c
        for j in js
            if keep[j];pred[j]+=p;v[j]+=p*(1-p);end
        end
    end
    cluster=mean(((T[keep].-pred[keep]).^2)./(v[keep].+1))
    [log(mean(K)+.1),mean(K.==0),mean(K.==1),log(mean(K.^2)+.1),log1p(cluster)]
end
const GRID=vcat([(model="Common",a=a,h=0.) for a in -4.:.5:4.],[(model=model,a=a,h=h) for model in ["Pair","Site"] for a in -4.:.5:4. for h in [.5,1.,2.,3.]])
function library(systems)
    file=joinpath(DEST,"simulation_library.csv")
    signature=bytes2hex(sha256("positive-two-state-v2-conditional-site-residual|$(REPS)|$(GRID)|"*join([string(s.name,s.N,s.O,s.keep) for s in systems],"|")))
    stamp=joinpath(DEST,"simulation_library_signature.txt")
    if isfile(file) && isfile(stamp)
        @assert strip(read(stamp,String))==signature "Simulation cache is from another specification; use a new output directory"
        return CSV.read(file,DataFrame)
    end
    rows=NamedTuple[]
    for (d,s) in enumerate(systems)
        rng=MersenneTwister(10000+d)
        for (g,p) in enumerate(GRID),rep in 1:REPS
            A=simulate(rng,s,p.model,p.a,p.h)
            for (scope,keep) in [("full",trues(s.N)),("training",s.keep)]
                x=summaries(s,A,keep)
                push!(rows,(dataset=s.name,grid=g,replicate=rep,scope=scope,x1=x[1],x2=x[2],x3=x[3],x4=x[4],x5=x[5]))
            end
        end
        println("Simulation library: ",s.name);flush(stdout)
    end
    tab=DataFrame(rows);CSV.write(file,tab);write(stamp,signature*"\n");tab
end
function makeref(lib)
    refs=Dict()
    for g in groupby(lib,[:dataset,:grid,:scope])
        dims=g.scope[1]=="full" ? [1,3,4,5] : collect(1:5)
        X=Matrix(g[:,Symbol.("x".*string.(dims))]);mu=vec(mean(X;dims=1));V=cov(X)
        V=.8V+.2Diagonal(diag(V))+1e-5I
        refs[(g.dataset[1],g.grid[1],g.scope[1])]=(dims,mu,inv(V),logdet(V))
    end
    refs
end
function score(ref,x)
    dims,mu,Q,ld=ref;z=x[dims].-mu
    (dot(z,Q*z)+ld)/length(dims)
end
function fitgrid(refs,observations,scope,model)
    inds=findall(p->p.model==model,GRID)
    losses=[mean(score(refs[(ds,g,scope)],x) for (ds,x) in observations) for g in inds]
    inds[argmin(losses)],minimum(losses)
end
function recovery(systems,refs)
    rows=NamedTuple[]
    for d in [2,4,3]
        s=systems[d];rng=MersenneTwister(30000+d)
        for generating in MODELS,h in (generating=="Common" ? [0.] : [.5,1.,2.]),rep in 1:24
            A=simulate(rng,s,generating,-1.,h);x=summaries(s,A,trues(s.N))
            fits=[fitgrid(refs,[(s.name,x)],"full",model) for model in MODELS];winner=argmin(last.(fits))
            for (j,model) in enumerate(MODELS)
                g,loss=fits[j];p=GRID[g]
                push!(rows,(dataset=s.name,generating=generating,replicate=rep,true_a=-1.,true_h=h,fitted_model=model,a=p.a,h=p.h,score=loss,selected=j==winner,correct=MODELS[winner]==generating))
            end
        end
    end
    tab=DataFrame(rows);CSV.write(joinpath(DEST,"recovery.csv"),tab)
    CSV.write(joinpath(DEST,"recovery_summary.csv"),combine(groupby(tab[tab.selected,:],[:dataset,:generating,:true_h]),:correct=>mean=>:recovery_fraction,:a=>mean=>:mean_a,:h=>mean=>:mean_h))
    truth=tab[tab.generating.==tab.fitted_model,:]
    truth.a_error=truth.a.-truth.true_a;truth.h_error=truth.h.-truth.true_h
    CSV.write(joinpath(DEST,"parameter_recovery.csv"),combine(groupby(truth,[:dataset,:generating,:true_h]),:a_error=>mean=>:a_bias,:h_error=>mean=>:h_bias,:a_error=>(x->sqrt(mean(x.^2)))=>:a_rmse,:h_error=>(x->sqrt(mean(x.^2)))=>:h_rmse))
end
function meanactivity(s,model,a,h,signs)
    if model!="Site"
        qs=model=="Common" ? [logistic(a)] : logistic.([a-h,a+h])
        return sum(mean(c*q/(-expm1(c*log1p(-q))) for q in qs) for c in s.C)
    end
    qlo,qhi=logistic(a-h),logistic(a+h)
    total=0.
    for c in unique(s.C)
        weight=2.0^(-c);value=0.
        for t in 0:c
            value+=weight*(t*qhi+(c-t)*qlo)/(-expm1(t*log1p(-qhi)+(c-t)*log1p(-qlo)))
            t<c && (weight*=(c-t)/(t+1))
        end
        total+=count(==(c),s.C)*value
    end
    total
end
function scenarios(systems)
    rows=NamedTuple[];curves=NamedTuple[]
    for (d,s) in enumerate(systems)
        rng=MersenneTwister(40000+d);signs=rand(rng,[-1.,1.],s.N,128);target=sum(s.K);H=survival(s.N);idx=mainidx(s.N)
        for model in MODELS,h in (model=="Common" ? [0.] : [1.,2.,3.])
            low=-14.;high=14.
            for _ in 1:38
                mid=(low+high)/2
                if meanactivity(s,model,mid,h,signs)>target;high=mid;else;low=mid;end
            end
            a=(low+high)/2;G=zeros(s.N+1,96)
            for rep in 1:96
                A=simulate(rng,s,model,a,h);K=vec(sum(A;dims=2));G[:,rep]=gap(K,s.N,H)
                lo,hi=extremes(s.C,sum(K));gl=mean(gap(lo,s.N,H)[idx]);gh=mean(gap(hi,s.N,H)[idx]);ig=mean(G[idx,rep])
                push!(rows,(dataset=s.name,model=model,h=h,a=a,replicate=rep,target_M=target,expected_M=meanactivity(s,model,a,h,signs),realised_M=sum(K),singleton=mean(K.==1),tail_ge5=mean(K.>=5),q90=quantile(K,.9),integrated_gap_pp=100ig,Z_B=gh-gl>1e-10 ? (ig-gl)/(gh-gl) : NaN))
            end
            for m in 0:s.N;push!(curves,(dataset=s.name,model=model,h=h,m=m,removal_fraction=m/s.N,mean_gap=mean(G[m+1,:]),lower=quantile(G[m+1,:],.025),upper=quantile(G[m+1,:],.975)));end
        end
        println("Controlled scenarios: ",s.name);flush(stdout)
    end
    CSV.write(joinpath(DEST,"controlled_scenarios.csv"),DataFrame(rows));CSV.write(joinpath(DEST,"scenario_curves.csv"),DataFrame(curves))
end
function rowloglik(s,i,p,keep)
    js=s.O[i];lt=0.;lzero=0.;ltestzero=0.;hit=false
    for (t,j) in enumerate(js)
        q=p[t];lzero+=log1p(-q)
        if keep[j];hit|=s.A[i,j];lt+=s.A[i,j] ? log(q) : log1p(-q)
        else;ltestzero+=log1p(-q);end
    end
    lt-log(-expm1(lzero))+(hit ? 0. : log(-expm1(ltestzero)))
end
function conditionaldraw(rng,s,pairp,sitep,model)
    A=copy(s.A);A[:,.!s.keep].=false
    for (i,js) in enumerate(s.O)
        test=filter(j->!s.keep[j],js);positive=!any(s.A[i,s.keep])
        p=model=="Site" ? sitep[test] : fill(pairp[i],length(test))
        drawrow!(rng,A,i,test,p;positive=positive)
    end
    A
end
function withinpredict(rng,s,p)
    matrices=BitMatrix[];chainrows=NamedTuple[]
    if p.model!="Site"
        qs=p.model=="Pair" ? logistic.([p.a-p.h,p.a+p.h]) : [logistic(p.a)]
        weights=[exp.(let x=[rowloglik(s,i,fill(q,s.C[i]),s.keep) for q in qs];x.-lse(x);end) for i in eachindex(s.C)]
        for _ in 1:PREDREPS
            pairp=[qs[length(qs)==1 || rand(rng)<w[1] ? 1 : 2] for w in weights]
            push!(matrices,conditionaldraw(rng,s,pairp,Float64[],p.model))
        end
        return matrices,chainrows
    end
    neighbors=[[i for i in eachindex(s.O) if j in s.O[i]] for j in 1:s.N]
    for chain in 1:2
        z=fill(chain==1 ? -1. : 1.,s.N);q=logistic.(p.a .+ p.h.*z)
        ll=[rowloglik(s,i,q[s.O[i]],s.keep) for i in eachindex(s.O)];accepted=0;trace=Float64[]
        @assert all(isfinite,ll)
        for sweep in 1:1800
            for j in randperm(rng,s.N)
                old=q[j];q[j]=logistic(p.a-p.h*z[j]);ii=neighbors[j]
                new=[rowloglik(s,i,q[s.O[i]],s.keep) for i in ii]
                if log(rand(rng))<sum(new)-sum(ll[ii]);z[j]*=-1;ll[ii]=new;accepted+=1;else;q[j]=old;end
            end
            for block in [collect(1:s.N),randperm(rng,s.N)[1:max(1,fld(s.N,2))]]
                proposal=copy(q);proposal[block]=logistic.(p.a .- p.h.*z[block])
                new=[rowloglik(s,i,proposal[s.O[i]],s.keep) for i in eachindex(s.O)]
                if log(rand(rng))<sum(new)-sum(ll);z[block].*=-1;q=proposal;ll=new;end
            end
            if sweep>800 && sweep%5==0
                push!(matrices,conditionaldraw(rng,s,Float64[],q,"Site"));push!(trace,mean(q))
            end
        end
        push!(chainrows,(dataset=s.name,chain=chain,acceptance=accepted/(1800s.N),mean_site_probability=mean(trace),sd_site_probability=std(trace)))
    end
    matrices,chainrows
end
function evaluate(systems,refs)
    fits=NamedTuple[];metrics=NamedTuple[];curves=NamedTuple[];diagnostics=NamedTuple[];conditional=NamedTuple[];chains=NamedTuple[];frequency=NamedTuple[]
    baseline=CSV.read(joinpath(SRC,"parameter_transfer/curves.csv"),DataFrame)
    fullobs=[(s.name,summaries(s,s.A,trues(s.N))) for s in systems]
    for (d,s) in enumerate(systems),scope in ["LODO","Within"]
        H=survival(s.N);idx=mainidx(s.N);emp=gap(s.K,s.N,H);meanemp=mean(emp[idx]);rng=MersenneTwister(50000+100d+(scope=="Within"))
        obs=scope=="LODO" ? filter(x->x[1]!=s.name,fullobs) : [(s.name,summaries(s,s.A,s.keep))]
        for model in MODELS
            g,loss=fitgrid(refs,obs,scope=="LODO" ? "full" : "training",model);p=GRID[g]
            push!(fits,(dataset=s.name,scope=scope,model=model,a=p.a,h=p.h,score=loss,boundary=abs(p.a)==4 || p.h==3,training_datasets=join(first.(obs),";")))
            mats,cc=scope=="LODO" ? ([simulate(rng,s,model,p.a,p.h) for _ in 1:PREDREPS],NamedTuple[]) : withinpredict(rng,s,p)
            append!(chains,cc);Ks=hcat([vec(sum(A;dims=2)) for A in mats]...);G=hcat([gap(Ks[:,b],s.N,H) for b in axes(Ks,2)]...);pred=vec(mean(G;dims=2));err=100(pred[idx].-emp[idx])
            cells=[(i,j) for (i,js) in enumerate(s.O) for j in js if scope=="LODO" || !s.keep[j]]
            counts=[sum(A[i,j] for A in mats) for (i,j) in cells];phat=(counts.+.5)./(length(mats)+1);y=[s.A[i,j] for (i,j) in cells]
            brier=mean((phat.-y).^2);logscore=-mean(y.*log.(phat).+(1 .- y).*log1p.(-phat))
                drawmeans=vec(mean(G[idx,:];dims=1));batchmeans=[mean(drawmeans[j:j+19]) for j in 1:20:length(drawmeans)-19]
                push!(metrics,(dataset=s.name,scope=scope,model=model,mae_pp=mean(abs.(err)),maximum_error_pp=maximum(abs.(err)),signed_error_pp=mean(err),relative_error_percent=100mean(abs.(err))/(100meanemp),brier=brier,marginal_log_loss=logscore,evaluated_cells=length(cells),prediction_draws=length(mats),mean_gap_mcse_pp=100std(batchmeans)/sqrt(length(batchmeans))))
            for m in 0:s.N;push!(curves,(dataset=s.name,scope=scope,model=model,m=m,removal_fraction=m/s.N,empirical=emp[m+1],predicted=pred[m+1],lower=quantile(G[m+1,:],.025),upper=quantile(G[m+1,:],.975)));end
            for k in 1:s.N;push!(frequency,(dataset=s.name,scope=scope,model=model,K=k,observed=mean(s.K.==k),predicted=mean(Ks.==k),observed_ccdf=mean(s.K.>=k),predicted_ccdf=mean(Ks.>=k)));end
            for c in sort(unique(s.C))
                ii=findall(==(c),s.C)
                push!(conditional,(dataset=s.name,scope=scope,model=model,C=c,links=length(ii),observed_mean_K=mean(s.K[ii]),predicted_mean_K=mean(Ks[ii,:]),observed_singleton=mean(s.K[ii].==1),predicted_singleton=mean(Ks[ii,:].==1)))
            end
            observedstats=[mean(s.K.==1),mean(s.K.>=5),quantile(s.K,.9),sum(s.K),summaries(s,s.A,trues(s.N))[5]]
            simstats=hcat([[mean(Ks[:,b].==1),mean(Ks[:,b].>=5),quantile(Ks[:,b],.9),sum(Ks[:,b]),summaries(s,mats[b],trues(s.N))[5]] for b in eachindex(mats)]...)
            for (j,stat) in enumerate(["singleton","tail_ge5","q90","M","log_site_residual_clustering"])
                vals=simstats[j,:];push!(diagnostics,(dataset=s.name,scope=scope,model=model,statistic=stat,observed=observedstats[j],predicted=mean(vals),lower=quantile(vals,.025),upper=quantile(vals,.975),covered=quantile(vals,.025)<=observedstats[j]<=quantile(vals,.975)))
            end
            # Positions use each simulated matrix's own M, not target activity.
            zz=Float64[]
            for b in 1:4:length(mats)
                lo,hi=extremes(s.C,sum(Ks[:,b]));gl=mean(gap(lo,s.N,H)[idx]);gh=mean(gap(hi,s.N,H)[idx])
                gh-gl>1e-10 && push!(zz,(mean(G[idx,b])-gl)/(gh-gl))
            end
            lo,hi=extremes(s.C,sum(s.K));gl=mean(gap(lo,s.N,H)[idx]);gh=mean(gap(hi,s.N,H)[idx]);oz=(meanemp-gl)/(gh-gl)
            push!(diagnostics,(dataset=s.name,scope=scope,model=model,statistic="Z_B",observed=oz,predicted=mean(zz),lower=quantile(zz,.025),upper=quantile(zz,.975),covered=quantile(zz,.025)<=oz<=quantile(zz,.975)))
        end
        println("Prediction: ",s.name," ",scope);flush(stdout)
    end
    for (file,rows) in [("fitted_parameters",fits),("prediction_metrics",metrics),("prediction_curves",curves),("prediction_diagnostics",diagnostics),("conditional_support",conditional),("site_sampler_checks",chains),("support_predictions",frequency)]
        CSV.write(joinpath(DEST,file*".csv"),DataFrame(rows))
    end
    t=DataFrame(metrics);CSV.write(joinpath(DEST,"performance_summary.csv"),combine(groupby(t,[:scope,:model]),:mae_pp=>mean=>:mean_mae_pp,:mae_pp=>median=>:median_mae_pp,:relative_error_percent=>mean=>:mean_relative_error_percent,:brier=>mean=>:mean_brier,:marginal_log_loss=>mean=>:mean_log_loss))
    comparisons=NamedTuple[]
    for s in systems
        prior=baseline[(baseline.dataset.==s.name).&(baseline.removal_fraction.<=.8),:]
        common=only(t.mae_pp[(t.dataset.==s.name).&(t.scope.=="LODO").&(t.model.=="Common")])
        for z in eachrow(t[(t.dataset.==s.name).&(t.scope.=="LODO"),:])
            push!(comparisons,(dataset=s.name,model=z.model,mae_pp=z.mae_pp,common_rule_mae_pp=common,transfer_baseline_mae_pp=mean(abs.(prior.transferred_error_pp)),own_fit_mae_pp=mean(abs.(prior.own_fit_error_pp)),change_vs_transfer_pp=z.mae_pp-mean(abs.(prior.transferred_error_pp))))
        end
    end
    CSV.write(joinpath(DEST,"baseline_comparison.csv"),DataFrame(comparisons))
end
function plots()
    set_theme!(Theme(font="Arial",fontsize=18,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=13,yticklabelsize=13)))
    env=CSV.read(joinpath(DEST,"envelope_curves.csv"),DataFrame);summ=CSV.read(joinpath(DEST,"envelope_summary.csv"),DataFrame)
    f=Figure(size=(1700,760))
    for (i,ds) in enumerate(DATA)
        ax=Axis(f[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Sites removed (%)",ylabel=mod1(i,5)==1 ? "Gap (pp)" : "")
        for (label,col) in [("A",(:gray,.18)),("B",(COLORS[i],.28))]
            v=env[(env.dataset.==ds).&(env.envelope.==label).&(env.removal_fraction.<=.8),:]
            band!(ax,100v.removal_fraction,100v.minimum,100v.maximum;color=col)
            label=="B" && lines!(ax,100v.removal_fraction,100v.empirical;color=COLORS[i],linewidth=3)
        end
        xlims!(ax,0,80);ylims!(ax,0,60)
    end
    Legend(f[3,1:5],[PolyElement(color=(:gray,.3)),PolyElement(color=(:blue,.3)),LineElement(color=:black)],["N, L, M fixed","Also constrain opportunities","Empirical"];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"01_envelopes.png"),f;px_per_unit=2)
    f=Figure(size=(1150,650))
    for (j,field,label) in [(1,:mean_width_pp,"Mean envelope width (pp)"),(2,:integrated_Z,"Empirical normalised position")]
        ax=Axis(f[1,j],yticks=(1:10,replace.(DATA,"_"=>" ")),yreversed=true,xlabel=label)
        for (i,ds) in enumerate(DATA)
            v=summ[summ.dataset.==ds,:];a=only(v[v.envelope.=="A",field]);b=only(v[v.envelope.=="B",field])
            lines!(ax,[a,b],[i,i];color=COLORS[i]);scatter!(ax,[a],[i];color=:gray,marker=:diamond);scatter!(ax,[b],[i];color=COLORS[i],markersize=12)
        end
        xlims!(ax,0,j==1 ? 30 : 1)
    end
    Legend(f[2,1:2],[MarkerElement(color=:gray,marker=:diamond),MarkerElement(color=:black,marker=:circle)],["N, L, M fixed","Also constrain opportunities"];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"02_widths_positions.png"),f;px_per_unit=2)
    scen=CSV.read(joinpath(DEST,"controlled_scenarios.csv"),DataFrame);f=Figure(size=(1250,650));styles=Dict("Common"=>:solid,"Pair"=>:dash,"Site"=>:dot)
    for (j,field,label) in [(1,:singleton,"Singleton fraction"),(2,:integrated_gap_pp,"Mean gap (pp)")]
        ax=Axis(f[1,j],xlabel="Heterogeneity h (log odds)",ylabel=label)
        for (i,ds) in enumerate(DATA),model in ["Pair","Site"]
            v=scen[(scen.dataset.==ds).&(scen.model.==model),:];c=mean(scen[(scen.dataset.==ds).&(scen.model.=="Common"),field])
            yy=[c;[mean(v[v.h.==h,field]) for h in [1.,2.,3.]]];lines!(ax,[0,1,2,3],yy;color=COLORS[i],linestyle=styles[model],linewidth=2)
        end
    end
    Legend(f[2,1:2],[LineElement(color=:black,linestyle=:dash),LineElement(color=:black,linestyle=:dot)],["Persistent pair differences","Shared site differences"];orientation=:horizontal,framevisible=false)
    Legend(f[3,1:2],[LineElement(color=c) for c in COLORS],replace.(DATA,"_"=>" ");orientation=:horizontal,nbanks=2,framevisible=false)
    save(joinpath(DEST,"03_controlled_scenarios.png"),f;px_per_unit=2)
    pred=CSV.read(joinpath(DEST,"prediction_curves.csv"),DataFrame);baseline=CSV.read(joinpath(SRC,"parameter_transfer/curves.csv"),DataFrame)
    for scope in ["LODO","Within"]
        f=Figure(size=(1700,790))
        for (i,ds) in enumerate(DATA)
            ax=Axis(f[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Sites removed (%)",ylabel=mod1(i,5)==1 ? "Gap (pp)" : "")
            for model in MODELS
                v=pred[(pred.dataset.==ds).&(pred.scope.==scope).&(pred.model.==model).&(pred.removal_fraction.<=.8),:]
                lines!(ax,100v.removal_fraction,100v.predicted;color=:black,linestyle=styles[model],linewidth=2)
                model=="Common" && lines!(ax,100v.removal_fraction,100v.empirical;color=COLORS[i],linewidth=3)
            end
            if scope=="LODO"
                v=baseline[(baseline.dataset.==ds).&(baseline.removal_fraction.<=.8),:];lines!(ax,100v.removal_fraction,100v.transferred_gap;color=:gray55,linestyle=:dashdot)
                lines!(ax,100v.removal_fraction,100v.own_fit_gap;color=:gray75,linestyle=:solid)
            end
            xlims!(ax,0,80);ylims!(ax,0,60)
        end
        elements=[LineElement(color=:black,linestyle=styles[m]) for m in MODELS];labels=copy(MODELS)
        push!(elements,LineElement(color=COLORS[1],linewidth=3));push!(labels,"Empirical (dataset colour)")
        if scope=="LODO"
            push!(elements,LineElement(color=:gray55,linestyle=:dashdot));push!(labels,"Parameter transfer")
            push!(elements,LineElement(color=:gray75));push!(labels,"Own fit (in-sample)")
        end
        Legend(f[3,1:5],elements,labels;orientation=:horizontal,framevisible=false)
        save(joinpath(DEST,scope=="LODO" ? "04_transfer_gaps.png" : "05_within_system_gaps.png"),f;px_per_unit=2)
    end
    diag=CSV.read(joinpath(DEST,"prediction_diagnostics.csv"),DataFrame);f=Figure(size=(1600,850))
    for (j,stat) in enumerate(["singleton","M","log_site_residual_clustering"]), (i,model) in enumerate(MODELS)
        label=stat=="singleton" ? "Singleton fraction" : stat=="M" ? "Total activity" : "Site clustering (log)"
        ax=Axis(f[i,j],title=model*": "*label,xlabel="Observed",ylabel="Held-out prediction")
        v=diag[(diag.scope.=="LODO").&(diag.model.==model).&(diag.statistic.==stat),:];lim=maximum(vcat(v.observed,v.upper))*1.05
        lines!(ax,[0,lim],[0,lim];color=:gray70,linestyle=:dash)
        for (d,ds) in enumerate(DATA)
            z=only(eachrow(v[v.dataset.==ds,:]));rangebars!(ax,[z.observed],[z.lower],[z.upper];color=COLORS[d]);scatter!(ax,[z.observed],[z.predicted];color=COLORS[d],markersize=10)
        end
    end
    Legend(f[4,1:3],[MarkerElement(color=c,marker=:circle) for c in COLORS],replace.(DATA,"_"=>" ");orientation=:horizontal,nbanks=2,framevisible=false)
    save(joinpath(DEST,"06_diagnostic_predictions.png"),f;px_per_unit=2)
    support=CSV.read(joinpath(DEST,"support_predictions.csv"),DataFrame);f=Figure(size=(1700,790))
    for (i,ds) in enumerate(DATA)
        ticks=filter(k->k<=maximum(support.K[support.dataset.==ds]),[1,2,5,10,20,50,100,200])
        ax=Axis(f[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Support K",ylabel=mod1(i,5)==1 ? "Proportion with support ≥ K" : "",xscale=log10,xticks=(ticks,string.(ticks)))
        for model in MODELS
            v=support[(support.dataset.==ds).&(support.scope.=="LODO").&(support.model.==model),:]
            lines!(ax,v.K,v.predicted_ccdf;color=:black,linestyle=styles[model],linewidth=2)
            model=="Common" && lines!(ax,v.K,v.observed_ccdf;color=COLORS[i],linewidth=3)
        end
        ylims!(ax,0,1)
    end
    Legend(f[3,1:5],vcat([LineElement(color=:black,linestyle=styles[m]) for m in MODELS],[LineElement(color=COLORS[1],linewidth=3)]),vcat(MODELS,["Empirical (dataset colour)"]);orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"07_support_predictions.png"),f;px_per_unit=2)
end
function main()
    stage=get(ENV,"OPPORTUNITY_STAGE","all")
    stage=="none" && return
    stage=="plots" && return plots()
    checkbounds();println("Exhaustive bounds passed");flush(stdout)
    systems=loadinputs();bounds(systems)
    lib=library(systems);refs=makeref(lib)
    if stage=="recovery";recovery(systems,refs);return;end
    if stage=="scenarios";scenarios(systems);return;end
    if stage=="all";recovery(systems,refs);scenarios(systems);end
    evaluate(systems,refs);plots()
    println("Complete: ",DEST)
end
main()
end

module NeutralAllocation
using CSV, DataFrames, Random, Statistics, SHA, CairoMakie
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"..")) : abspath(ARGS[1])
const SOURCE=joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory")
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(SOURCE,"neutral_allocation")
const B=5000
const COLORS=["#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499","#661100"]
mkpath(DEST)
function probabilities(c,p)
    logw=[sum(log((c-j+1)/j) for j in 1:k)+k*log(p)+(c-k)*log1p(-p) for k in 1:c]
    w=exp.(logw.-maximum(logw));w/sum(w)
end
function table(C,M)
    low=1e-8;high=1-1e-8
    for _ in 1:55
        p=(low+high)/2
        if sum(c*p/(-expm1(c*log1p(-p))) for c in C)>M;high=p;else;low=p;end
    end
    P=[probabilities(c,(low+high)/2) for c in C]
    D=zeros(length(C)+1,M+1);D[1,1]=1
    for i in eachindex(C)
        for total in i:min(sum(@view C[1:i]),M),k in 1:min(C[i],total-i+1)
            D[i+1,total+1]+=P[i][k]*D[i,total-k+1]
        end
        D[i+1,:]./=maximum(@view D[i+1,:])
    end
    @assert D[end,M+1]>0
    D,P
end
function draw(rng,C,M,D,P)
    K=zeros(Int,length(C));left=M
    for i in length(C):-1:1
        upto=min(C[i],left-i+1);w=[P[i][k]*D[i,left-k+1] for k in 1:upto]
        mass=sum(w);@assert mass>0
        u=rand(rng)*mass;v=0.;chosen=upto
        for k in 1:upto;v+=w[k];if u<=v;chosen=k;break;end;end
        K[i]=chosen;left-=chosen
    end
    @assert left==0 && sum(K)==M && all(1 .<= K .<= C)
    K
end
function check()
    tested=0;err=0.
    for C in [[2,3,4],[1,3,2],[4,4,4]],M in length(C):sum(C)
        D,P=table(C,M)
        configurations=[collect(k) for k in Iterators.product((1:c for c in C)...) if sum(k)==M]
        exact=[prod(binomial(c,k) for (c,k) in zip(C,K)) for K in configurations];exact=exact/sum(exact)
        for (K,truth) in zip(configurations,exact)
            left=M;prob=1.
            for i in length(C):-1:1
                w=[P[i][k]*D[i,left-k+1] for k in 1:min(C[i],left-i+1)]
                prob*=w[K[i]]/sum(w);left-=K[i]
            end
            err=max(err,abs(prob-truth));tested+=1
        end
    end
    @assert err<1e-12
    CSV.write(joinpath(DEST,"sampler_check.csv"),DataFrame(configurations=[tested],maximum_probability_error=[err],passed=[true]))
end
function main()
    check()
    versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
    for z in eachrow(versions);@assert bytes2hex(sha256(read(joinpath(ROOT,z.path))))==z.sha256;end
    qc=CSV.read(joinpath(SOURCE,"quality_control.csv"),DataFrame);@assert all(qc.valid)
    envelopes=CSV.read(joinpath(SOURCE,"envelope_summary.csv"),DataFrame)
    versions.path=String.(versions.path);versions.sha256=String.(versions.sha256)
    for f in ["quality_control.csv","envelope_summary.csv"]
        path=joinpath(SOURCE,f);push!(versions,(relpath(path,ROOT),bytes2hex(sha256(read(path)))))
    end
    CSV.write(joinpath(DEST,"input_versions.csv"),versions)
    summary=NamedTuple[];draws=NamedTuple[]
    for (d,g) in enumerate(groupby(qc,:dataset))
        ds=g.dataset[1];N=g.N[1];C=g.C;K=g.K;M=sum(K);L=length(K);idx=0:fld(4N,5)
        coefficients=zeros(N)
        for k in 1:N
            coefficients[k]=mean(m/N-(m<k ? 0. : prod((m-j)/(N-j) for j in 0:k-1)) for m in idx)
        end
        e=only(eachrow(envelopes[(envelopes.dataset.==ds).&(envelopes.envelope.=="B"),:]));emp=100mean(coefficients[K])
        @assert isapprox(emp,e.integrated_empirical_pp;atol=1e-10)
        lo=e.integrated_minimum_pp;hi=e.integrated_maximum_pp
        D,P=table(C,M);rng=MersenneTwister(731113+d);values=zeros(B);singletons=zeros(B)
        for b in 1:B
            sim=draw(rng,C,M,D,P);values[b]=100mean(coefficients[sim]);singletons[b]=mean(sim.==1)
            push!(draws,(dataset=ds,replicate=b,gap_pp=values[b],Z=(values[b]-lo)/(hi-lo),singleton_fraction=singletons[b]))
        end
        lower,upper=quantile(values,[.025,.975]);pbelow=(1+count(x->x<=emp,values))/(B+1);pabove=(1+count(x->x>=emp,values))/(B+1)
        push!(summary,(dataset=ds,N=N,L=L,M=M,replicates=B,empirical_gap_pp=emp,null_mean_gap_pp=mean(values),null_lower_gap_pp=lower,null_upper_gap_pp=upper,
            difference_pp=emp-mean(values),empirical_Z=(emp-lo)/(hi-lo),null_mean_Z=(mean(values)-lo)/(hi-lo),null_lower_Z=(lower-lo)/(hi-lo),null_upper_Z=(upper-lo)/(hi-lo),
            empirical_percentile=100mean(values.<=emp),two_sided_MC_p=min(1.,2min(pbelow,pabove)),position=emp<lower ? "below" : emp>upper ? "above" : "inside",
            empirical_singleton=mean(K.==1),null_singleton=mean(singletons),null_mean_mcse_pp=std(values)/sqrt(B)))
        println(ds,": empirical=",round(emp;digits=2),", null=",round(mean(values);digits=2)," [",round(lower;digits=2),", ",round(upper;digits=2),"]");flush(stdout)
    end
    tab=DataFrame(summary);CSV.write(joinpath(DEST,"summary.csv"),tab);CSV.write(joinpath(DEST,"draws.csv"),DataFrame(draws))
    for f in ["metadata_relationships.csv","metadata_source_versions.csv"];cp(joinpath(SOURCE,f),joinpath(DEST,f);force=true);end
    set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false)))
    fig=Figure(size=(1250,670))
    for (j,field,label) in [(1,"Z","Position within opportunity-constrained envelope"),(2,"gap_pp","Mean gap (percentage points)")]
        ax=Axis(fig[1,j],yticks=(1:nrow(tab),replace.(tab.dataset,"_"=>" ")),yreversed=true,xlabel=label)
        for (i,z) in enumerate(eachrow(tab))
            l=z[Symbol("null_lower_"*field)];u=z[Symbol("null_upper_"*field)];mu=z[Symbol("null_mean_"*field)];obs=z[Symbol("empirical_"*field)]
            lines!(ax,[l,u],[i,i];color=COLORS[i],linewidth=7);scatter!(ax,[mu],[i];color=:white,strokecolor=COLORS[i],strokewidth=2,markersize=10)
            scatter!(ax,[obs],[i];color=:black,marker=:diamond,markersize=13)
        end
        j==1 && xlims!(ax,0,1)
    end
    Legend(fig[2,1:2],[LineElement(color=:gray,linewidth=7),MarkerElement(color=:white,strokecolor=:gray,strokewidth=2,marker=:circle),MarkerElement(color=:black,marker=:diamond)],
        ["Random allocation: central 95%","Random mean","Empirical"];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"neutral_comparison.png"),fig;px_per_unit=2)
end
main()
end

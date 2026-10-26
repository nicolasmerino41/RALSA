module SiteAccumulation
using CSV,DataFrames,Statistics,Random,LinearAlgebra,SHA,CairoMakie,Printf
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"../..")) : abspath(ARGS[1])
const SOURCE=joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory")
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(SOURCE,"site_accumulation")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const B=300
const MODELS=["Empirical","Common","Pair"]
const COLORS=Dict("Empirical"=>"#111111","Common"=>"#999999","Pair"=>"#CC6677")
logistic(x)=1/(1+exp(-x));logz(x)=log(-expm1(x));lse(x)=maximum(x)+log(sum(exp.(x.-maximum(x))))
include(joinpath(ROOT,"scripts/beta_pair_conditioned_M/conditional_support.jl"))
mkpath(DEST)
function features(K)
 v=sort(filter(>(0),K));L=length(v);M=sum(v)
 L==0 && return (links=0,support=0,mean_support=NaN,gini=NaN,singleton=NaN,top10=NaN)
 q=.9L;a=floor(Int,q);bottom=sum(@view v[1:a])+(a<L ? (q-a)*v[a+1] : 0.)
 (links=L,support=M,mean_support=M/L,gini=2dot(1:L,v)/(L*M)-(L+1)/L,singleton=mean(v.==1),top10=1-bottom/M)
end
function main()
 versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
 for r in eachrow(versions);@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256;end
 used=NamedTuple[]
 function load(path;kwargs...)
  push!(used,(path=relpath(path,ROOT),sha256=bytes2hex(sha256(read(path)))));CSV.read(path,DataFrame;kwargs...)
 end
 qc=load(joinpath(SOURCE,"quality_control.csv"));fits=load(joinpath(SOURCE,"fitted_parameters.csv"));old=load(joinpath(SOURCE,"conditioned_M/support_frequencies.csv"))
 sites=load(joinpath(ROOT,"outputs/statistical description of local support/site_universe.csv");types=Dict(:site=>String))
 rows=NamedTuple[];exactrows=NamedTuple[];checks=NamedTuple[];inputs=NamedTuple[]
 for (di,ds) in enumerate(DATA)
  g=qc[qc.dataset.==ds,:];K=Int.(g.K);C=Int.(g.C);N=only(unique(g.N));L=length(K);M=sum(K)
  ss=String.(sites.site[sites.dataset.==ds]);@assert length(ss)==N;si=Dict(s=>j for (j,s) in enumerate(ss));li=Dict((r.consumer,r.resource)=>i for (i,r) in enumerate(eachrow(g)))
  occ=load(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv");types=String);ints=load(joinpath(ROOT,"data/processed",ds*"_interactions.csv");types=String)
  pres=Dict((v.trophic_level[1],v.species[1])=>Set(v.site) for v in groupby(occ,[:trophic_level,:species]))
  O=[sort([si[s] for s in intersect(pres[("consumer",r.consumer)],pres[("resource",r.resource)])]) for r in eachrow(g)]
  empirical=[Int[] for _ in 1:N];seen=Set{Tuple{Int,Int}}()
  for r in eachrow(ints);i=li[(r.consumer,r.resource)];j=si[r.site];if !((i,j) in seen);push!(seen,(i,j));push!(empirical[j],i);end;end
  @assert length.(O)==C && [count(t->first(t)==i,seen) for i in 1:L]==K && all(j in O[i] for (i,j) in seen)
  push!(inputs,(dataset=ds,N=N,L=L,M=M,total_known_opportunities=sum(C)))
  cache=Dict{String,Any}();H=gapmatrix(N)
  for model in MODELS
   if model=="Empirical";f=[mean(K.==k) for k in 1:N]
   else
    p=only(eachrow(fits[(fits.dataset.==ds).&(fits.scope.=="LODO").&(fits.model.==model),:]));@assert Set(split(p.training_datasets,";"))==Set(filter(!=(ds),DATA))
    D,P=dp(pairweights(C,p.a,p.h),M);cache[model]=(D,P);f=exactfreq(M,D,P,N)
    @assert maximum(abs.(f.-old.predicted_frequency[(old.dataset.==ds).&(old.model.==model)]))<1e-8
   end
   for n in 0:N
    expectedL=L*(dot(H[N-n+1,:],f)+n/N)
    push!(exactrows,(dataset=ds,model=model,n=n,N=N,site_fraction=n/N,expected_links=expectedL,expected_link_fraction=expectedL/L,expected_total_support=n*M/N))
   end
  end
  for rep in 1:B
   order=randperm(MersenneTwister(991000+di*1000+rep),N)
   for (mi,model) in enumerate(MODELS)
    if model=="Empirical";cells=empirical
    else
     rng=MersenneTwister(1991000+di*10000+mi*1000+rep);D,P=cache[model];kk=drawK(rng,M,D,P);@assert sum(kk)==M && all(1 .<=kk .<=C)
     cells=[Int[] for _ in 1:N]
     for i in 1:L,j in O[i][randperm(rng,C[i])[1:kk[i]]];push!(cells[j],i);end
    end
    k=zeros(Int,L)
    for (n,j) in enumerate(order)
     for i in cells[j];k[i]+=1;end
     z=features(k);z.links>0 && @assert isapprox(z.mean_support*z.links,z.support;atol=1e-8)
     push!(rows,merge((dataset=ds,model=model,replicate=rep,n=n,N=N,site_fraction=n/N,regional_links=L,total_support=M),z))
    end
    @assert sum(k)==M && all(k.>0) && maximum(k)<=N
   end
  end
  println(ds," accumulation complete");flush(stdout)
  CSV.write(joinpath(DEST,"accumulation_draws.csv"),DataFrame(rows));CSV.write(joinpath(DEST,"exact_expectations.csv"),DataFrame(exactrows))
 end
 tab=DataFrame(rows);ex=DataFrame(exactrows);summ=NamedTuple[];diff=NamedTuple[]
 grouped=Dict((v.dataset[1],v.model[1],v.n[1])=>v for v in groupby(tab,[:dataset,:model,:n]))
 for v in groupby(tab,[:dataset,:model,:n])
  ds=v.dataset[1];model=v.model[1];n=v.n[1];N=v.N[1]
  for field in [:links,:support,:mean_support,:gini,:singleton,:top10]
   val=Float64.(v[!,field]);val=filter(isfinite,val);lo,hi=isempty(val) ? (NaN,NaN) : quantile(val,[.025,.975]);mu=isempty(val) ? NaN : mean(val);se=isempty(val) ? NaN : std(val)/sqrt(length(val))
   push!(summ,(dataset=ds,model=model,n=n,N=N,site_fraction=n/N,metric=String(field),mean=mu,lower=lo,upper=hi,mcse=se,valid_replicates=length(val),empty_replicates=count(v.links.==0)))
  end
 end
 st=DataFrame(summ)
 for ds in DATA,model in MODELS
  for field in ["links","support"]
   v=st[(st.dataset.==ds).&(st.model.==model).&(st.metric.==field),:];truth=ex[(ex.dataset.==ds).&(ex.model.==model).&(ex.n.>0),field=="links" ? :expected_links : :expected_total_support]
   @assert abs(v.mean[end]-truth[end])<1e-8
   push!(checks,(dataset=ds,model=model,quantity=field,maximum_mean_vs_exact_error=maximum(abs.(v.mean.-truth)),full_endpoint_verified=true))
  end
  model=="Empirical" && continue
  for n in unique(tab.n[tab.dataset.==ds]),field in [:links,:support,:mean_support,:gini,:singleton,:top10]
   ob=grouped[(ds,"Empirical",n)];mo=grouped[(ds,model,n)];@assert ob.replicate==mo.replicate
   values=filter(isfinite,Float64.(ob[!,field].-mo[!,field]));lo,hi=isempty(values) ? (NaN,NaN) : quantile(values,[.025,.975])
   push!(diff,(dataset=ds,model=model,n=n,site_fraction=n/ob.N[1],metric=String(field),mean_empirical_minus_model=isempty(values) ? NaN : mean(values),lower=lo,upper=hi,valid_replicates=length(values)))
  end
 end
 CSV.write(joinpath(DEST,"accumulation_summary.csv"),st);CSV.write(joinpath(DEST,"paired_differences.csv"),DataFrame(diff));CSV.write(joinpath(DEST,"checks.csv"),DataFrame(checks));CSV.write(joinpath(DEST,"dataset_inputs.csv"),DataFrame(inputs))
 CSV.write(joinpath(DEST,"input_versions.csv"),versions)
 helper=joinpath(ROOT,"scripts/beta_pair_conditioned_M/conditional_support.jl");push!(used,(path=relpath(helper,ROOT),sha256=bytes2hex(sha256(read(helper)))))
 CSV.write(joinpath(DEST,"source_versions.csv"),DataFrame(used));CSV.write(joinpath(DEST,"parameters_used.csv"),fits[(fits.scope.=="LODO").&in.(fits.model,Ref(["Common","Pair"])),:])
 for file in ["metadata_relationships.csv","metadata_source_versions.csv"];cp(joinpath(SOURCE,file),joinpath(DEST,file);force=true);end
end
include("figures.jl")
if !(length(ARGS)>=3 && ARGS[3]=="plot");main();end
plots();report()
end

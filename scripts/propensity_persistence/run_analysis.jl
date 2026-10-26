module PropensityPersistence
using CSV,DataFrames,Statistics,Random,LinearAlgebra,SHA,CairoMakie,Printf
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"../..")) : abspath(ARGS[1])
const SOURCE=joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory")
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(SOURCE,"propensity_persistence")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const GRID=[0.,.1,.3,.6,.9]
const B=300
logistic(x)=1/(1+exp(-x));logz(x)=log(-expm1(x));lse(x)=maximum(x)+log(sum(exp.(x.-maximum(x))))
include(joinpath(ROOT,"scripts/beta_pair_conditioned_M/conditional_support.jl"))
mkpath(DEST)
function weights(C,family,rho)
 [begin
  lc=[sum(log((c-j+1)/j) for j in 1:k) for k in 1:c]
  if rho==0
   lc.-c*log(2)
  elseif family=="Two-state"
   p=(1-sqrt(rho))/2
   [lc[k]+lse([k*log(p)+(c-k)*log1p(-p),k*log1p(-p)+(c-k)*log(p)])-log(2) for k in 1:c]
  else
   a=(1/rho-1)/2
   [lc[k]+sum(log(a+j) for j in 0:k-1)+sum(log(a+j) for j in 0:c-k-1;init=0.)-sum(log(2a+j) for j in 0:c-1) for k in 1:c]
  end
 end for c in C]
end
function features(K,C)
 v=sort(K);L=length(K);M=sum(K);q=.9L;j=floor(Int,q)
 groups=[findall(==(c),C) for c in unique(C) if count(==(c),C)>=3]
 within=sum(length(ix)*var(K[ix]./C[ix];corrected=false) for ix in groups)/sum(length,groups)
 (singleton=mean(K.==1),gini=2dot(1:L,v)/(L*M)-(L+1)/L,top10=1-(sum(v[1:j])+(q-j)*v[j+1])/M,within_C_variance=within)
end
function tests()
 C=[2,3,4];M=6
 for fam in ["Two-state","Beta"],rho in GRID
  w=weights(C,fam,rho);D,P=dp(w,M);f=exactfreq(M,D,P,4)
  brute=zeros(4);z=0.
  for a in 1:2,b in 1:3,c in 1:4
   a+b+c==M || continue
   mass=exp(w[1][a]+w[2][b]+w[3][c]);z+=mass
   for k in (a,b,c);brute[k]+=mass/3;end
  end
  @assert maximum(abs.(f.-brute/z))<1e-10
  c=7;positive=exp.(weights([c],fam,rho)[1])
  mu=dot(1:c,positive);variance=dot((1:c).^2,positive)-mu^2
  @assert isapprox(mu,c/2;atol=1e-10)
  @assert isapprox(variance,c/4+c*(c-1)*rho/4;atol=1e-9)
 end
 @assert weights(C,"Beta",0.)==weights(C,"Two-state",0.)
end
function run()
 versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
 for r in eachrow(versions);@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256;end
 qc=CSV.read(joinpath(SOURCE,"quality_control.csv"),DataFrame)
 prior=CSV.read(joinpath(SOURCE,"conditioned_M/support_frequencies.csv"),DataFrame)
 summaries=NamedTuple[];freqs=NamedTuple[];draws=NamedTuple[];checks=NamedTuple[]
 for (di,ds) in enumerate(DATA)
  g=qc[qc.dataset.==ds,:];K=Int.(g.K);C=Int.(g.C);N=only(unique(g.N));M=sum(K);L=length(K)
  @assert all(1 .<= K .<= C .<= N)
  emp=features(K,C)
  for (fi,fam) in enumerate(["Two-state","Beta"]),rho in GRID
   D,P=dp(weights(C,fam,rho),M);f=exactfreq(M,D,P,N)
   if rho==0
    old=prior[(prior.dataset.==ds).&(prior.model.=="Common"),:];sort!(old,:K)
    @assert maximum(abs.(f.-old.predicted_frequency))<1e-8
   end
   for k in 1:N;push!(freqs,(dataset=ds,family=fam,rho=rho,K=k,empirical=mean(K.==k),expected=f[k]));end
   vals=NamedTuple[];rng=MersenneTwister(721000+1000di+100fi+round(Int,100rho))
   for b in 1:B
    kk=drawK(rng,M,D,P);@assert sum(kk)==M && all(1 .<=kk .<=C)
    v=features(kk,C);push!(vals,v);push!(draws,merge((dataset=ds,family=fam,rho=rho,replicate=b),v))
   end
   for metric in keys(emp)
    v=[x[metric] for x in vals];lo,hi=quantile(v,[.025,.975])
    push!(summaries,(dataset=ds,N=N,L=L,M=M,family=fam,rho=rho,metric=String(metric),empirical=emp[metric],mean=mean(v),lower=lo,upper=hi,mcse=std(v)/sqrt(B)))
   end
   push!(checks,(dataset=ds,family=fam,rho=rho,probability_sum=sum(f),expected_M=L*dot(1:N,f),constraints_pass=true))
   println(ds," ",fam," rho=",rho);flush(stdout)
  end
  CSV.write(joinpath(DEST,"summary.csv"),DataFrame(summaries));CSV.write(joinpath(DEST,"frequencies.csv"),DataFrame(freqs));CSV.write(joinpath(DEST,"draws.csv"),DataFrame(draws));CSV.write(joinpath(DEST,"checks.csv"),DataFrame(checks))
 end
 CSV.write(joinpath(DEST,"input_versions.csv"),versions)
 paths=[joinpath(SOURCE,"quality_control.csv"),joinpath(SOURCE,"conditioned_M/support_frequencies.csv"),joinpath(ROOT,"scripts/beta_pair_conditioned_M/conditional_support.jl")]
 CSV.write(joinpath(DEST,"source_versions.csv"),DataFrame(path=relpath.(paths,ROOT),sha256=[bytes2hex(sha256(read(p))) for p in paths]))
 for file in ["metadata_relationships.csv","metadata_source_versions.csv"];cp(joinpath(SOURCE,file),joinpath(DEST,file);force=true);end
end
include("figures.jl")
tests()
if !(length(ARGS)>2 && ARGS[3]=="plot");run();end
plots()
end

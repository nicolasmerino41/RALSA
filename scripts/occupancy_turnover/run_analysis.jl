module OccupancyTurnover
using CSV,DataFrames,Statistics,SHA,CairoMakie,Printf
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"../..")) : abspath(ARGS[1])
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory/occupancy_turnover")
const SOURCE=joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
mkpath(DEST)
function components(a,r,t)
 den=2a+r+t;den==0 && return (NaN,NaN,NaN,false)
 wn=(r+t)/den;defined=2a+r>0;os=defined ? r/(2a+r) : 0.
 (wn,os,wn-os,defined)
end
function main()
 @assert components(2,0,1)[1:3]==(.2,0.,.2)
 @assert components(2,1,0)[1:3]==(.2,.2,0.)
 @assert components(0,1,4)[1:3]==(1.,1.,0.)
 @assert components(0,0,4)[1:3]==(1.,0.,1.)
 versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
 for r in eachrow(versions);@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256;end
 used=NamedTuple[]
 function load(path;kwargs...)
  push!(used,(path=relpath(path,ROOT),sha256=bytes2hex(sha256(read(path)))));CSV.read(path,DataFrame;kwargs...)
 end
 qc=load(joinpath(SOURCE,"quality_control.csv"));sites=load(joinpath(ROOT,"outputs/statistical description of local support/site_universe.csv");types=Dict(:site=>String))
 stages=NamedTuple[];species=NamedTuple[];links=NamedTuple[];pairs=NamedTuple[];sums=NamedTuple[];distributions=NamedTuple[]
 for ds in DATA
  g=qc[qc.dataset.==ds,:];N=only(unique(g.N));L=nrow(g);ss=sites.site[sites.dataset.==ds];@assert length(ss)==N
  si=Dict(s=>j for (j,s) in enumerate(ss));li=Dict((r.consumer,r.resource)=>i for (i,r) in enumerate(eachrow(g)))
  occ=load(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv");types=String)
  ints=load(joinpath(ROOT,"data/processed",ds*"_interactions.csv");types=String)
  pres=Dict((v.trophic_level[1],v.species[1])=>Set(si[s] for s in v.site) for v in groupby(occ,[:trophic_level,:species]))
  degree=Dict(key=>0 for key in keys(pres));O=BitSet[];I=[BitSet() for _ in 1:N];S=[Set{Tuple{String,String}}() for _ in 1:N]
  for (key,locs) in pres,j in locs;push!(S[j],key);end
  endpoints=Tuple{Tuple{String,String},Tuple{String,String}}[]
  for r in eachrow(g)
   a=("consumer",r.consumer);b=("resource",r.resource);degree[a]+=1;degree[b]+=1;push!(endpoints,(a,b));push!(O,BitSet(intersect(pres[a],pres[b])))
  end
  for r in eachrow(ints);push!(I[si[r.site]],li[(r.consumer,r.resource)]);end
  @assert [length(o) for o in O]==g.C
  @assert [count(i in ix for ix in I) for i in 1:L]==g.K
  @assert all(j in O[i] for j in 1:N for i in I[j])
  for (key,locs) in pres;push!(species,(dataset=ds,role=key[1],species=key[2],occupancy=length(locs),degree=degree[key]));end
  s=DataFrame(species);s=s[s.dataset.==ds,:]
  for role in ["consumer","resource","combined"]
   v=role=="combined" ? s : s[s.role.==role,:];linked=v[v.degree.>0,:]
   a=mean(v.occupancy.==1);b=mean(linked.occupancy.==1);w=sum(v.degree[v.occupancy.==1])/sum(v.degree)
   fc=mean(g.C.==1);fk=mean(g.K.==1)
   push!(stages,(dataset=ds,role=role,N=N,L=L,species=nrow(v),species_with_links=nrow(linked),singleton_species=a,singleton_linked_species=b,singleton_degree_weighted=w,singleton_C=fc,singleton_K=fk,high_low_high=a>fc && fk>fc,high_low_high_linked=b>fc && fk>fc))
   for n in 1:N
    push!(distributions,(dataset=ds,role=role,n=n,unweighted=mean(v.occupancy.==n),linked_only=mean(linked.occupancy.==n),degree_weighted=sum(v.degree[v.occupancy.==n])/sum(v.degree)))
   end
  end
  category=zeros(Int,3)
  for (i,r) in enumerate(eachrow(g))
   a,b=endpoints[i];oa=length(pres[a]);ob=length(pres[b])
   cat=r.K!=1 ? "Repeated" : min(oa,ob)==1 ? "Single-site species" : r.C==1 ? "Single-site overlap" : "Repeated opportunities"
   cat!="Repeated" && (category[cat=="Single-site species" ? 1 : cat=="Single-site overlap" ? 2 : 3]+=1)
   push!(links,(dataset=ds,consumer=r.consumer,resource=r.resource,consumer_occupancy=oa,resource_occupancy=ob,C=r.C,K=r.K,category=cat))
  end
  @assert sum(category)==count(g.K.==1)
  @assert category[1]+category[2]==count(g.C.==1)
  localpairs=NamedTuple[]
  for j in 1:N-1,k in j+1:N
   shared=intersect(I[j],I[k]);unique_links=symdiff(I[j],I[k]);a=length(shared);r=0
   for idx in unique_links
    u,v=endpoints[idx]
    r+=(u in S[j] && u in S[k] && v in S[j] && v in S[k])
   end
   t=length(unique_links)-r;wn,os,st,defined=components(a,r,t)
   ns=length(S[j])+length(S[k]);bs=ns>0 ? 1-2length(intersect(S[j],S[k]))/ns : NaN
   @assert !isfinite(wn) || (abs(wn-os-st)<1e-10 && st>=-1e-10)
   push!(localpairs,(dataset=ds,site_a=ss[j],site_b=ss[k],shared_links=a,shared_species_link_differences=r,species_turnover_link_differences=t,beta_species=bs,beta_WN=wn,beta_OS=os,beta_ST=st,OS_defined=defined))
  end
  append!(pairs,localpairs);p=DataFrame(localpairs);valid=isfinite.(p.beta_WN);q=p[valid,:];def=q[q.OS_defined,:]
  totaldiff=sum(q.shared_species_link_differences)+sum(q.species_turnover_link_differences)
  push!(sums,(dataset=ds,N=N,L=L,singletons=sum(category),single_site_species=category[1],single_site_overlap=category[2],repeated_opportunities=category[3],site_pairs=nrow(p),nonempty_pairs=nrow(q),OS_undefined_pairs=count(.!q.OS_defined),mean_beta_species=mean(filter(isfinite,p.beta_species)),mean_beta_WN=mean(q.beta_WN),mean_beta_OS=mean(q.beta_OS),mean_beta_ST=mean(q.beta_ST),mean_beta_WN_OS_defined=isempty(def) ? NaN : mean(def.beta_WN),mean_beta_OS_OS_defined=isempty(def) ? NaN : mean(def.beta_OS),mean_beta_ST_OS_defined=isempty(def) ? NaN : mean(def.beta_ST),turnover_fraction_of_differing_links=totaldiff>0 ? sum(q.species_turnover_link_differences)/totaldiff : NaN))
  println(ds," complete");flush(stdout)
 end
 for (name,data) in [("singleton_stages",stages),("species",species),("link_categories",links),("site_pair_dissimilarity",pairs),("dataset_summary",sums),("occupancy_distributions",distributions)];CSV.write(joinpath(DEST,name*".csv"),DataFrame(data));end
 CSV.write(joinpath(DEST,"input_versions.csv"),versions);CSV.write(joinpath(DEST,"source_versions.csv"),DataFrame(used))
 for file in ["metadata_relationships.csv","metadata_source_versions.csv"];cp(joinpath(SOURCE,file),joinpath(DEST,file);force=true);end
end
include("figures.jl")
if !(length(ARGS)>2 && ARGS[3]=="plot");main();end
plots()
end

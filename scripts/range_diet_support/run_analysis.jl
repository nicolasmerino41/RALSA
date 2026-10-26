using CSV, DataFrames, Statistics, Random, LinearAlgebra, SHA, CairoMakie, Printf
root=length(ARGS)>0 ? abspath(ARGS[1]) : normpath(joinpath(@__DIR__,"../.."))
dest=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(root,"outputs/statistical description of local support/range_diet_support")
mkpath(dest)
base=joinpath(root,"outputs/statistical description of local support/opportunity_theory")
datasets=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
sources=NamedTuple[]
function load(path; kwargs...)
 push!(sources,(path=relpath(path,root),sha256=bytes2hex(sha256(read(path)))))
 CSV.read(path,DataFrame;kwargs...)
end
function ranks(x)
 p=sortperm(x);r=zeros(length(x));i=1
 while i<=length(x)
  j=i;while j<length(x) && x[p[j+1]]==x[p[i]];j+=1;end
  r[p[i:j]].=(i+j)/2;i=j+1
 end
 r
end
rho(x,y)=length(x)>2 && std(x)>0 && std(y)>0 ? cor(ranks(x),ranks(y)) : NaN
function partial(x,y,z)
 X=hcat(ones(length(x)),ranks(z));a=ranks(x);b=ranks(y)
 ra=a-X*(X\a);rb=b-X*(X\b)
 norm(ra)>1e-8 && norm(rb)>1e-8 ? cor(ra,rb) : NaN
end
qc=load(joinpath(base,"quality_control.csv"))
universe=load(joinpath(root,"outputs/statistical description of local support/site_universe.csv");types=Dict(:site=>String))
versions=load(joinpath(base,"input_versions.csv"))
for r in eachrow(versions)
 @assert bytes2hex(sha256(read(joinpath(root,r.path))))==r.sha256
end
species=NamedTuple[];links=NamedTuple[];stats=NamedTuple[];groups=NamedTuple[];origins=NamedTuple[];tests=NamedTuple[];sensitivity=NamedTuple[]
rng=MersenneTwister(20261019);B=1999
for ds in datasets
 g=qc[qc.dataset.==ds,:];N=only(unique(g.N));L=nrow(g)
 occ=load(joinpath(root,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv");types=String)
 ints=load(joinpath(root,"data/processed",ds*"_interactions.csv");types=String)
 ss=Set(universe.site[universe.dataset.==ds]);@assert length(ss)==N
 @assert all(in(ss),occ.site) && all(in(ss),ints.site)
 pres=Dict((v.trophic_level[1],v.species[1])=>Set(v.site) for v in groupby(occ,[:trophic_level,:species]))
 observed=Dict((v.consumer[1],v.resource[1])=>Set(v.site) for v in groupby(ints,[:consumer,:resource]))
 @assert length(observed)==L
 deg=Dict(k=>0 for k in keys(pres))
 for r in eachrow(g)
  a=("consumer",r.consumer);b=("resource",r.resource)
  @assert length(intersect(pres[a],pres[b]))==r.C
  @assert length(observed[(r.consumer,r.resource)])==r.K
  @assert issubset(observed[(r.consumer,r.resource)],intersect(pres[a],pres[b]))
  @assert 1<=r.K<=r.C<=N
  deg[a]+=1;deg[b]+=1
 end
 local_species=NamedTuple[]
 for ((role,sp),sites) in pres
  ix=role=="consumer" ? g.consumer.==sp : g.resource.==sp
  d=deg[(role,sp)];linked=d>0
  push!(local_species,(dataset=ds,role=role,species=sp,N=N,occupancy=length(sites),occupancy_fraction=length(sites)/N,degree=d,mean_K=linked ? mean(g.K[ix]) : NaN,mean_C=linked ? mean(g.C[ix]) : NaN,mean_q=linked ? mean(g.K[ix]./g.C[ix]) : NaN,singleton_fraction=linked ? mean(g.K[ix].==1) : NaN))
 end
 append!(species,local_species);s=DataFrame(local_species)
 oc=[length(pres[("consumer",r.consumer)]) for r in eachrow(g)]
 ore=[length(pres[("resource",r.resource)]) for r in eachrow(g)]
 dc=[deg[("consumer",r.consumer)] for r in eachrow(g)]
 dr=[deg[("resource",r.resource)] for r in eachrow(g)]
 for (i,r) in enumerate(eachrow(g))
  category=min(oc[i],ore[i])==1 ? "One-site endpoint" : r.C==1 ? "One-site overlap" : "Repeated opportunities"
  push!(links,(dataset=ds,consumer=r.consumer,resource=r.resource,N=N,K=r.K,C=r.C,q=r.K/r.C,consumer_occupancy=oc[i],resource_occupancy=ore[i],consumer_degree=dc[i],resource_degree=dr[i],minimum_occupancy=min(oc[i],ore[i]),singleton_origin=r.K==1 ? category : "Repeated link"))
 end
 for role in ["consumer","resource"]
  t=s[(s.role.==role).&(s.degree.>0),:];o=role=="consumer" ? oc : ore;d=role=="consumer" ? dc : dr
  push!(stats,(dataset=ds,role=role,N=N,recorded_species=count(s.role.==role),linked_species=nrow(t),rho_degree_occupancy=rho(t.degree,t.occupancy),rho_degree_mean_K=rho(t.degree,t.mean_K),rho_occupancy_mean_K=rho(t.occupancy,t.mean_K),rho_degree_mean_K_given_occupancy=partial(t.degree,t.mean_K,t.occupancy),rho_degree_mean_K_given_mean_C=partial(t.degree,t.mean_K,t.mean_C),rho_degree_mean_q=rho(t.degree,t.mean_q)))
  for threshold in [1,2,3]
   low=t.degree.<=threshold;edge=d.<=threshold;single=g.K.==1
   push!(sensitivity,(dataset=ds,role=role,maximum_degree=threshold,species_fraction=mean(low),link_fraction=mean(edge),singleton_share=sum(edge.&single)/sum(single),singleton_rate=any(edge) ? mean(single[edge]) : NaN,other_singleton_rate=any(.!edge) ? mean(single[.!edge]) : NaN,median_occupancy=any(low) ? median(t.occupancy[low]) : NaN,other_median_occupancy=any(.!low) ? median(t.occupancy[.!low]) : NaN))
  end
  labels=["1","2","3–4","5–8","9+"]
  bin(x)=x==1 ? 1 : x==2 ? 2 : x<=4 ? 3 : x<=8 ? 4 : 5
  for j in 1:5
   a=bin.(t.degree).==j;b=bin.(d).==j
   push!(groups,(dataset=ds,role=role,degree_class=labels[j],class_index=j,species_count=sum(a),species_fraction=mean(a),link_count=sum(b),link_fraction=mean(b),median_occupancy_percent=any(a) ? 100median(t.occupancy[a])/N : NaN,mean_K=any(b) ? mean(g.K[b]) : NaN,singleton_rate=any(b) ? mean(g.K[b].==1) : NaN,singleton_share=sum(b.&(g.K.==1))/sum(g.K.==1)))
  end
  strata=[findall(g.C.==c) for c in unique(g.C)]
  eligible=reduce(vcat,[ix for ix in strata if length(unique(g.K[ix]))>1];init=Int[])
  # Conditional shuffle preserves each C class and the full support distribution.
  low=d.<=2
  contrast(k)=any(low) && any(.!low) ? mean(k[low].==1)-mean(k[.!low].==1) : NaN
  obs=[rho(d,g.K),contrast(g.K)];null=zeros(B,2);kp=copy(g.K)
  for b in 1:B
   for ix in strata;kp[ix]=shuffle(rng,g.K[ix]);end
   null[b,:]=[rho(d,kp),contrast(kp)]
  end
  for j in 1:2
   z=null[:,j];mu=mean(z);p=isfinite(obs[j]) ? (1+count(abs.(z.-mu).>=abs(obs[j]-mu)-1e-12))/(B+1) : NaN
   push!(tests,(dataset=ds,role=role,metric=j==1 ? "degree_K_rho" : "low_degree_singleton_excess",observed=obs[j],null_mean=mu,null_lower=quantile(z,.025),null_upper=quantile(z,.975),effect=obs[j]-mu,p=p,movable_link_fraction=length(eligible)/L,permutations=B))
  end
 end
 for threshold in [1,2,3]
  rare=g.K.<=threshold;a=min.(oc,ore).<=threshold;b=(g.C.<=threshold).&(.!a);c=g.C.>threshold
  @assert all((a.|b.|c)) && !any(a.&b) && all(g.K[a.|b].<=threshold)
  for (label,ix) in zip(["Restricted endpoint","Restricted overlap","Beyond opportunity restriction"],[a,b,c])
   push!(origins,(dataset=ds,maximum_K=threshold,category=label,rare_links=count(rare),count=sum(rare.&ix),share_of_rare=sum(rare.&ix)/sum(rare),share_of_all=sum(rare.&ix)/L))
  end
 end
 println(ds," verified and analysed");flush(stdout)
end
td=DataFrame(tests);td.adjusted_p=fill(NaN,nrow(td));valid=findall(isfinite,td.p);order=valid[sortperm(td.p[valid])]
let q=1.
 for j in length(order):-1:1;q=min(q,td.p[order[j]]*length(order)/j);td.adjusted_p[order[j]]=q;end
end
for (name,data) in [("species",species),("links",links),("species_associations",stats),("degree_groups",groups),("rare_link_origins",origins),("specialist_thresholds",sensitivity)]
 CSV.write(joinpath(dest,name*".csv"),DataFrame(data))
end
CSV.write(joinpath(dest,"conditional_tests.csv"),td);CSV.write(joinpath(dest,"input_versions.csv"),versions)
CSV.write(joinpath(dest,"source_versions.csv"),DataFrame(sources))
for file in ["metadata_relationships.csv","metadata_source_versions.csv"]
 cp(joinpath(base,file),joinpath(dest,file);force=true)
end
include("plots.jl")
include("report.jl")

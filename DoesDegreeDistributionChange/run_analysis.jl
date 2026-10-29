module DegreeChange
using CSV,DataFrames,Statistics,Random,SHA,CairoMakie,Printf
const ROOT=length(ARGS)>0 ? abspath(ARGS[1]) : dirname(@__DIR__)
const OUT=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(@__DIR__,"outputs")
const B=length(ARGS)>2 ? parse(Int,ARGS[3]) : 2000
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const BASE=joinpath(ROOT,"outputs/statistical description of local support/K_CO_SP_distributions")
const GRID=collect(0.:.05:5.)
function gini(x)
 v=sort(x[x.>0]);n=length(v)
 n==0 ? NaN : 2sum((1:n).*v)/(n*sum(v))-(n+1)/n
end
function ranks(x)
 p=sortperm(x);r=zeros(length(x));i=1
 while i<=length(x)
  j=i;while j<length(x)&&x[p[j+1]]==x[p[i]];j+=1;end
  r[p[i:j]].=(i+j)/2;i=j+1
 end;r
end
rho(x,y)=std(x)>0&&std(y)>0 ? cor(ranks(x),ranks(y)) : NaN
function survival(N,m,K)
 e=ones(N+1)
 for k in 1:N;e[k+1]=k>m ? 0. : e[k]*(m-k+1)/(N-k+1);end
 1 .- e[K.+1]
end
function degrees(ids,keep,S)
 d=zeros(Int,S)
 for j in eachindex(ids);d[ids[j]]+=keep[j];end
 d
end
function main()
 mkpath(OUT);rng=MersenneTwister(20261028);sources=NamedTuple[]
 function load(path;kwargs...)
  push!(sources,(path=relpath(path,ROOT),sha256=bytes2hex(sha256(read(path)))));CSV.read(path,DataFrame;kwargs...)
 end
 qc=load(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/pairs.csv"));qc=qc[qc.K.>0,:]
 versions=load(joinpath(BASE,"input_versions.csv"))
 hashchecks=NamedTuple[]
 for r in eachrow(versions)
  path=joinpath(ROOT,r.path);present=isfile(path)
  required=endswith(r.path,"_interactions.csv")||endswith(r.path,"_occupancy.csv")
  @assert !required||present
  matches=present ? bytes2hex(sha256(read(path)))==r.sha256 : false
  @assert !present||matches
  push!(hashchecks,(path=r.path,required_input=required,present=present,hash_matches=matches))
 end
 CSV.write(joinpath(OUT,"historical_manifest_checks.csv"),DataFrame(hashchecks))
 @assert isapprox(gini([1,1]),0) && isapprox(gini([1,3]),.25)
 @assert survival(10,1,[1,5])≈[.9,1.]
 sp=NamedTuple[];sim=NamedTuple[];cc=NamedTuple[];full=NamedTuple[];checks=NamedTuple[]
 for ds in DATA
  v=qc[qc.dataset.==ds,:];N=only(unique(v.N));L=nrow(v)
  occupancy=load(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv");types=String)
  ss=sort(unique(occupancy.site));@assert length(ss)==N
  si=Dict(s=>i for (i,s) in enumerate(ss));li=Dict((r.consumer,r.resource)=>i for (i,r) in enumerate(eachrow(v)))
  ints=load(joinpath(ROOT,"data/processed",ds*"_interactions.csv");types=String)
  mat=falses(N,L)
  for r in eachrow(ints);mat[si[r.site],li[(r.consumer,r.resource)]]=true;end
  @assert vec(sum(mat;dims=1))==v.K
  names=[sort(unique(v.consumer)),sort(unique(v.resource))]
  ids=[[findfirst(==(x),names[1]) for x in v.consumer],[findfirst(==(x),names[2]) for x in v.resource]]
  initial=[degrees(ids[h],trues(L),length(names[h])) for h in 1:2]
  masks=[initial[h].<=median(initial[h]) for h in 1:2]
  fractions=[.4,.8];ms=[round(Int,f*N) for f in fractions]
  probs=[survival(N,m,v.K) for m in ms]
  # Every row is a site; shuffling columns preserves all site totals and K values.
  for (j,m) in enumerate(ms),h in 1:2
   d=initial[h];ex=zeros(length(d))
   for l in 1:L;ex[ids[h][l]]+=probs[j][l];end
   for i in eachindex(d)
    push!(sp,(dataset=ds,role=h==1 ? "consumer" : "resource",species=names[h][i],N=N,target_removal=fractions[j],m=m,actual_removal=m/N,initial_degree=d[i],expected_degree=ex[i],expected_fraction=ex[i]/d[i],shuffle_expected_fraction=mean(probs[j]),initial_degree_group=masks[h][i] ? "lower" : "higher"))
   end
  end
  # Streaming moments for scaled CCDF; complete integer-degree histograms are saved separately.
  acc=Dict((j,h,s)=>(zeros(length(GRID)),zeros(length(GRID)),Ref(0)) for j in 1:2 for h in 1:2 for s in ["Observed","Shuffled"])
  hist=Dict((j,h,s)=>zeros(Int,maximum(initial[h])+1) for j in 1:2 for h in 1:2 for s in ["Observed","Shuffled"])
  degsum=Dict((j,h)=>zeros(length(initial[h])) for j in 1:2 for h in 1:2)
  for b in 1:B
   order=randperm(rng,N);perm=randperm(rng,L)
   if b==1
    @assert vec(sum(mat[:,perm];dims=2))==vec(sum(mat;dims=2))
    @assert sort(vec(sum(mat[:,perm];dims=1)))==sort(v.K)
   end
   for (j,m) in enumerate(ms)
    kept=vec(any(mat[order[m+1:end],:];dims=1));shuffled=kept[perm]
    @assert sum(kept)==sum(shuffled)
    for h in 1:2,(scenario,which) in [("Observed",kept),("Shuffled",shuffled)]
     d=degrees(ids[h],which,length(initial[h]));active=d[d.>0];lo=masks[h];hi=.!lo
     scenario=="Observed" && (degsum[(j,h)].+=d)
     push!(sim,(dataset=ds,role=h==1 ? "consumer" : "resource",target_removal=fractions[j],m=m,actual_removal=m/N,replicate=b,scenario=scenario,gini_active=gini(d),gini_initial=gini(initial[h]),disconnected_fraction=mean(d.==0),low_retained_fraction=mean(d[lo]./initial[h][lo]),high_retained_fraction=any(hi) ? mean(d[hi]./initial[h][hi]) : NaN,low_disconnected_fraction=mean(d[lo].==0),high_disconnected_fraction=any(hi) ? mean(d[hi].==0) : NaN,active_species=length(active),regional_links=sum(kept)))
     for k in d;hist[(j,h,scenario)][k+1]+=1;end
     if !isempty(active)
      y=[mean(active./mean(active).>=x) for x in GRID];a,a2,n=acc[(j,h,scenario)];a.+=y;a2.+=y.^2;n[]+=1
     end
    end
   end
  end
  for j in 1:2,h in 1:2
   d=initial[h];expected=zeros(length(d))
   for l in 1:L;expected[ids[h][l]]+=probs[j][l];end
   push!(checks,(dataset=ds,role=h==1 ? "consumer" : "resource",target_removal=fractions[j],N=N,m=ms[j],actual_removal=ms[j]/N,linked_species=length(d),L=L,M=sum(v.K),replicates=B,max_simulated_expected_fraction_error=maximum(abs.(degsum[(j,h)]/B.-expected)./d),rho_initial_degree_expected_fraction=rho(d,expected./d),low_degree_cutoff=median(d)))
   for scenario in ["Observed","Shuffled"]
    a,a2,n=acc[(j,h,scenario)]
    for (ix,x) in enumerate(GRID)
     avg=n[]>0 ? a[ix]/n[] : NaN;se=n[]>1 ? sqrt(max(0.,(a2[ix]-a[ix]^2/n[])/(n[]-1))/n[]) : NaN
     push!(cc,(dataset=ds,role=h==1 ? "consumer" : "resource",target_removal=fractions[j],actual_removal=ms[j]/N,scenario=scenario,scaled_degree=x,ccdf=avg,monte_carlo_se=se,nonempty_replicates=n[],initial_ccdf=mean(d./mean(d).>=x)))
    end
    for k in 0:maximum(d);push!(full,(dataset=ds,role=h==1 ? "consumer" : "resource",target_removal=fractions[j],scenario=scenario,degree=k,mean_species_count=hist[(j,h,scenario)][k+1]/B));end
   end
  end
  println(ds," complete");flush(stdout)
 end
 sf=DataFrame(sim);summary=NamedTuple[];effects=NamedTuple[]
 for v in groupby(sf,[:dataset,:role,:target_removal,:scenario])
  z=v.gini_active[isfinite.(v.gini_active)];delta=z.-v.gini_initial[1]
  push!(summary,(dataset=v.dataset[1],role=v.role[1],target_removal=v.target_removal[1],actual_removal=v.actual_removal[1],scenario=v.scenario[1],initial_gini=v.gini_initial[1],mean_gini=mean(z),mean_change=mean(delta),change_lower=quantile(delta,.025),change_upper=quantile(delta,.975),mean_disconnected=mean(v.disconnected_fraction),mean_low_retention=mean(v.low_retained_fraction),mean_high_retention=mean(v.high_retained_fraction),mean_low_disconnected=mean(v.low_disconnected_fraction),mean_high_disconnected=mean(v.high_disconnected_fraction),empty_replicates=count(v.active_species.==0),replicates=nrow(v)))
 end
 for v in groupby(sf,[:dataset,:role,:target_removal])
  a=sort(v[v.scenario.=="Observed",:],:replicate);b=sort(v[v.scenario.=="Shuffled",:],:replicate)
  d=a.gini_active.-b.gini_active;d=d[isfinite.(d)];se=std(d)/sqrt(length(d))
  push!(effects,(dataset=v.dataset[1],role=v.role[1],target_removal=v.target_removal[1],actual_removal=v.actual_removal[1],mean_gini_observed_minus_shuffled=mean(d),paired_sd=std(d),monte_carlo_se=se,mc_lower=mean(d)-1.96se,mc_upper=mean(d)+1.96se,paired_lower=quantile(d,.025),paired_upper=quantile(d,.975),valid_pairs=length(d)))
 end
 for (name,data) in [("species_expected_retention",sp),("removal_replicates",sim),("scaled_ccdf",cc),("raw_degree_histograms",full),("checks",checks),("summary",summary),("paired_effects",effects)]
  CSV.write(joinpath(OUT,name*".csv"),DataFrame(data))
 end
 CSV.write(joinpath(OUT,"source_versions.csv"),DataFrame(sources));CSV.write(joinpath(OUT,"input_versions.csv"),versions)
 for f in ["metadata_relationships.csv","metadata_source_versions.csv"];cp(joinpath(BASE,f),joinpath(OUT,f);force=true);end
 issues=joinpath(ROOT,"outputs/data_audit_20261026/issues_and_actions.csv");isfile(issues)&&cp(issues,joinpath(OUT,"outstanding_data_issues.csv");force=true)
 include("plots_and_report.jl")
end
main()
end

using CSV,DataFrames,Statistics,LinearAlgebra,Random,SHA,CairoMakie,Printf
const ROOT=normpath(joinpath(@__DIR__,"../.."))
const OUT=joinpath(@__DIR__,"outputs")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const MODELS=["Power_law","Cutoff_power_law","Discrete_lognormal","Exponential"]
const NAMES=["Power law","Truncated power law","Lognormal","Exponential"]
const COLORS=["#0072B2","#AA4499","#009E73","#D55E00"]
include(joinpath(ROOT,"scripts/K_CO_SP_distributions/distribution_helpers.jl"))
mkpath(OUT)
function main()
 manifest=CSV.read(joinpath(ROOT,"DoesDegreeDistributionChange/outputs/source_versions.csv"),DataFrame)
 for r in eachrow(manifest);@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256;end
 CSV.write(joinpath(OUT,"input_versions.csv"),manifest)
 write(joinpath(OUT,"helper_sha256.txt"),bytes2hex(sha256(read(joinpath(ROOT,"scripts/K_CO_SP_distributions/distribution_helpers.jl")))))
 pairs=CSV.read(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/pairs.csv"),DataFrame);pairs=pairs[pairs.K.>0,:]
 rng=MersenneTwister(20261028);details=NamedTuple[];states=NamedTuple[];orders=NamedTuple[]
 for ds in DATA
  p=pairs[pairs.dataset.==ds,:];N=only(unique(p.N));L=nrow(p)
  occ=CSV.read(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv"),DataFrame;types=String)
  sites=sort(unique(occ.site));@assert length(sites)==N
  si=Dict(s=>i for (i,s) in enumerate(sites));li=Dict((r.consumer,r.resource)=>i for (i,r) in enumerate(eachrow(p)))
  mat=falses(N,L)
  ints=CSV.read(joinpath(ROOT,"data/processed",ds*"_interactions.csv"),DataFrame;types=String)
  for r in eachrow(ints);mat[si[r.site],li[(r.consumer,r.resource)]]=true;end
  @assert vec(sum(mat;dims=1))==p.K
  consumers=sort(unique(p.consumer));ci=Dict(s=>i for (i,s) in enumerate(consumers));ids=[ci[s] for s in p.consumer];U=length(unique(p.resource))
  levels=sort(unique([round(Int,f*N) for f in 0:.05:.95]));filter!(m->m<N,levels)
  cache=Dict{Tuple,Any}()
  for b in 1:100
   order=randperm(rng,N)
   for (rank,j) in enumerate(order);push!(orders,(dataset=ds,replicate=b,rank=rank,site=sites[j]));end
   for m in levels
    keep=vec(any(mat[order[m+1:end],:];dims=1));degree=zeros(Int,length(consumers))
    for j in eachindex(ids);degree[ids[j]]+=keep[j];end
    degree=degree[degree.>0];n=length(degree);valid=n>=3 && length(unique(degree))>1
    winner=0;w=zeros(4);margin=NaN
    if valid
     counts=[count(==(k),degree) for k in 1:U]
     fs=get!(cache,Tuple(counts)) do;[fit(counts,model) for model in MODELS];end
     a=[z.aic for z in fs];@assert all(isfinite,a)
     w=exp.(-0.5 .*(a.-minimum(a)));w./=sum(w);winner=argmax(w);margin=sort(a)[2]-minimum(a)
     for (j,z) in enumerate(fs)
      push!(details,(dataset=ds,replicate=b,m=m,N=N,removal=m/N,n=n,model=NAMES[j],AIC=z.aic,weight=w[j],KS_D=z.D,alpha=z.alpha,lambda=z.lambda,mu=z.mu,sigma=z.sigma,boundary=z.boundary,kkt=z.kkt))
     end
    end
    push!(states,(dataset=ds,replicate=b,m=m,N=N,removal=m/N,connected_species=n,distinct_degrees=length(unique(degree)),valid=valid,winner=winner,winner_weight=maximum(w),delta_AIC_runner_up=margin))
   end
  end
  println(ds," complete");flush(stdout)
 end
 det=DataFrame(details);st=DataFrame(states)
 CSV.write(joinpath(OUT,"fits.csv"),det);CSV.write(joinpath(OUT,"states.csv"),st);CSV.write(joinpath(OUT,"removal_orders.csv"),DataFrame(orders))
 av=NamedTuple[]
 for g in groupby(st,[:dataset,:m])
  ds=g.dataset[1];m=g.m[1];fits=det[(det.dataset.==ds).&(det.m.==m),:]
  for name in NAMES
   sub=fits[fits.model.==name,:]
   push!(av,(dataset=ds,m=m,removal=g.removal[1],model=name,mean_weight=isempty(sub) ? NaN : mean(sub.weight),valid_sequences=count(g.valid),mean_species=mean(g.connected_species),low_species=quantile(g.connected_species,.025),high_species=quantile(g.connected_species,.975)))
  end
 end
 avg=DataFrame(av);CSV.write(joinpath(OUT,"mean_model_support.csv"),avg)
 transitions=NamedTuple[]
 for g in groupby(st,[:dataset,:replicate])
  sort!(g,:m);base=g.winner[1];idx=0
  for j in 2:nrow(g)-2
   if base>0 && g.winner[j]>0 && g.winner[j]!=base && all(g.winner[j:j+2].==g.winner[j]);idx=j;break;end
  end
  switches=sum(g.winner[j]>0 && g.winner[j-1]>0 && g.winner[j]!=g.winner[j-1] for j in 2:nrow(g))
  returnbase=idx>0 && any(g.winner[idx+3:end].==base)
  push!(transitions,(dataset=g.dataset[1],replicate=g.replicate[1],baseline=base,threshold=idx>0 ? g.removal[idx] : NaN,new_family=idx>0 ? g.winner[idx] : 0,confirmation=idx>0 ? g.removal[idx+2] : NaN,weight_at_threshold=idx>0 ? g.winner_weight[idx] : NaN,n_at_threshold=idx>0 ? g.connected_species[idx] : 0,switches=switches,returns_to_baseline=returnbase))
 end
 tr=DataFrame(transitions);CSV.write(joinpath(OUT,"sequence_transitions.csv"),tr)
 ts=DataFrame([(dataset=g.dataset[1],detected=count(isfinite,g.threshold),median_threshold=any(isfinite,g.threshold) ? median(filter(isfinite,g.threshold)) : NaN,low_threshold=any(isfinite,g.threshold) ? quantile(filter(isfinite,g.threshold),.25) : NaN,high_threshold=any(isfinite,g.threshold) ? quantile(filter(isfinite,g.threshold),.75) : NaN,median_switches=median(g.switches),return_count=count(g.returns_to_baseline)) for g in groupby(tr,:dataset)])
 CSV.write(joinpath(OUT,"transition_summary.csv"),ts)
 # Check original baseline, 40% and 80% likelihoods for identical sequences.
 old=CSV.read(joinpath(ROOT,"DoesDegreeDistributionChange/degree_model_selection/outputs/fit_diagnostics.csv"),DataFrame)
 checks=NamedTuple[]
 for r in eachrow(old)
  model=r.model=="Log-normal" ? "Lognormal" : r.model
  q=det[(det.dataset.==r.dataset).&(det.replicate.==r.replicate).&(abs.(det.removal.-r.actual_removal).<1e-10).&(det.model.==model),:]
  isempty(q) && continue
  push!(checks,(dataset=r.dataset,removal=r.removal,replicate=r.replicate,model=model,AIC_difference=only(q.AIC)-r.AIC))
 end
 ck=DataFrame(checks);CSV.write(joinpath(OUT,"previous_analysis_checks.csv"),ck);@assert maximum(abs.(ck.AIC_difference))<1e-6
 plots(st,avg,tr,ts)
 open(joinpath(OUT,"FINDINGS.md"),"w") do io
  println(io,"# Degree family transitions under removal\n\nConsumer degrees, same four bounded MLE models and fixed initial resource-count bound as the previous analysis. Same seed and 100 nested removal orders; saved original 0/40/80 AICs agree within 1e-6. Grid: nearest integer site counts to 0:5:95%, with duplicates and complete removal excluded. Actual fractions are used throughout.\n\nFits require at least three connected species and two distinct degrees; other states are grey and cannot establish a transition. No new KS bootstrap is run: this analysis concerns relative AIC support, not absolute goodness of fit. All fitted parameters, KS D, optimiser residuals and sample sizes are retained.\n\nA sustained departure is the first new family winning three successive evaluated levels. It can later reverse; returns and all switches are reported. Confirmation spans differ by site count. Threshold summaries condition on detected transitions; absence is not assigned 95%. Labels alone do not establish ecological regimes. Ties choose the first listed family; AIC weights expose uncertainty.\n")
  println(io,"| Dataset | Sustained departure /100 | Median onset % | IQR % | Median switches | Later returns /100 |\n|---|---:|---:|---|---:|---:|")
  for r in eachrow(ts);println(io,@sprintf("| %s | %d | %.1f | %.1f–%.1f | %.1f | %d |",r.dataset,r.detected,100r.median_threshold,100r.low_threshold,100r.high_threshold,r.median_switches,r.return_count));end
  println(io,"\nMean weights average valid sequences equally. Raster rows retain replicate order; pale versions of family colours indicate a runner-up within 2 AIC. Species bands show central 95% removal-sequence ranges, not confidence intervals. Small samples can produce unstable selections. Existing sampling, site-universe, occurrence and shared-collection limitations remain. Analysis is conditional on observed networks, not independent landscape validation.")
 end
 show(stdout,ts;allcols=true);println()
end
function plots(st,avg,tr,ts)
 set_theme!(Theme(fontsize=17,Axis=(topspinevisible=false,rightspinevisible=false)))
 f=Figure(size=(1550,750));h=Figure(size=(1550,750));s=Figure(size=(1550,750))
 cmap=vcat([RGBf(.8,.8,.8)],parse.(RGBf,COLORS))
 for (i,ds) in enumerate(DATA)
  row=div(i-1,5)+1;col=mod(i-1,5)+1;g=st[st.dataset.==ds,:];a=avg[avg.dataset.==ds,:]
  ax=Axis(f[row,col],title=replace(ds,"_"=>" "),xlabel=row==2 ? "Sites removed (%)" : "",ylabel=col==1 ? "Mean AIC weight" : "",xticks=0:25:100)
  for j in 1:4;q=a[a.model.==NAMES[j],:];lines!(ax,100q.removal,q.mean_weight,color=COLORS[j],linewidth=2.5);end
  xlims!(ax,0,100);ylims!(ax,0,1)
  ax2=Axis(h[row,col],title=replace(ds,"_"=>" "),xlabel=row==2 ? "Sites removed (%)" : "",ylabel=col==1 ? "Removal sequence" : "",xticks=0:25:100,yreversed=true)
  levels=sort(unique(g.removal));edges=vcat(0.,(levels[1:end-1].+levels[2:end])./2,min(1.,levels[end]+(levels[end]-levels[end-1])/2))
  matrix=fill(RGBAf(.8,.8,.8,1),length(levels),100)
  for r in eachrow(g)
   j=findfirst(==(r.removal),levels);c=cmap[r.winner+1]
   matrix[j,r.replicate]=r.winner>0 && r.delta_AIC_runner_up<2 ? RGBAf(.45c.r+.55,.45c.g+.55,.45c.b+.55,1) : RGBAf(c,1)
  end
  heatmap!(ax2,100edges,.5:1:100.5,matrix);xlims!(ax2,0,100)
  ax3=Axis(s[row,col],title=replace(ds,"_"=>" "),xlabel=row==2 ? "Sites removed (%)" : "",ylabel=col==1 ? "Connected consumers" : "",xticks=0:25:100)
  q=a[a.model.==NAMES[1],:];band!(ax3,100q.removal,q.low_species,q.high_species,color=(:gray,.25));lines!(ax3,100q.removal,q.mean_species,color=:black,linewidth=2);xlims!(ax3,0,100)
 end
 for fig in [f,h];Legend(fig[3,1:5],[PolyElement(color=c) for c in COLORS],NAMES,orientation=:horizontal,framevisible=false);end
 save(joinpath(OUT,"01_AIC_support.png"),f,px_per_unit=2);save(joinpath(OUT,"02_sequence_transitions.png"),h,px_per_unit=2);save(joinpath(OUT,"03_remaining_species.png"),s,px_per_unit=2)
 t=Figure(size=(1050,650));ax=Axis(t[1,1],xlabel="First sustained departure (% sites removed)",yticks=(1:10,replace.(DATA,"_"=>" ")),yreversed=true)
 for (i,ds) in enumerate(DATA)
  q=tr[(tr.dataset.==ds).&isfinite.(tr.threshold),:]
  scatter!(ax,100q.threshold,i.+range(-.22,.22,length=nrow(q)),color=(:gray,.3),markersize=6)
  r=ts[i,:]
  if isfinite(r.median_threshold);lines!(ax,[100r.low_threshold,100r.high_threshold],[i,i],color=COLORS[1],linewidth=5);scatter!(ax,[100r.median_threshold],[i],color=:black,markersize=11);end
  text!(ax,102,i,text="$(r.detected)/100",align=(:left,:center),fontsize=16)
 end
 xlims!(ax,0,116);save(joinpath(OUT,"04_transition_thresholds.png"),t,px_per_unit=2)
end
main()


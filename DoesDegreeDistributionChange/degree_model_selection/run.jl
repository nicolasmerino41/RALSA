using CSV,DataFrames,Statistics,LinearAlgebra,Random,SHA,CairoMakie,Printf
const ROOT=length(ARGS)>0 ? ARGS[1] : dirname(dirname(@__DIR__))
const OUT=joinpath(@__DIR__,"outputs")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const MODELS=["Power_law","Cutoff_power_law","Discrete_lognormal","Exponential"]
const NAMES=["Power law","Truncated power law","Log-normal","Exponential"]
const COLOURS=["#0072B2","#AA4499","#009E73","#D55E00"]
const REMOVALS=100
const B=199
const HELPER=joinpath(ROOT,"scripts/K_CO_SP_distributions/distribution_helpers.jl")
include(HELPER)
function degree_sequence(ids,keep,S)
 d=zeros(Int,S)
 for j in eachindex(ids);d[ids[j]]+=keep[j];end
 d[d.>0]
end
function parameter_text(model,fs)
 if model=="Power_law"
  @sprintf("α = %.2f",median(z.alpha for z in fs))
 elseif model=="Cutoff_power_law"
  @sprintf("α = %.2f, λ = %.3f",median(z.alpha for z in fs),median(z.lambda for z in fs))
 elseif model=="Exponential"
  @sprintf("λ = %.3f",median(z.lambda for z in fs))
 elseif count(z->z.boundary,fs)>length(fs)/2
  "Power-law boundary (σ → ∞)"
 else
  finite=filter(z->!z.boundary,fs)
  @sprintf("μ = %.2f, σ = %.2f",median(z.mu for z in finite),median(z.sigma for z in finite))
 end
end
function main()
 mkpath(OUT)
 manifest=CSV.read(joinpath(ROOT,"DoesDegreeDistributionChange/outputs/source_versions.csv"),DataFrame)
 for row in eachrow(manifest)
  @assert bytes2hex(sha256(read(joinpath(ROOT,row.path))))==row.sha256 "Input changed: $(row.path)"
 end
 CSV.write(joinpath(OUT,"input_versions.csv"),manifest)
 write(joinpath(OUT,"fitting_helper_version.txt"),HELPER*"\nSHA256: "*bytes2hex(sha256(read(HELPER))))
 pairs=CSV.read(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/pairs.csv"),DataFrame)
 pairs=pairs[pairs.K.>0,:]
 rng=MersenneTwister(20261028);bootrng=MersenneTwister(20261029)
 details=NamedTuple[];selected=NamedTuple[];scores=NamedTuple[]
 for (di,ds) in enumerate(DATA)
  p=pairs[pairs.dataset.==ds,:];N=only(unique(p.N));L=nrow(p)
  occ=CSV.read(joinpath(ROOT,"outputs/support_prediction_zero/confirmed/inputs",ds*"_occupancy.csv"),DataFrame;types=String)
  sites=sort(unique(occ.site));@assert length(sites)==N
  index=Dict(s=>i for (i,s) in enumerate(sites));links=Dict((r.consumer,r.resource)=>i for (i,r) in enumerate(eachrow(p)))
  mat=falses(N,L)
  ints=CSV.read(joinpath(ROOT,"data/processed",ds*"_interactions.csv"),DataFrame;types=String)
  for r in eachrow(ints);mat[index[r.site],links[(r.consumer,r.resource)]]=true;end
  @assert vec(sum(mat;dims=1))==p.K
  consumers=sort(unique(p.consumer));ids=[findfirst(==(s),consumers) for s in p.consumer]
  U=length(unique(p.resource));initial=degree_sequence(ids,trues(L),length(consumers))
  sequences=Dict(0.0 => [initial],0.4 => Vector{Int}[],0.8 => Vector{Int}[])
  for b in 1:REMOVALS
   order=randperm(rng,N)
   for f in [.4,.8]
    m=round(Int,f*N);keep=vec(any(mat[order[m+1:end],:];dims=1))
    d=degree_sequence(ids,keep,length(consumers));@assert !isempty(d) && maximum(d)<=maximum(initial)
    push!(sequences[f],d)
   end
  end
  for f in [0.,.4,.8]
   seqs=sequences[f];counts=[[count(==(k),d) for k in 1:U] for d in seqs]
   fitted=[[fit(c,m) for m in MODELS] for c in counts]
   weights=[begin
    a=[z.aic for z in fs];w=exp.(-0.5 .* (a .- minimum(a)));w./sum(w)
   end for fs in fitted]
   avgweights=[mean(w[mi] for w in weights) for mi in eachindex(MODELS)]
   winner=argmax(avgweights);ordered=sort(avgweights;rev=true)
   winshare=mean(argmax(w)==winner for w in weights)
   fs=[z[winner] for z in fitted];Ds=[z.D for z in fs]
   boots=[bootstrap(counts[i],MODELS[winner],fs[i],bootrng) for i in eachindex(counts)]
   ps=[z.p for z in boots]
   @assert maximum(z.kkt for z in fs)<1e-6
   for (b,allfits) in enumerate(fitted),(mi,z) in enumerate(allfits)
    push!(details,(community="C$di",dataset=ds,removal=f,actual_removal=round(Int,f*N)/N,replicate=b,n=sum(counts[b]),upper_degree=U,model=NAMES[mi],log_likelihood=z.ll,AIC=z.aic,AIC_weight=weights[b][mi],KS_D=z.D,KS_p=mi==winner ? ps[b] : NaN,bootstrap_replicates=mi==winner ? B : 0,alpha=z.alpha,lambda=z.lambda,mu=z.mu,sigma=z.sigma,boundary=z.boundary,kkt=z.kkt,bootstrap_max_kkt=mi==winner ? boots[b].maxkkt : NaN))
   end
   for mi in eachindex(MODELS)
    push!(scores,(community="C$di",dataset=ds,removal=f,model=NAMES[mi],mean_AIC_weight=avgweights[mi],replicate_win_fraction=mean(argmax(w)==mi for w in weights)))
   end
   push!(selected,(community="C$di",dataset=ds,removal=f,actual_removal=round(Int,f*N)/N,best_model=NAMES[winner],parameters=parameter_text(MODELS[winner],fs),KS_D=median(Ds),KS_p=median(ps),KS_p_low=minimum(ps),KS_p_high=maximum(ps),fraction_p_below_005=mean(ps.<.05),mean_AIC_weight=avgweights[winner],runner_up_weight=ordered[2],replicate_win_fraction=winshare,boundary_fraction=mean(z.boundary for z in fs),sequences=length(seqs)))
   println(ds," ",Int(100f),"%: ",NAMES[winner],"; weight=",round(avgweights[winner];digits=2));flush(stdout)
  end
 end
 det=DataFrame(details);sel=DataFrame(selected);sc=DataFrame(scores)
 CSV.write(joinpath(OUT,"fit_diagnostics.csv"),det);CSV.write(joinpath(OUT,"selected_results.csv"),sel);CSV.write(joinpath(OUT,"model_weights.csv"),sc)
 proportions=DataFrame([(removal=f,model=NAMES[mi],communities=count((sel.removal.==f).&(sel.best_model.==NAMES[mi])),percentage=10count((sel.removal.==f).&(sel.best_model.==NAMES[mi]))) for f in [0.,.4,.8] for mi in eachindex(MODELS)])
 CSV.write(joinpath(OUT,"summary_proportions.csv"),proportions)
 set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false)))
 overview=Figure(size=(1000,650))
 Label(overview[0,1],"Distribution Model Selection Across Network States\n(N=10 Communities)",fontsize=23)
 ax=Axis(overview[1,1],xlabel="Site removal",ylabel="Communities (%)",xticks=(1:3,["0%","40%","80%"]),yticks=0:20:100)
 for (j,f) in enumerate([0.,.4,.8])
  base=0.
  for mi in eachindex(MODELS)
   height=only(proportions.percentage[(proportions.removal.==f).&(proportions.model.==NAMES[mi])])
   if height>0
    poly!(ax,Rect2f(j-.3,base,.6,height);color=COLOURS[mi])
    text!(ax,j,base+height/2;text="$(height)%",align=(:center,:center),color=:white,fontsize=20)
   end
   base+=height
  end
 end
 xlims!(ax,.4,3.6);ylims!(ax,0,100)
 Legend(overview[2,1],[PolyElement(color=c) for c in COLOURS],NAMES;orientation=:horizontal,nbanks=2,framevisible=false)
 save(joinpath(OUT,"01_model_selection_overview.png"),overview;px_per_unit=2)
 table=Figure(size=(1600,1500))
 Label(table[1,1],"Community",font=:bold,fontsize=21)
 for (j,f) in enumerate([0.,.4,.8]);Label(table[1,j+1],f==0 ? "0% (Baseline)" : "$(Int(100f))% Removal\n(medians; 100 removals)",font=:bold,fontsize=18);end
 for (i,ds) in enumerate(DATA)
  Label(table[i+1,1],"C$i\n"*replace(ds,"_"=>" "),fontsize=18,halign=:left)
  for (j,f) in enumerate([0.,.4,.8])
   r=only(eachrow(sel[(sel.dataset.==ds).&(sel.removal.==f),:]))
   txt=r.best_model*"\nParam: "*r.parameters*@sprintf("\nKS D: %.3f\nKS p: %.3f",r.KS_D,r.KS_p)
   Label(table[i+1,j+1],txt,fontsize=17)
  end
 end
 rowgap!(table.layout,16);colgap!(table.layout,32);resize_to_layout!(table)
 save(joinpath(OUT,"02_reference_table.png"),table;px_per_unit=2)
 md=["# Distribution model selection across network states","","Consumers; ten communities. At 40% and 80%, cells report median parameters, KS D and KS p for the selected family across 100 removals. Median p is a summary of separate tests, not a combined test.","","| Community | 0% (Baseline) | 40% Removal | 80% Removal |","| :--- | :--- | :--- | :--- |"]
 for (i,ds) in enumerate(DATA)
  cells=String[]
  for f in [0.,.4,.8]
   r=only(eachrow(sel[(sel.dataset.==ds).&(sel.removal.==f),:]))
   push!(cells,r.best_model*"<br>Param: "*r.parameters*@sprintf("<br>KS D: %.3f<br>KS p: %.3f",r.KS_D,r.KS_p))
  end
  push!(md,"| **C$i — $(replace(ds,"_"=>" "))** | "*join(cells," | ")*" |")
 end
 append!(md,["","## What changes?",""])
 for f in [.4,.8]
  changed=String[];repeatability=Float64[]
  for ds in DATA
   a=only(eachrow(sel[(sel.dataset.==ds).&(sel.removal.==0),:]));b=only(eachrow(sel[(sel.dataset.==ds).&(sel.removal.==f),:]))
   if a.best_model!=b.best_model;push!(changed,replace(ds,"_"=>" "));push!(repeatability,b.replicate_win_fraction);end
  end
  push!(md,"- At $(Int(100f))% removal, the best-supported family changes in $(length(changed))/10 communities"*(isempty(changed) ? "." : ": "*join(changed,", ")*". In those communities, the selected new family wins a median $(round(Int,100median(repeatability)))% of individual removals."))
 end
 append!(md,["","At 40% removal the family is unchanged in all ten communities, while its parameters change. In all seven communities initially described by an exponential or truncated power law, the fitted decline/cutoff rate λ increases. This supports a shorter degree scale within the same family at moderate removal. At 80%, label changes become more common, but do not establish a universal regime shift.","","A change in the winning label is evidence of a change in relative model support, not by itself a structural regime shift. Overlapping candidate shapes, sampling variation and model uncertainty remain. The overview gives each community equal weight and assigns it to its highest mean AIC-weight family; it does not display that uncertainty.","","## How the fits work","","1. For each community, calculate the original positive consumer degrees. Generate 100 uniform, fixed-count site removals at each level (nearest integer number of sites); a link survives whenever at least one supporting site remains. Use nested 40%/80% removals within each random ordering. Fit each resulting degree sequence separately, excluding disconnected species.","2. Fit four discrete distributions by maximum likelihood on degrees 1 through U, where U is the initial number of resource species with an observed link. This same biological upper bound is used at all removal levels. No fitted tail threshold and no pooled replication.","3. Unnormalised probabilities: power law k^(-α); truncated power law k^(-α) exp(-λk); log-normal k^(-1) exp(-(log(k)-μ)^2/(2σ^2)); exponential exp(-λk). Each is normalised by summing over k=1,…,U. α and λ are nonnegative. The log-normal is a discrete density-kernel model, matching the previous implementation. A zero-curvature boundary is a power-law limit, not a finite μ/σ estimate.","4. MLE finds the parameters that maximise Σ log P(observed degree). α is estimated from the individual degree observations, not from the slope of a plotted cumulative line. Larger α means faster decline; λ controls exponential decline or cutoff. μ and σ describe position and spread in log-degree space.","5. AIC = 2 × parameter count − 2 × maximised log likelihood. Power law and exponential use one parameter; the other models use two. Within each sequence, weight = exp(−ΔAIC/2) / sum exp(−ΔAIC/2). At baseline select the highest weight; after removal select the family with the highest arithmetic mean weight across the 100 sequences. All communities have equal weight in the overview.","6. KS D is the largest absolute difference between the empirical and fitted cumulative probabilities. To calculate p, simulate 199 degree samples of the same size from the fitted family and refit each sample. p = (1 + number of simulated D values at least as large as the observed D) / 200. This handles discreteness and estimated parameters. It is a model-conditional diagnostic, not adjusted for choosing the winning family. Species degrees are dependent, so it is not a dependence-robust ecological test. A large p does not establish that a family is true.","7. The removal columns show medians of the selected family's separately fitted parameters and KS results. They are not fits to the averaged degree distribution. The full per-sequence diagnostics, weights, p-value ranges, failure fractions and optimisation checks are in the CSVs. No R² is used for fitting or selection.","","Inputs are hash-checked against the preceding analysis. Existing site-universe and occurrence limitations are unchanged; see DoesDegreeDistributionChange/outputs/outstanding_data_issues.csv. The same Göttingen sites and other collection dependencies limit independence between communities."])
 write(joinpath(OUT,"RESULTS_AND_METHODS.md"),join(md,"\n"))
 println("Complete. Maximum bootstrap optimisation residual: ",maximum(filter(isfinite,det.bootstrap_max_kkt)))
end
main()

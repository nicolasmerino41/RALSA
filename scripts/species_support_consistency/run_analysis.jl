module SpeciesSupportConsistency
using CSV,DataFrames,Statistics,Random,SHA,CairoMakie,Printf
const ROOT=isempty(ARGS) ? normpath(joinpath(@__DIR__,"../..")) : abspath(ARGS[1])
const SOURCE=joinpath(ROOT,"outputs/statistical description of local support/opportunity_theory")
const DEST=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(SOURCE,"species_support_consistency")
const DATA=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
const B=1999
mkpath(DEST)
dispersion(x,groups)=mean(2var(x[ix];corrected=true) for ix in groups)
function bh(p)
 order=sortperm(p);out=zeros(length(p));v=1.
 for j in length(p):-1:1;v=min(v,p[order[j]]*length(p)/j);out[order[j]]=v;end
 out
end
function run()
 versions=CSV.read(joinpath(SOURCE,"input_versions.csv"),DataFrame)
 for r in eachrow(versions);@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256;end
 qc=CSV.read(joinpath(SOURCE,"quality_control.csv"),DataFrame)
 rows=NamedTuple[];species=NamedTuple[];audit=NamedTuple[]
 for (di,ds) in enumerate(DATA)
  g=qc[qc.dataset.==ds,:];K=Float64.(g.K);C=Int.(g.C);q=K./C;N=only(unique(g.N));L=length(K)
  @assert all(1 .<=K .<=C .<=N) && all(g.valid)
  strata=[findall(==(c),C) for c in unique(C)]
  informative=sum(length(ix) for ix in strata if length(ix)>1 && length(unique(K[ix]))>1)
  for role in [:consumer,:resource]
   ids=String.(g[!,role]);allgroups=[findall(==(s),ids) for s in unique(ids)];groups=filter(ix->length(ix)>=2,allgroups)
   for ix in allgroups
    push!(species,(dataset=ds,role=String(role),species=ids[ix[1]],partners=length(ix),mean_K=mean(K[ix]),mean_C=mean(C[ix]),mean_q=mean(q[ix]),pooled_q=sum(K[ix])/sum(C[ix]),sd_q=length(ix)>1 ? std(q[ix]) : NaN))
   end
   isempty(groups) && error("No eligible species: $ds $role")
   observed=[dispersion(q,groups),dispersion(K,groups)]
   for (ni,null) in enumerate(["Unrestricted","Within_C"])
    rng=MersenneTwister(935000+100di+10ni+(role==:consumer ? 1 : 2));values=zeros(B,2)
    for b in 1:B
     perm=collect(1:L)
     if null=="Within_C"
      for ix in strata;perm[ix]=shuffle(rng,ix);end
      @assert C[perm]==C
     else;shuffle!(rng,perm);end
     for (j,x) in enumerate([q,K]);values[b,j]=dispersion(x[perm],groups);end
    end
    for (j,metric) in enumerate(["q","K"])
     v=values[:,j];mu=mean(v);lo,hi=quantile(v,[.025,.975]);p=(1+count(v.<=observed[j]))/(B+1)
     push!(rows,(dataset=ds,role=String(role),metric=metric,null=null,N=N,links=L,eligible_species=length(groups),eligible_links=sum(length,groups),movable_link_fraction=informative/L,observed_dispersion=observed[j],null_mean=mu,null_lower=lo,null_upper=hi,similarity_effect=mu>0 ? 1-observed[j]/mu : NaN,p_similarity=p,null_mcse=std(v)/sqrt(B)))
    end
   end
  end
  push!(audit,(dataset=ds,N=N,L=L,total_support=sum(K),opportunity_classes=length(strata),movable_links=informative,validated=true))
  println(ds," complete");flush(stdout)
 end
 t=DataFrame(rows);t.p_adjusted=zeros(nrow(t))
 for null in unique(t.null);ix=findall(==(null),t.null);t.p_adjusted[ix]=bh(t.p_similarity[ix]);end
 CSV.write(joinpath(DEST,"consistency_tests.csv"),t);CSV.write(joinpath(DEST,"species_summaries.csv"),DataFrame(species));CSV.write(joinpath(DEST,"checks.csv"),DataFrame(audit));CSV.write(joinpath(DEST,"input_versions.csv"),versions)
 paths=[joinpath(SOURCE,"quality_control.csv"),@__FILE__]
 CSV.write(joinpath(DEST,"source_versions.csv"),DataFrame(path=paths,sha256=[bytes2hex(sha256(read(p))) for p in paths]))
 for file in ["metadata_relationships.csv","metadata_source_versions.csv"];cp(joinpath(SOURCE,file),joinpath(DEST,file);force=true);end
end
function plots()
 t=CSV.read(joinpath(DEST,"consistency_tests.csv"),DataFrame)
 set_theme!(Theme(fontsize=18,Axis=(xgridvisible=false,ygridvisible=false,)))
 for (null,file) in [("Within_C","01_opportunity_controlled"),("Unrestricted","02_unadjusted")]
  fig=Figure(size=(1150,680))
  for (j,metric) in enumerate(["q","K"])
   ax=Axis(fig[1,j],xlabel="Within-species similarity beyond random allocation (%)",ylabel="",yticks=(1:10,replace.(DATA,"_"=>" ")),yreversed=true,title=metric=="q" ? "Realisation rate K/C" : "Absolute support K")
   vlines!(ax,[0];color=:gray,linestyle=:dash)
   for (ri,role) in enumerate(["consumer","resource"])
    col=ri==1 ? "#0072B2" : "#D55E00";offset=ri==1 ? -.13 : .13
    for (i,ds) in enumerate(DATA)
     r=only(eachrow(t[(t.dataset.==ds).&(t.role.==role).&(t.metric.==metric).&(t.null.==null),:]))
     lo=100*(1-r.null_upper/r.null_mean);hi=100*(1-r.null_lower/r.null_mean)
     lines!(ax,[lo,hi],[i+offset,i+offset];color=(col,.35),linewidth=5)
     scatter!(ax,[100r.similarity_effect],[i+offset];color=col,marker=r.p_adjusted<.05 ? :circle : :xcross,markersize=13)
    end
   end
  end
  Legend(fig[2,1:2],[MarkerElement(color="#0072B2",marker=:circle),MarkerElement(color="#D55E00",marker=:circle),MarkerElement(color=:black,marker=:circle),MarkerElement(color=:black,marker=:xcross)],["Consumers","Resources","Adjusted p < 0.05","Adjusted p ≥ 0.05"],orientation=:horizontal,framevisible=false)
  save(joinpath(DEST,file*".png"),fig;px_per_unit=2)
 end
 lines=["# Species consistency of local support","","Primary test: shuffle observed K only among links with exactly the same C, retaining pair identities and network topology. Each eligible species has at least two partners and receives equal weight. Statistic: average squared difference between two partners' q or K, averaged equally over species. Positive effect = smaller differences than the null; 100 × (1 − observed dispersion / mean null dispersion). 1,999 permutations; one-sided p for similarity; Benjamini–Hochberg correction across all 40 primary tests. Unrestricted tests have a separate 40-test correction.","","Plot segments show central 95% null intervals centred relative to the null mean; they are not confidence intervals around observed effects. Crosses denote lack of adjusted evidence for similarity, not evidence of no structure. K=1,C=1 links have q=1 by selection, and are handled by exact-C conditioning. Only observed regional links are tested.",""]
 for null in ["Within_C","Unrestricted"],metric in ["q","K"],role in ["consumer","resource"]
  v=t[(t.null.==null).&(t.metric.==metric).&(t.role.==role),:];sig=v.dataset[(v.p_adjusted.<.05).&(v.similarity_effect.>0)]
  push!(lines,"- $null, $role, $metric: $(length(sig))/10 datasets show adjusted evidence for similarity. Median effect $(round(100median(v.similarity_effect);digits=1))%. Significant datasets: $(isempty(sig) ? "none" : join(sig,", ")).")
 end
 append!(lines,["","The unrestricted comparison moves (K,C) together and asks whether species group similar observed supports/rates overall. Exact-C shuffling asks whether grouping exceeds what each species' opportunity profile already implies, while preserving the empirical support distribution at every C and total M. Sparse C classes cannot be shuffled meaningfully; movable fractions are reported. These are two marginal role tests, not simultaneous consumer/resource variance components: a species signal can reflect its partners, shared sites, observation effort or other network structure. It is not causal evidence for abundance or compatibility.","","Confirmed input versions and existing occurrence flags are retained. Uncertain excluded records: Gottin PP 474 missing-identity records, Garraf PP2 four zero columns, Montseny one zero row. Much occupancy is interaction-derived. Prior Gottin HP V1 correction, Garraf HP/Olot one-link discrepancies and undocumented sites remain limitations. Gottin datasets share sites and other collection overlaps remain unresolved; counts of significant datasets are descriptive, not independent replication probabilities."])
 write(joinpath(DEST,"FINDINGS.md"),join(lines,"\n"))
end
if !(length(ARGS)>2 && ARGS[3]=="plot");run();end
plots()
end

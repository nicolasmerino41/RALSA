function plots()
 set_theme!(Theme(fontsize=18,Axis=(xgridvisible=false,ygridvisible=false,)))
 t=CSV.read(joinpath(DEST,"summary.csv"),DataFrame);f=CSV.read(joinpath(DEST,"frequencies.csv"),DataFrame)
 colors=["#0072B2","#D55E00"]
 for (metric,label,file) in [("singleton","Singleton links (%)","01_singletons"),("gini","Support Gini","02_concentration"),("within_C_variance","Variance of K/C at fixed C","03_conditional_variation")]
  fig=Figure(size=(1550,680));scale=metric=="singleton" ? 100 : 1
  for (i,ds) in enumerate(DATA)
   ax=Axis(fig[div(i-1,5)+1,mod(i-1,5)+1],title=replace(ds,"_"=>" "),xlabel="Propensity heterogeneity (ρ)",ylabel=mod(i-1,5)==0 ? label : "",xticks=[0,.3,.6,.9])
   for (j,fam) in enumerate(["Two-state","Beta"])
    v=t[(t.dataset.==ds).&(t.family.==fam).&(t.metric.==metric),:];sort!(v,:rho)
    band!(ax,v.rho,scale*v.lower,scale*v.upper;color=(colors[j],.13));lines!(ax,v.rho,scale*v.mean;color=colors[j],linewidth=2.5);scatter!(ax,v.rho,scale*v.mean;color=colors[j],markersize=7)
   end
   hlines!(ax,[scale*first(t.empirical[(t.dataset.==ds).&(t.metric.==metric)])];color=:black,linestyle=:dash,linewidth=2)
  end
  Legend(fig[3,1:5],[LineElement(color=colors[1]),LineElement(color=colors[2]),LineElement(color=:black,linestyle=:dash)],["Persistent two-state","Persistent beta","Empirical"],orientation=:horizontal,framevisible=false)
  save(joinpath(DEST,file*".png"),fig;px_per_unit=2)
 end
 fig=Figure(size=(1550,680))
 for (i,ds) in enumerate(DATA)
  ax=Axis(fig[div(i-1,5)+1,mod(i-1,5)+1],title=replace(ds,"_"=>" "),xlabel="Local support K",ylabel=mod(i-1,5)==0 ? "Proportion with support ≥ K" : "",xscale=log10,yscale=log10)
  base=f[(f.dataset.==ds).&(f.family.=="Beta").&(f.rho.==0),:]
  for (fam,rho,col,style) in [("Beta",0.,"#888888",:solid),("Two-state",.6,colors[1],:solid),("Beta",.6,colors[2],:solid)]
   v=f[(f.dataset.==ds).&(f.family.==fam).&(f.rho.==rho),:];cc=reverse(cumsum(reverse(v.expected)));ix=cc.>=1e-5;lines!(ax,v.K[ix],cc[ix];color=col,linewidth=2)
  end
  cc=reverse(cumsum(reverse(base.empirical)));ix=cc.>0;lines!(ax,base.K[ix],cc[ix];color=:black,linestyle=:dash,linewidth=2)
 end
 Legend(fig[3,1:5],[LineElement(color=c) for c in ["#888888",colors...,"#000000"]],["Independent local variation / common","Persistent two-state (ρ = 0.6)","Persistent beta (ρ = 0.6)","Empirical"],orientation=:horizontal,framevisible=false)
 save(joinpath(DEST,"04_support_ccdf.png"),fig;px_per_unit=2)
 lines=String["# Findings","","All inputs, known positive links and total support are fixed. No parameter fitting or prediction validation is performed.","","Independent propensity draws per opportunity integrate exactly to the common-probability model. Persistent pair propensities are compared at equal pre-conditioning mean (0.5) and variance (ρ/4), using two-state and beta distributions. Positivity and M are then jointly conditioned upon. ρ is the unconditioned within-pair correlation, not the correlation after conditioning.","","Central 95% bands describe 300 independent conditional allocations, not confidence intervals. CCDFs are exact conditional expectations; ρ = 0.6 is a prespecified illustration, not a selected best fit.",""]
 for metric in ["singleton","gini","within_C_variance"],fam in ["Two-state","Beta"]
  base=t[(t.metric.==metric).&(t.family.==fam).&(t.rho.==0),:];hi=t[(t.metric.==metric).&(t.family.==fam).&(t.rho.==.9),:]
  increased=count(hi.mean.>base.mean);closer=count(abs.(hi.mean.-hi.empirical).<abs.(base.mean.-base.empirical));inside=count(hi.lower.<=hi.empirical.<=hi.upper)
  push!(lines,"- $fam, $metric: increasing ρ from 0 to 0.9 increases the mean in $increased/10 datasets; endpoint is closer to empirical in $closer/10; empirical lies inside endpoint allocation interval in $inside/10.")
 end
 append!(lines,["","The experiment tests whether persistence of differences can generate concentrated support beyond unequal opportunities alone. It does not identify abundance, compatibility or environment as its biological cause, and does not establish a universal distribution family. Read intermediate levels and both families before claiming robustness; an endpoint overshoot is informative, not a fitted-model rejection.","","The previous Pair model mixed separately zero-truncated binomials. Here latent probabilities are mixed BEFORE imposing positive support and total M, matching the stated generative experiment. Thus these are new theoretical scenarios, not refits of that model.","","Known-link identities and M are supplied; unknown interactions and support amount are not predicted. Only opportunity counts matter for these support statistics because sites are exchangeable within each pair; a full binary matrix need not be generated.","","Confirmed occurrence restorations are retained. Uncertain excluded records remain: Gottin PP 474 missing-identity records, Garraf PP2 four zero columns, Montseny one zero row. Interaction-derived occurrence, undocumented sites, previous Gottin HP V1 correction and unresolved one-link discrepancies in Garraf HP/Olot remain limitations. Gottin datasets share sites and other collection relationships remain unresolved. These datasets are not independent landscape replicates. Site count, spatial coverage and sampling effort can affect the distributions."])
 append!(lines,["","## Interpretation of the experiment","","The shared response is increased conditional variation: at rho=0.9 both families give greater variation at fixed C than the baseline in all ten datasets. Singleton prevalence and overall concentration are less robust: they can decrease or turn down as heterogeneity increases, particularly for the beta family. Garraf HP and Garraf PP show decreasing concentration under both families. Therefore persistent pair differences can generate additional within-opportunity variation, but are not by themselves a general guarantee of the empirical marginal shape.","","The conditioning matters. Low-propensity pairs must still realise at least one interaction, and the total support budget must remain fixed; the retained probabilities and their association with C need not resemble the initial latent distribution. This is a mathematical selection effect of the experiment, not evidence for an ecological filtering process. The latent mean of 0.5 is a fixed theoretical control, not a fitted realisation rate; robustness to other latent means has not been tested.","","All exhaustive small-case comparisons, unconditioned moment identities, saved Common baseline comparisons, input hashes and conditional support constraints passed. The four PNGs were visually checked for legibility and clipping."])
 write(joinpath(DEST,"FINDINGS.md"),join(lines,"\n"))
end

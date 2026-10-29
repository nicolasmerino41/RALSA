function deliver()
 sp=CSV.read(joinpath(OUT,"species_expected_retention.csv"),DataFrame)
 su=CSV.read(joinpath(OUT,"summary.csv"),DataFrame)
 ef=CSV.read(joinpath(OUT,"paired_effects.csv"),DataFrame)
 cc=CSV.read(joinpath(OUT,"scaled_ccdf.csv"),DataFrame)
 checks=CSV.read(joinpath(OUT,"checks.csv"),DataFrame)
 set_theme!(Theme(font="Arial",fontsize=18,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=14,yticklabelsize=14)))
 blue="#0072B2";orange="#D55E00";purple="#AA4499"
 for (h,role) in enumerate(["consumer","resource"])
  f=Figure(size=(1800,850))
  for (i,ds) in enumerate(DATA)
   ax=Axis(f[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Initial regional partners",ylabel=mod1(i,5)==1 ? "Expected partners retained (%)" : "",xscale=log10)
   for (fraction,col) in zip([.4,.8],[blue,purple])
    v=sp[(sp.dataset.==ds).&(sp.role.==role).&(sp.target_removal.==fraction),:]
    scatter!(ax,v.initial_degree,100v.expected_fraction;color=(col,.25),markersize=6)
    t=combine(groupby(v,:initial_degree),:expected_fraction=>mean=>:retained);sort!(t,:initial_degree)
    lines!(ax,t.initial_degree,100t.retained;color=col,linewidth=2.5)
    hlines!(ax,[100v.shuffle_expected_fraction[1]];color=col,linestyle=:dash,linewidth=1.5)
   end
   maxd=maximum(sp.initial_degree[(sp.dataset.==ds).&(sp.role.==role)]);ticks=filter(x->x<=maxd,[1,2,5,10,20,50,100,200]);ax.xticks=(ticks,string.(ticks));ylims!(ax,0,102)
  end
  Legend(f[3,1:5],[LineElement(color=blue,linewidth=3),LineElement(color=purple,linewidth=3),LineElement(color=:gray,linestyle=:dash)],["≈40% sites removed","≈80% sites removed","Random support assignment: expected baseline"],orientation=:horizontal,framevisible=false)
  save(joinpath(OUT,"0$(h)_$(role)_expected_retention.png"),f;px_per_unit=2)
 end
 f=Figure(size=(1500,1050));fe=Figure(size=(1500,1050));fd=Figure(size=(1500,1050))
 for (h,role) in enumerate(["consumer","resource"]),(j,frac) in enumerate([.4,.8])
  title=uppercasefirst(role)*"s · ≈$(Int(100frac))% removal"
  ax=Axis(f[h,j],title=title,xlabel="Change in degree inequality (Gini)",yticks=(1:10,replace.(DATA,"_"=>" ")),yreversed=true)
  ae=Axis(fe[h,j],title=title,xlabel="Observed − shuffled degree inequality",yticks=(1:10,replace.(DATA,"_"=>" ")),yreversed=true)
  ad=Axis(fd[h,j],title=title,xlabel="Species disconnected (%)",yticks=(1:10,replace.(DATA,"_"=>" ")),yreversed=true)
  vlines!(ax,[0.];color=:gray,linestyle=:dash);vlines!(ae,[0.];color=:gray,linestyle=:dash)
  for (i,ds) in enumerate(DATA)
   for (scenario,col,offset) in [("Observed",blue,-.12),("Shuffled",orange,.12)]
    r=only(eachrow(su[(su.dataset.==ds).&(su.role.==role).&(su.target_removal.==frac).&(su.scenario.==scenario),:]))
    lines!(ax,[r.change_lower,r.change_upper],[i+offset,i+offset];color=(col,.4),linewidth=3)
    scatter!(ax,[r.mean_change],[i+offset];color=col,markersize=10)
    scatter!(ad,[100r.mean_low_disconnected],[i+offset];color=col,marker=:circle,markersize=10)
    scatter!(ad,[100r.mean_high_disconnected],[i+offset];color=col,marker=:rect,markersize=10)
   end
   r=only(eachrow(ef[(ef.dataset.==ds).&(ef.role.==role).&(ef.target_removal.==frac),:]))
   lines!(ae,[r.paired_lower,r.paired_upper],[i,i];color=(:gray,.4),linewidth=3)
   lines!(ae,[r.mc_lower,r.mc_upper],[i,i];color=blue,linewidth=3)
   scatter!(ae,[r.mean_gini_observed_minus_shuffled],[i];color=blue,markersize=10)
  end
  xlims!(ax,-.32,.12);xlims!(ae,-.20,.20);xlims!(ad,0,100)
 end
 Legend(f[3,1:2],[MarkerElement(color=c,marker=:circle) for c in [blue,orange]],["Observed support assignment","Shuffled support assignment"],orientation=:horizontal,framevisible=false)
 Legend(fe[3,1:2],[LineElement(color=:gray,linewidth=3),LineElement(color=blue,linewidth=3)],["Central 95% of paired removal/allocation outcomes","95% Monte Carlo precision of the mean"],orientation=:horizontal,framevisible=false)
 Legend(fd[3,1:2],[MarkerElement(color=blue,marker=:circle),MarkerElement(color=orange,marker=:circle),MarkerElement(color=:gray,marker=:circle),MarkerElement(color=:gray,marker=:rect)],["Observed","Shuffled","Initially lower degree","Initially higher degree"],orientation=:horizontal,framevisible=false)
 save(joinpath(OUT,"03_degree_shape_change.png"),f;px_per_unit=2)
 save(joinpath(OUT,"04_support_assignment_effect.png"),fe;px_per_unit=2)
 save(joinpath(OUT,"05_species_disconnection.png"),fd;px_per_unit=2)
 for (h,role) in enumerate(["consumer","resource"]),(j,frac) in enumerate([.4,.8])
  f=Figure(size=(1800,850))
  for (i,ds) in enumerate(DATA)
   v=cc[(cc.dataset.==ds).&(cc.role.==role).&(cc.target_removal.==frac),:]
   obs=sort(v[v.scenario.=="Observed",:],:scaled_degree);nul=sort(v[v.scenario.=="Shuffled",:],:scaled_degree)
   ax=Axis(f[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" ")*"\n"*@sprintf("%.1f%% sites removed",100obs.actual_removal[1]),titlesize=17,xlabel="Degree / mean connected degree",ylabel=mod1(i,5)==1 ? "Connected species with degree ≥x" : "")
   lines!(ax,obs.scaled_degree,obs.initial_ccdf;color=:black,linewidth=2)
   lines!(ax,obs.scaled_degree,obs.ccdf;color=blue,linewidth=2.5)
   lines!(ax,nul.scaled_degree,nul.ccdf;color=orange,linewidth=2.5,linestyle=:dash)
   xlims!(ax,0,5);ylims!(ax,0,1)
  end
  Legend(f[3,1:5],[LineElement(color=:black),LineElement(color=blue),LineElement(color=orange,linestyle=:dash)],["Original network","Observed assignment after removal","Shuffled assignment after removal"],orientation=:horizontal,framevisible=false)
  save(joinpath(OUT,"$(5+2(h-1)+j)_$(role)_shape_$(Int(100frac)).png"),f;px_per_unit=2)
 end
 open(joinpath(OUT,"START_HERE.md"),"w") do io
  println(io,"# Does the degree distribution change?\n\nRead figures 01–02 for proportional retention, 03 for shape changes, and 04 for the effect of assigning support to particular species. Figure 05 separates species disconnection. Figures 6–9 compare mean-scaled cumulative distributions.\n")
  println(io,"## What was tested\n\nFor each dataset, consumers and resources are analysed separately. Exact expected fractions of partners retained use every link's K. Then ",B," paired replicates compare the observed site–link matrix with randomly reassigned complete columns. Both scenarios use the same removed sites, and the same assignment is used at both removal levels within each replicate. Removal levels are the nearest integer to 40% and 80%; exact fractions are reported.\n")
  println(io,"The shuffle preserves the regional graph, the entire K distribution, total support, every site's interaction count, and the number of surviving regional links for each removal sequence. It changes which species receive those support patterns, including their grouping within species. It is a structural null, not an occupancy-feasible alternative landscape. Thus differences cannot be attributed exclusively to the slope of mean K versus degree.\n")
  println(io,"## Results\n\n|Role|Target removal|Positive degree–retention correlation|Observed mean Gini increases / decreases|Observed mean Gini higher / lower than shuffle|\n|---|---:|---:|---:|---:|")
  for role in ["consumer","resource"],frac in [.4,.8]
   c=checks[(checks.role.==role).&(checks.target_removal.==frac),:];v=su[(su.role.==role).&(su.target_removal.==frac).&(su.scenario.=="Observed"),:];e=ef[(ef.role.==role).&(ef.target_removal.==frac),:]
   println(io,"|$role|$(Int(100frac))%|$(count(c.rho_initial_degree_expected_fraction.>0))/10|$(count(v.mean_change.>0)) / $(count(v.mean_change.<0))|$(count(e.mean_gini_observed_minus_shuffled.>0)) / $(count(e.mean_gini_observed_minus_shuffled.<0))|")
  end
  println(io,"\nCounts describe directions, not independent significant replications. Gini is the prespecified shape statistic: higher values indicate greater inequality among connected species, not necessarily a power-law exponent change or a heavier upper tail specifically. Mean-scaled curves show where distributions differ.\n")
  println(io,"|Dataset|Role|Actual removal|Degree–retention rho|Observed ΔGini|Shuffled ΔGini|Observed − shuffled|\n|---|---|---:|---:|---:|---:|---:|")
  for ds in DATA,role in ["consumer","resource"],frac in [.4,.8]
   c=only(eachrow(checks[(checks.dataset.==ds).&(checks.role.==role).&(checks.target_removal.==frac),:]));v=su[(su.dataset.==ds).&(su.role.==role).&(su.target_removal.==frac),:];a=only(v.mean_change[v.scenario.=="Observed"]);b=only(v.mean_change[v.scenario.=="Shuffled"])
   println(io,@sprintf("|%s|%s|%.1f%%|%.3f|%.3f|%.3f|%.3f|",ds,role,100c.actual_removal,c.rho_initial_degree_expected_fraction,a,b,a-b))
  end
  println(io,"\n## Reading the comparison\n\n- If initially high-degree species retain a greater fraction of links, that is more than merely retaining more links in absolute numbers.\n- If observed and shuffled shapes differ, assignment of support to the graph matters beyond the marginal support distribution.\n- Similar shapes can coexist with unequal species-level retention; distributional stability is not evidence of equal erosion.\n- Higher inequality is not synonymous with flattening every part of a cumulative curve. No distribution-family fits or exponent tests were selected after viewing results.\n")
  println(io,"## Definitions and uncertainty\n\nDegree zero is excluded from the connected-species shape calculation but retained in disconnection statistics and full raw histograms. Each replicate's degrees are divided by its own connected-species mean before averaging CCDFs. This is not the distribution of expected degrees. Original connected species define the denominator of disconnection. Low degree means initial degree <= the within-dataset/guild median, with ties kept together; high degree is above that median.\n\nFigure 03 intervals are central 95% removal/allocation ranges. Figure 04 grey intervals are central 95% paired outcome ranges; blue intervals are mean +/-1.96 Monte Carlo standard errors. Those show simulation precision only, not ecological confidence or independent-dataset inference. No permutation significance claim is made from one shuffle per removal sequence.\n")
  println(io,"Empty-remnant replicate count across scenarios: ",sum(su.empty_replicates),". Gini/CCDF condition on nonempty remnants; disconnection includes all replicates. CCDF display ends at five times the mean; unscaled integer-degree histograms include the full tail.\n")
  println(io,"## Input version and scope\n\nAll K values were reconstructed from current validated interactions and checked against the retained confirmed/pairs.csv; prior input hashes were verified. The site universe is reconstructed from confirmed occupancy and checked against the saved N for every dataset; the removed intermediate site-universe table is not required. Outstanding audit issues are copied, including Salix's 33 additional sampled zero-interaction localities and unincorporated Gottin_HP parasite slots. These results do not silently correct those inputs. The two Göttingen datasets share sites; other collection dependencies remain. Sampling and interaction-derived occupancy limitations still apply. No inference of intrinsic diet specialism or independent prediction of a new system is made.\n\nThis is conditional prediction from observed spatial support. The exact formula determines species' expected degrees; the full site patterns determine joint loss and realised distributions. A K distribution plus a correlation coefficient alone does not uniquely determine these outcomes.\n")
 end
end
deliver()


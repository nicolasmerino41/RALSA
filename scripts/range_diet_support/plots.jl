set_theme!(Theme(font="Arial",fontsize=18,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=14,yticklabelsize=14)))
colors=["#0072B2","#D55E00"];s=DataFrame(species);g=DataFrame(groups);o=DataFrame(origins);a=DataFrame(stats)
function legend_roles(f,row)
 Legend(f[row,1:5],[MarkerElement(color=c,marker=:circle) for c in colors],["Consumers","Resources"],orientation=:horizontal,framevisible=false)
end
for (role,col) in zip(["consumer","resource"],colors)
 f=Figure(size=(1800,850))
 for (i,ds) in enumerate(datasets)
  t=s[(s.dataset.==ds).&(s.role.==role).&(s.degree.>0),:]
  ax=Axis(f[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" ")*"\n"*@sprintf("Spearman ρ = %.2f",rho(t.degree,t.occupancy)),titlesize=17,xscale=log10,xlabel="Regional partners",ylabel=mod1(i,5)==1 ? "Sites occupied (%)" : "")
  ticks=filter(x->x<=maximum(t.degree),[1,2,5,10,20,50,100,200]);ax.xticks=(ticks,string.(ticks))
  scatter!(ax,t.degree,100t.occupancy_fraction;color=(col,.4),markersize=8)
  med=combine(groupby(t,:degree),:occupancy_fraction=>median=>:median);sort!(med,:degree)
  lines!(ax,med.degree,100med.median;color=col,linewidth=2);ylims!(ax,0,102)
 end
 save(joinpath(dest,role=="consumer" ? "01_consumer_degree_occupancy.png" : "02_resource_degree_occupancy.png"),f;px_per_unit=2)
end
f=Figure(size=(1800,850))
for (i,ds) in enumerate(datasets)
 ax=Axis(f[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Regional partner count",ylabel=mod1(i,5)==1 ? "Links with K = 1 (%)" : "",xticks=(1:5,["1","2","3–4","5–8","9+"]))
 for (role,col) in zip(["consumer","resource"],colors)
  v=g[(g.dataset.==ds).&(g.role.==role),:];scatterlines!(ax,v.class_index,100v.singleton_rate;color=col,markersize=9,linewidth=2)
 end
 ylims!(ax,0,102)
end
legend_roles(f,3);save(joinpath(dest,"03_singleton_rate_by_degree.png"),f;px_per_unit=2)
f=Figure(size=(1500,700))
for (j,role) in enumerate(["consumer","resource"])
 ax=Axis(f[1,j],title=uppercasefirst(role)*"s",xlabel="Percentage",yticks=(1:10,replace.(datasets,"_"=>" ")),yreversed=true)
 v=DataFrame(sensitivity);v=v[(v.role.==role).&(v.maximum_degree.==2),:]
 for (col,field,offset) in zip(["#999999","#0072B2","#D55E00"],[:species_fraction,:link_fraction,:singleton_share],[-.2,0.,.2])
  scatter!(ax,100v[!,field],collect(1:10).+offset;color=col,markersize=13)
 end
 xlims!(ax,0,100)
end
Legend(f[2,1:2],[MarkerElement(color=c,marker=:circle) for c in ["#999999","#0072B2","#D55E00"]],["Species with ≤2 partners","Links involving those species","Singleton links involving those species"],orientation=:horizontal,framevisible=false)
save(joinpath(dest,"04_specialists_species_vs_links.png"),f;px_per_unit=2)
f=Figure(size=(1700,700));cats=["Restricted endpoint","Restricted overlap","Beyond opportunity restriction"];cs=["#0072B2","#E69F00","#CC79A7"]
for (j,h) in enumerate([1,2,3])
 ax=Axis(f[1,j],title="Links with K ≤ $h",xlabel="Share of these links (%)",yticks=(1:10,replace.(datasets,"_"=>" ")),yreversed=true)
 for (i,ds) in enumerate(datasets)
  left=0.
  for (cat,col) in zip(cats,cs)
   w=100only(o.share_of_rare[(o.dataset.==ds).&(o.maximum_K.==h).&(o.category.==cat)])
   barplot!(ax,[i],[left+w];fillto=left,direction=:x,color=col,width=.7);left+=w
  end
 end
 xlims!(ax,0,100)
end
Legend(f[2,1:3],[PolyElement(color=c) for c in cs],["An endpoint occupies ≤h sites","Overlap C≤h; both endpoints occupy >h sites","C>h, but K≤h"],orientation=:horizontal,framevisible=false)
save(joinpath(dest,"05_rare_support_decomposition.png"),f;px_per_unit=2)
f=Figure(size=(1600,750))
for (j,metric) in enumerate(["degree_K_rho","low_degree_singleton_excess"])
 ax=Axis(f[1,j],xlabel=j==1 ? "Degree–K correlation minus null mean" : "Singleton-rate contrast minus null mean (pp)",yticks=(1:10,replace.(datasets,"_"=>" ")),yreversed=true)
 vlines!(ax,[0.];color=:gray,linestyle=:dash)
 for (role,col,offset) in zip(["consumer","resource"],colors,[-.12,.12])
  v=td[(td.role.==role).&(td.metric.==metric),:];scale=j==1 ? 1. : 100.
  for i in 1:10
   y=i+offset;lines!(ax,scale.*[v.null_lower[i]-v.null_mean[i],v.null_upper[i]-v.null_mean[i]],[y,y];color=(col,.5),linewidth=4)
   scatter!(ax,[scale*v.effect[i]],[y];color=col,marker=v.adjusted_p[i]<.05 ? :circle : :xcross,markersize=12)
  end
 end
end
Legend(f[2,1:2],[MarkerElement(color=c,marker=:circle) for c in colors],["Consumers","Resources"],orientation=:horizontal,framevisible=false)
save(joinpath(dest,"06_beyond_cooccurrence.png"),f;px_per_unit=2)
l=DataFrame(links);ccdfs=NamedTuple[]
for mode in ["occupancy","degree"]
 fig=Figure(size=(1800,850))
 for (i,ds) in enumerate(datasets)
  v=l[l.dataset.==ds,:]
  low=mode=="occupancy" ? min.(v.consumer_occupancy,v.resource_occupancy).<=2 : min.(v.consumer_degree,v.resource_degree).<=2
  ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Support K (sites)",ylabel=mod1(i,5)==1 ? "Links with support ≥K (%)" : "",xscale=log10)
  ticks=filter(x->x<=maximum(v.K)+1,[1,2,5,10,20,50,100,200]);ax.xticks=(ticks,string.(ticks))
  for (ix,label,col) in [(.!low,"Both endpoints above threshold","#0072B2"),(low,"At least one low endpoint","#D55E00")]
   any(ix) || continue
   k=collect(1:maximum(v.K)+1);p=[mean(v.K[ix].>=z) for z in k]
   lines!(ax,k,100p;color=col,linewidth=2.5,linestyle=col=="#D55E00" ? :dash : :solid)
   for (z,w) in zip(k,p);push!(ccdfs,(dataset=ds,split=mode,threshold=2,group=label,links=sum(ix),K=z,ccdf=w));end
  end
  ylims!(ax,0,102)
 end
 labels=mode=="occupancy" ? ["An endpoint occupies ≤2 sites","Both endpoints occupy >2 sites"] : ["An endpoint has ≤2 partners","Both endpoints have >2 partners"]
 Legend(fig[3,1:5],[LineElement(color=c,linewidth=3) for c in ["#D55E00","#0072B2"]],labels,orientation=:horizontal,framevisible=false)
 save(joinpath(dest,mode=="occupancy" ? "07_support_by_occupancy.png" : "08_support_by_partner_breadth.png"),fig;px_per_unit=2)
end
CSV.write(joinpath(dest,"conditional_support_ccdfs.csv"),DataFrame(ccdfs))

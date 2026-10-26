using CSV,DataFrames,CairoMakie,SHA
root=isempty(ARGS) ? "C:/Users/MM-1/OneDrive/PhD/GitHub/RALSA" : abspath(ARGS[1])
source=joinpath(root,"outputs/statistical description of local support/opportunity_theory/opportunity_support_organisation/distributions.csv")
dest=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(dirname(dirname(source)),"support_CO_ccdf")
mkpath(dest)
d=CSV.read(source,DataFrame)
datasets=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
colors=["#0072B2","#D55E00"]
set_theme!(Theme(font="Arial",fontsize=20,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=17,yticklabelsize=17)))
for logarithmic in [true,false]
 fig=Figure(size=(1800,850))
 for (i,ds) in enumerate(datasets)
  ax=Axis(fig[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "),xlabel="Number of sites",ylabel=mod1(i,5)==1 ? "Proportion of pairs with ≥ x sites" : "",xscale=logarithmic ? log10 : identity,yscale=logarithmic ? log10 : identity)
  maximum_x=1;minimum_y=1.
  for (j,quantity) in enumerate(["Support","Opportunity"])
   v=d[(d.dataset.==ds).&(d.quantity.==quantity),:];sort!(v,:K)
   @assert isapprox(sum(v.frequency),1;atol=1e-8)
   @assert maximum(abs.(v.ccdf.-reverse(cumsum(reverse(v.frequency)))))<1e-8
   keep=v.frequency.>0
   maximum_x=max(maximum_x,maximum(v.K[keep]));minimum_y=min(minimum_y,minimum(v.ccdf[keep]))
   lines!(ax,v.K[keep],v.ccdf[keep];color=colors[j],linewidth=3)
  end
  xlims!(ax,logarithmic ? .95 : 0,maximum_x*1.04)
  ylims!(ax,logarithmic ? minimum_y*.7 : 0,logarithmic ? 1.1 : 1.02)
  if logarithmic
   ticks=filter(<=(maximum_x),[1,2,5,10,20,50,100,200,500]);ax.xticks=(ticks,string.(ticks))
  end
 end
 Legend(fig[3,1:5],[LineElement(color=c,linewidth=3) for c in colors],["Local interaction support K","Co-occurrence opportunities C"],orientation=:horizontal,framevisible=false)
 save(joinpath(dest,logarithmic ? "01_support_CO_ccdf_log.png" : "02_support_CO_ccdf_linear.png"),fig;px_per_unit=2)
end
CSV.write(joinpath(dest,"curves.csv"),d[in.(d.quantity,Ref(["Support","Opportunity"])),:])
write(joinpath(dest,"README.md"),"# Support and co-occurrence distributions\n\nEmpirical complementary cumulative distributions P(K≥x) and P(C≥x) for the SAME known regional interaction pairs (K>0), separately for ten datasets. Co-occurring pairs never observed interacting are excluded. Log-log and linear-axis versions are provided. Lines join observed integer values, matching the requested historical degree-distribution style; no model fitting. Panel x ranges adapt to each dataset. A horizontal shift is a scale difference, not necessarily a shape difference. Since K≤C for every pair, the C curve must be at least as high at a common threshold. Existing confirmed-occurrence provenance and flags from opportunity_support_organisation apply unchanged.\n\nSource: "*relpath(source,root)*"\nSHA256: "*bytes2hex(sha256(read(source)))*"\n")

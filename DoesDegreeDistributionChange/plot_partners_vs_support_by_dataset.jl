using CSV, DataFrames, CairoMakie, SHA
root=length(ARGS)>0 ? ARGS[1] : @__DIR__
input=joinpath(root,"partners_vs_support","points.csv")
out=joinpath(@__DIR__,"partners_vs_support_by_dataset")
mkpath(out)
t=CSV.read(input,DataFrame)
datasets=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
colours=["#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#332288","#882255","#44AA99","#666666"]
set_theme!(Theme(font="Arial",fontsize=20,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
for role in ["consumer","resource","both"]
 f=Figure(size=(1800,850))
 for (i,ds) in enumerate(datasets)
  all=t[t.dataset.==ds,:]
  v=all[(all.role.==role).|(role=="both"),:]
  ax=Axis(f[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" "))
  scatter!(ax,v.regional_partners,v.K;color=(colours[i],.5),markersize=7)
  xlims!(ax,0,maximum(all.regional_partners)*1.05)
  ylims!(ax,0,maximum(all.K)*1.07)
 end
 Label(f[3,1:5],"Number of regional partners",fontsize=23)
 Label(f[1:2,0],"Sites supporting the interaction (K)",rotation=pi/2,fontsize=23)
 save(joinpath(out,role*"_partners_vs_support.png"),f;px_per_unit=2)
end
open(joinpath(out,"README.md"),"w") do io
 println(io,"# Regional partners versus local support, by dataset\n\nThree 2×5 figures: consumers, resources, and both combined. Each point is one species–interaction combination. X is the species' regional degree; Y is that interaction's K. No averaging or jitter: identical points overlap. The combined figure includes both endpoints of each regional interaction.\n\nLinear axes adapt to each dataset but remain identical for that dataset across the three figures. Dataset order and colours match the previous plot.\n\nSource: ",input,"\n\nSHA256: ",bytes2hex(sha256(read(input))),"\n\nThese are the same confirmed positive-support records as the previous plot; no inputs or analyses were changed. Existing data-quality flags remain applicable.")
end

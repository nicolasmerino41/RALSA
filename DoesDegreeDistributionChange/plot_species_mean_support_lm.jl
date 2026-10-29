using CSV, DataFrames, CairoMakie, SHA, Statistics, Printf
root=length(ARGS)>0 ? ARGS[1] : @__DIR__
input=joinpath(root,"partners_vs_support","points.csv")
out=joinpath(@__DIR__,"species_mean_support_lm")
mkpath(out)
links=CSV.read(input,DataFrame)
t=combine(groupby(links,[:dataset,:role,:species,:regional_partners]),:K=>mean=>:K,nrow=>:link_count)
@assert all(t.link_count.==t.regional_partners)
CSV.write(joinpath(out,"species_points.csv"),rename(copy(t),:K=>:mean_K))
datasets=["Quercus","Nahuel","Salix_Galpar","Gottin_HP","Gottin_PP","Garraf_HP","Garraf_PP","Garraf_PP2","Olot","Montseny"]
colours=["#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#332288","#882255","#44AA99","#666666"]
set_theme!(Theme(font="Arial",fontsize=20,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
fits=NamedTuple[]
for role in ["consumer","resource","both"]
 f=Figure(size=(1800,850))
 for (i,ds) in enumerate(datasets)
  all=t[t.dataset.==ds,:]
  v=all[(all.role.==role).|(role=="both"),:]
  x=Float64.(v.regional_partners);y=Float64.(v.K)
  slope=sum((x.-mean(x)).*(y.-mean(y)))/sum(abs2,x.-mean(x))
  intercept=mean(y)-slope*mean(x)
  predicted=intercept.+slope.*x
  r2=1-sum(abs2,y.-predicted)/sum(abs2,y.-mean(y))
  @assert isfinite(r2) && isapprox(mean(predicted),mean(y))
  push!(fits,(dataset=ds,role=role,n=nrow(v),intercept=intercept,slope=slope,r_squared=r2))
  equation=@sprintf("Mean K = %.2f %s %.3f d\nR² = %.3f",intercept,slope<0 ? "−" : "+",abs(slope),r2)
  ax=Axis(f[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" ")*"\n"*equation,titlesize=17)
  scatter!(ax,v.regional_partners,v.K;color=(colours[i],.5),markersize=7)
  line_x=[minimum(x),maximum(x)]
  lines!(ax,line_x,intercept.+slope.*line_x;color=:black,linewidth=2)
  xlims!(ax,0,maximum(all.regional_partners)*1.05)
  ylims!(ax,0,maximum(all.K)*1.07)
 end
 Label(f[3,1:5],"Number of regional partners",fontsize=23)
 Label(f[1:2,0],"Mean supporting sites per link (mean K)",rotation=pi/2,fontsize=23)
 save(joinpath(out,role*"_partners_vs_support.png"),f;px_per_unit=2)
end
CSV.write(joinpath(out,"linear_fits.csv"),DataFrame(fits))
open(joinpath(out,"README.md"),"w") do io
 println(io,"# Species mean local support versus regional degree\n\nEach point is one species within one dataset and guild. X is its number of regional partners, d. Y is the arithmetic mean of K across its observed regional links. Species sharing a degree remain separate points; identical coordinates overlap. The combined panel includes all consumer and resource species, retaining guild identity.\n\nEach species receives equal weight in an ordinary least-squares regression with intercept: mean K = intercept + slope × d. R² = 1 − residual sum of squares / total sum of squares. Lines cover the observed degree range. These are descriptive associations; no significance tests are reported. Species can share links.\n\nAxes adapt to each dataset but are identical across guild panels for that dataset. Full species points and coefficients are saved as CSVs.\n\nSource: ",input,"\n\nSHA256: ",bytes2hex(sha256(read(input))),"\n\nExisting data-quality flags remain applicable. Earlier figures are preserved.")
end

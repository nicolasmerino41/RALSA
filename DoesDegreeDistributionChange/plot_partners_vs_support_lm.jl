using CSV, DataFrames, CairoMakie, SHA, Statistics, Printf
root=length(ARGS)>0 ? ARGS[1] : @__DIR__
input=joinpath(root,"partners_vs_support","points.csv")
out=joinpath(@__DIR__,"partners_vs_support_lm")
mkpath(out)
t=CSV.read(input,DataFrame)
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
  equation=@sprintf("K = %.2f %s %.3f d   |   R² = %.3f",intercept,slope<0 ? "−" : "+",abs(slope),r2)
  ax=Axis(f[cld(i,5),mod1(i,5)],title=replace(ds,"_"=>" ")*"\n"*equation,titlesize=17)
  scatter!(ax,v.regional_partners,v.K;color=(colours[i],.5),markersize=7)
  line_x=[minimum(x),maximum(x)]
  lines!(ax,line_x,intercept.+slope.*line_x;color=:black,linewidth=2)
  xlims!(ax,0,maximum(all.regional_partners)*1.05)
  ylims!(ax,0,maximum(all.K)*1.07)
 end
 Label(f[3,1:5],"Number of regional partners",fontsize=23)
 Label(f[1:2,0],"Sites supporting the interaction (K)",rotation=pi/2,fontsize=23)
 save(joinpath(out,role*"_partners_vs_support.png"),f;px_per_unit=2)
end
CSV.write(joinpath(out,"linear_fits.csv"),DataFrame(fits))
open(joinpath(out,"README.md"),"w") do io
 println(io,"# Regional partners versus local support: linear fits\n\nSame points, colours, dataset order and axis limits as the preceding three figures. Each panel fits K = intercept + slope × d, where d is regional degree, by ordinary least squares with an intercept. R² = 1 − residual sum of squares / total sum of squares. Lines cover the observed degree range. Full coefficients are in linear_fits.csv.\n\nAll species–link rows enter the fit, including coincident points. No degree averaging: species with more partners contribute more rows. Both combined includes each link once for each endpoint. Fits describe associations; rows are not independent ecological replicates, and no significance tests are reported.\n\nSource: ",input,"\n\nSHA256: ",bytes2hex(sha256(read(input))),"\n\nExisting data-quality flags remain applicable. Earlier figures are preserved.")
end

using CSV,DataFrames,Statistics,LinearAlgebra,CairoMakie,SHA,Printf
const ROOT=dirname(@__DIR__)
const OUT=joinpath(@__DIR__,"outputs")
mkpath(OUT)
const S="outputs/statistical description of local support"
const DD="DoesDegreeDistributionChange"
const BLUE="#0072B2";const ORANGE="#D55E00";const GREEN="#009E73"
const FAMILY=["Power law","Truncated power law","Log-normal","Exponential"]
const FC=[BLUE,"#AA4499",GREEN,ORANGE]
const STAGES=["Species weighted","Co-occurrence","Local support"]
versions=NamedTuple[]
function load(p;kwargs...)
 full=joinpath(ROOT,p);push!(versions,(path=p,sha256=bytes2hex(sha256(read(full)))));CSV.read(full,DataFrame;kwargs...)
end
raw=load("$S/support_skewness/observations_and_weights.csv")
sk=load("$S/support_skewness/skewness_summary.csv")
points=load("$DD/species_mean_support_lm/species_points.csv")
lms=load("$DD/species_mean_support_lm/linear_fits.csv")
props=load("$DD/degree_model_selection/outputs/summary_proportions.csv")
orders=load("$DD/removal_transitions/outputs/removal_orders.csv";types=Dict(:site=>String))
oldcurves=load("$S/opportunity_theory/regional_local_retention_panels/curves.csv")
const DATA=unique(sk.dataset)
set_theme!(Theme(font="Arial",fontsize=20,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,titlesize=20)))
panel(cell,t)=Label(cell,t,font=:bold,fontsize=25,halign=:left,padding=(0,0,0,0))
savefig(f,name)=(save(joinpath(OUT,name*".png"),f,px_per_unit=2);save(joinpath(OUT,name*".svg"),f))
function survival(N,m,k)
 m<k && return 1.
 z=1.
 for i in 0:k-1;z*=(m-i)/(N-i);end
 1-z
end
function interp(x,y,g)
 [begin j=searchsortedlast(x,t);j==length(x) ? y[end] : y[j]+(y[j+1]-y[j])*(t-x[j])/(x[j+1]-x[j]);end for t in g]
end
# Figure 4: schematic retention scenarios and empirical support regressions.
f=Figure(size=(1850,800));l=GridLayout(f[1,1]);r=GridLayout(f[1,2]);colsize!(f.layout,1,Relative(.46));colgap!(f.layout,70)
ax=Axis(l[1,1],xlabel="Initial regional degree",ylabel="Relative partner retention",xticks=([0,1],["Low","High"]),yticks=([0,1],["Low","High"]),ylabelsize=20)
scenario_colors=[BLUE,"#AA4499",ORANGE]
scenario_names=["Specialists favoured","Neutral","Generalists favoured"]
for j in 1:3
 x=[0.,1.];y=j==1 ? [.8,.2] : j==2 ? [.5,.5] : [.2,.8]
 lines!(ax,x,y,color=scenario_colors[j],linewidth=3.5)
end
xlims!(ax,-.03,1.03);ylims!(ax,0,1)
axislegend(ax,[LineElement(color=z,linewidth=3.5) for z in scenario_colors],scenario_names,position=:ct,framevisible=false,labelsize=19)
panel(l[1,1,TopLeft()],"a")
for (i,ds) in enumerate(DATA)
 ax=Axis(r[div(i-1,5)+1,mod(i-1,5)+1],title=replace(ds,"_"=>" "),titlesize=15,xticklabelsize=13,yticklabelsize=13,xticks=WilkinsonTicks(3),yticks=WilkinsonTicks(3))
 v=points[points.dataset.==ds,:];z=only(eachrow(lms[(lms.dataset.==ds).&(lms.role.=="both"),:]))
 x=Float64.(v.regional_partners);y=v.mean_K;coef=hcat(ones(length(x)),x)\y
 @assert isapprox(coef[1],z.intercept;atol=1e-8) && isapprox(coef[2],z.slope;atol=1e-8)
 scatter!(ax,x,y,color=(BLUE,.4),markersize=4.5)
 xx=collect(extrema(x));lines!(ax,xx,z.intercept.+z.slope.*xx,color=ORANGE,linewidth=2.3)
end
Label(r[3,1:5],"Initial regional degree",fontsize=19);Label(r[1:2,0],"Mean local support per partner",rotation=pi/2,fontsize=19);panel(r[1,1,TopLeft()],"b")
savefig(f,"Figure4")


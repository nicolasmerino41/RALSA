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
curves=NamedTuple[];gaps=NamedTuple[]
for ds in DATA
 n=only(unique(sk.N[sk.dataset.==ds]));k=raw.value[(raw.dataset.==ds).&(raw.stage.=="Local support")];sp=raw.value[(raw.dataset.==ds).&(raw.stage.=="Species unweighted")]
 for m in 0:n
  r=mean(survival(n,m,x) for x in k);s=mean(survival(n,m,x) for x in sp)
  push!(curves,(dataset=ds,m=m,N=n,removal=m/n,regional=r,species=s,local_support=1-m/n,gap_pp=100(r-(1-m/n))))
 end
 c=DataFrame(curves);v=c[c.dataset.==ds,:];old=oldcurves[oldcurves.dataset.==ds,:]
 @assert length(old.m)==nrow(v) && maximum(abs.(old.regional_retention_percent.-100v.regional))<1e-7
 push!(gaps,(dataset=ds,N=n,mean_gap=mean(v.gap_pp[v.removal.<=.8])))
end
c=DataFrame(curves);gp=sort(DataFrame(gaps),:mean_gap);ordered=gp.dataset
CSV.write(joinpath(OUT,"figure2_curves.csv"),c);CSV.write(joinpath(OUT,"figure2_mean_gaps.csv"),gp)
# Figure 2: left half, two right quarters.
f=Figure(size=(1600,1050));left=GridLayout(f[1:2,1]);right=GridLayout(f[1:2,2]);colsize!(f.layout,1,Relative(.5));colgap!(f.layout,65);colsize!(left,1,Relative(1));colsize!(right,1,Relative(1));rowsize!(left,1,Relative(.91))
a=Axis(left[1,1],xlabel="Sites removed (%)",ylabel="Retained (%)",xticks=0:20:100,yticks=0:20:100)
grid=collect(0:.01:1)
for (col,color) in [(:species,GREEN),(:regional,BLUE),(:local_support,ORANGE)]
 ys=Vector{Float64}[]
 for ds in ordered
  v=c[c.dataset.==ds,:];lines!(a,100v.removal,100v[!,col],color=(color,.20),linewidth=1.5);push!(ys,interp(v.removal,v[!,col],grid))
 end
 lines!(a,100grid,100vec(mean(reduce(hcat,ys),dims=2)),color=color,linewidth=4)
end
xlims!(a,0,100);ylims!(a,0,101);panel(left[1,1,TopLeft()],"a")
Legend(left[2,1],[LineElement(color=z,linewidth=3) for z in [GREEN,BLUE,ORANGE]],["Species","Regional interactions","Local support"],framevisible=false,orientation=:vertical,halign=:left,tellwidth=false)
b=Axis(right[1,1],xlabel="Regional links (%)",yticks=(1:10,replace.(ordered,"_"=>" ")),yreversed=true,xticks=0:25:100,yticklabelsize=17)
bins=[(1,1),(2,2),(3,4),(5,8),(9,typemax(Int))];bc=["#D55E00","#E69F00","#56B4E9","#0072B2","#332288"]
histrows=NamedTuple[]
for (i,ds) in enumerate(ordered)
 k=raw.value[(raw.dataset.==ds).&(raw.stage.=="Local support")];base=0.
 for (j,(lo,hi)) in enumerate(bins)
  v=100mean(lo.<=k.<=hi);barplot!(b,[i],[v],fillto=base,direction=:x,color=bc[j],width=.7);base+=v
  push!(histrows,(dataset=ds,bin=["1","2","3–4","5–8","≥9"][j],percent=v))
 end
end
# Explicit rectangles avoid cumulative bar-height interpretation.
empty!(b)
for (i,ds) in enumerate(ordered)
 base=0.
 for (j,r) in enumerate(filter(r->r.dataset==ds,histrows))
  poly!(b,Rect2f(base,i-.35,r.percent,.7),color=bc[j]);base+=r.percent
 end
end
xlims!(b,0,100);ylims!(b,10.6,.4);panel(right[1,1,TopLeft()],"b")
Legend(right[2,1],[PolyElement(color=z) for z in bc],["K = 1","2","3–4","5–8","≥9"],orientation=:horizontal,framevisible=false,labelsize=15)
d=Axis(right[3,1],xlabel="Mean gap (percentage points)",yticks=(1:10,replace.(ordered,"_"=>" ")),yreversed=true,yticklabelsize=17)
barplot!(d,1:10,gp.mean_gap,direction=:x,color=BLUE,width=.65)
for i in 1:10;text!(d,gp.mean_gap[i]+.2,i,text=@sprintf("%.1f",gp.mean_gap[i]),align=(:left,:center),fontsize=16);end
xlims!(d,0,maximum(gp.mean_gap)*1.2);ylims!(d,10.6,.4);panel(right[3,1,TopLeft()],"c")
rowsize!(right,1,Relative(.43));rowsize!(right,3,Relative(.43));savefig(f,"Figure2")
CSV.write(joinpath(OUT,"figure2_support_bins.csv"),DataFrame(histrows))
# Figure 3: three frequency panels above; condensed skewness below left.
f=Figure(size=(1100,1100));top=GridLayout(f[1,1:2]);rowsize!(f.layout,1,Relative(.5));colsize!(f.layout,1,Relative(.5));example="Quercus"
for (j,stage) in enumerate(STAGES)
 v=raw[(raw.dataset.==example).&(raw.stage.==stage),:];n=maximum(v.value);freq=[sum(v.weight[v.value.==k])/sum(v.weight) for k in 1:22]
 ax=Axis(top[1,j],title=["Species occupancy","Co-occurrence","Local support"][j],xlabel="Sites",ylabel=j==1 ? "Proportion" : "",xticks=[1,10,20],yticks=0:.1:.5)
 barplot!(ax,1:22,freq,color=[GREEN,"#AA4499",ORANGE][j],gap=.1);xlims!(ax,.4,22.6);ylims!(ax,0,.55)
end
Label(top[0,1:3],"Quercus",font=:bold,fontsize=22);panel(top[1,1,TopLeft()],"b")
ax=Axis(f[2,1],xticks=(1:3,["Species","CO","Support"]),ylabel="Skewness",xticklabelsize=18)
for ds in DATA
 vals=[only(sk.skewness[(sk.dataset.==ds).&(sk.stage.==stage)]) for stage in STAGES]
 color=ds==example ? BLUE : :gray65
 lines!(ax,1:3,vals,color=(color,ds==example ? 1 : .7),linewidth=ds==example ? 4 : 1.5)
 scatter!(ax,1:3,vals,color=color,markersize=ds==example ? 13 : 6)
end
hlines!(ax,[0.],color=:gray85);xlims!(ax,.8,3.2);panel(f[2,1,TopLeft()],"c")
axislegend(ax,[LineElement(color=BLUE,linewidth=3),LineElement(color=:gray65,linewidth=2)],["Quercus","Other datasets"],framevisible=false,position=:lt,labelsize=16)
Box(f[2,2],color=:white,strokewidth=0);savefig(f,"Figure3_coded_panels")
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
# Figure 5: same 100 observed removal orders used for the model selection.
example="Garraf_PP2"
ints=load("data/processed/$(example)_interactions.csv";types=String)
ord=orders[orders.dataset.==example,:];sites=ord.site[ord.replicate.==1];N=length(sites)
pairs=unique(ints[:,[:consumer,:resource]]);cons=unique(pairs.consumer);maxdegree=maximum(count(==(s),pairs.consumer) for s in cons)
cc=NamedTuple[]
for frac in [0.,.4,.8]
 vals=Vector{Float64}[]
 for rep in 1:100
  order=ord[ord.replicate.==rep,:];sort!(order,:rank);kept=Set(order.site[round(Int,frac*N)+1:end])
  links=unique(ints[in.(ints.site,Ref(kept)),[:consumer,:resource]])
  deg=[count(==(s),links.consumer) for s in cons];filter!(>(0),deg)
  push!(vals,[mean(deg.>=k) for k in 1:maxdegree])
 end
 yy=vec(mean(reduce(hcat,vals),dims=2))
 for k in 1:maxdegree;push!(cc,(dataset=example,removal=frac,degree=k,ccdf=yy[k]));end
end
cd=DataFrame(cc);CSV.write(joinpath(OUT,"figure5_degree_curves.csv"),cd)
f=Figure(size=(1500,700));ax=Axis(f[1,1],title="Garraf PP2",xlabel="Regional partners",ylabel="P(degree ≥ x)",yscale=log10)
for (j,frac) in enumerate([0.,.4,.8])
 v=cd[cd.removal.==frac,:];keep=[i for i in 1:nrow(v) if v.ccdf[i]>0 && (i==1 || i==findlast(>(0),v.ccdf) || v.ccdf[i]!=v.ccdf[i-1])]
 lines!(ax,v.degree[keep],v.ccdf[keep],color=[BLUE,ORANGE,GREEN][j],linewidth=3)
end
panel(f[1,1,TopLeft()],"a");Legend(f[2,1],[LineElement(color=z,linewidth=3) for z in [BLUE,ORANGE,GREEN]],["0%","40%","80%"],orientation=:horizontal,framevisible=false,title="Sites removed")
ax2=Axis(f[1,2],xlabel="Sites removed (%)",ylabel="Datasets (%)",xticks=(1:3,["0","40","80"]),yticks=0:20:100)
for (i,frac) in enumerate([0.,.4,.8])
 base=0.
 for (j,name) in enumerate(FAMILY)
  h=only(props.percentage[(props.removal.==frac).&(props.model.==name)])
  if h>0
   poly!(ax2,Rect2f(i-.3,base,.6,h),color=FC[j]);text!(ax2,i,base+h/2,text="$(Int(h))%",align=(:center,:center),color=:white,fontsize=20)
  end
  base+=h
 end
end
xlims!(ax2,.4,3.6);ylims!(ax2,0,100);panel(f[1,2,TopLeft()],"b")
Legend(f[2,2],[PolyElement(color=z) for z in FC],FAMILY,nbanks=2,orientation=:horizontal,framevisible=false,labelsize=17)
colgap!(f.layout,70);savefig(f,"Figure5")
CSV.write(joinpath(OUT,"input_versions.csv"),DataFrame(versions))
println("All four figures saved to ",OUT)



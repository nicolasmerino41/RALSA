using CSV, DataFrames, CairoMakie, Statistics, Printf
include(joinpath(@__DIR__,"..","src","config.jl"))
const DEST=joinpath(OUT,"support_prediction_zero")
const CASES=["Original","Confirmed","Inclusive"]
const DIRS=[joinpath(OUT,"support_prediction"),joinpath(DEST,"confirmed"),joinpath(DEST,"inclusive")]
const COLORS=["#777777","#0072B2","#D55E00"]
tables=DataFrame[]
for (label,path) in zip(CASES,DIRS)
    d=CSV.read(joinpath(path,"metrics.csv"),DataFrame);d.scenario=fill(label,nrow(d));push!(tables,d)
end
m=vcat(tables...);CSV.write(joinpath(DEST,"comparison_metrics.csv"),m)
summary=combine(groupby(m,[:scenario,:scheme,:model]),:nlpd=>mean=>:nlpd,
    :total_relative_error=>(x->mean(abs.(x)))=>:total_error,
    :links_relative_error=>(x->mean(abs.(x)))=>:links_error,
    :gap_error=>(x->mean(abs.(x)))=>:gap_error)
CSV.write(joinpath(DEST,"comparison_summary.csv"),summary)
for d in tables[2:3]
    joined=innerjoin(select(tables[1],:scheme,:dataset,:model,:observed_total,:observed_links,:observed_gap),
        select(d,:scheme,:dataset,:model,:observed_total,:observed_links,:observed_gap);on=[:scheme,:dataset,:model],makeunique=true)
    @assert all(joined.observed_total.==joined.observed_total_1)
    @assert all(joined.observed_links.==joined.observed_links_1)
    @assert all(isapprox.(joined.observed_gap,joined.observed_gap_1;atol=1e-12))
end
set_theme!(Theme(font="Arial",fontsize=20,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false)))
function output(f,name)
    save(joinpath(DEST,name*".png"),f;px_per_unit=2);save(joinpath(DEST,name*".pdf"),f)
end
changes=CSV.read(joinpath(DEST,"input_changes.csv"),DataFrame)
f=Figure(size=(1400,650))
for (j,metric,label) in [(1,:opportunity,"Co-occurrence opportunities added (%)"),(2,:probability,"Empirical realisation, K / n")]
    ax=Axis(f[1,j],yticks=(1:10,[DISPLAY[d] for d in DATA]),yreversed=true,xlabel=label)
    for (s,c,off) in [("confirmed",COLORS[2],-.12),("inclusive",COLORS[3],.12)]
        d=changes[changes.scenario.==s,:]
        v=[only(d[d.dataset.==ds,metric==:opportunity ? :new_opportunities : :new_p_emp]) for ds in DATA]
        old=[only(d[d.dataset.==ds,metric==:opportunity ? :old_opportunities : :old_p_emp]) for ds in DATA]
        val=metric==:opportunity ? 100 .* (v./old.-1) : v
        for i in 1:10;lines!(ax,metric==:opportunity ? [0,val[i]] : [old[i],val[i]],[i+off,i+off],color=c);end
        scatter!(ax,val,(1:10).+off,color=c,markersize=12)
        if metric==:probability && s=="confirmed";scatter!(ax,old,1:10,color=COLORS[1],marker=:diamond,markersize=10);end
    end
end
Legend(f[2,1:2],[MarkerElement(color=c,marker=i==1 ? :diamond : :circle) for (i,c) in enumerate(COLORS)],CASES;orientation=:horizontal,framevisible=false)
output(f,"comparison_01_inputs")
f=Figure(size=(1450,850))
for (j,metric,label,factor) in [(1,:nlpd,"Mean negative log predictive probability",1),(2,:total_error,"Mean absolute support error (%)",100),
    (3,:links_error,"Mean absolute richness error (%)",100),(4,:gap_error,"Mean absolute gap error (percentage points)",100)]
    ax=Axis(f[cld(j,2),mod1(j,2)],xticks=(1:3,CASES),ylabel=label,ylabelsize=18)
    for (model,c) in [("M0","#D55E00"),("M1","#0072B2")]
        vals=[factor*only(summary[(summary.scenario.==s).&(summary.scheme.=="LODO").&(summary.model.==model),metric]) for s in CASES]
        scatterlines!(ax,1:3,vals,color=c,markersize=14,linewidth=2)
        for i in 1:3;text!(ax,i,vals[i],text=@sprintf("%.2f",vals[i]),offset=(0,10),align=(:center,:bottom),fontsize=16);end
    end
    ymax=factor*maximum(summary[summary.scheme.=="LODO",metric])
    ylims!(ax,0,1.18ymax)
end
Legend(f[3,1:2],[LineElement(color=c) for c in ["#D55E00","#0072B2"]],["M0: constant probability","M1: heterogeneous probabilities"];orientation=:horizontal,framevisible=false)
output(f,"comparison_02_transfer_accuracy")
f=Figure(size=(1800,900));lo=m[m.scheme.=="LODO",:]
cols=[(s,model) for s in CASES for model in ["M0","M1"]]
for (j,metric,label) in [(1,:total_relative_error,"Support error (%)"),(2,:links_relative_error,"Richness error (%)"),(3,:gap_error,"Gap error (percentage points)")]
    z=[100only(lo[(lo.scenario.==s).&(lo.model.==model).&(lo.dataset.==ds),metric]) for (s,model) in cols,ds in DATA]
    bound=maximum(abs.(z));ax=Axis(f[1,j],xticks=(1:6,[s*" "*v for (s,v) in cols]),xticklabelrotation=pi/3,yticks=(1:10,[DISPLAY[d] for d in DATA]),yreversed=true)
    hm=heatmap!(ax,1:6,1:10,z;colormap=[:steelblue4,:white,:darkorange3],colorrange=(-bound,bound))
    for a in 1:6,b in 1:10;text!(ax,a,b;text=@sprintf("%.1f",z[a,b]),align=(:center,:center),fontsize=14,color=abs(z[a,b])>.65bound ? :white : :black);end
    Colorbar(f[2,j],hm;vertical=false,label=label)
end
output(f,"comparison_03_dataset_errors")
println(summary)

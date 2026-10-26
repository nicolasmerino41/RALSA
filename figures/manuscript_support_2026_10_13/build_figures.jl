using CSV, DataFrames, CairoMakie, Statistics, SHA

const ROOT = length(ARGS) >= 1 ? abspath(ARGS[1]) : normpath(joinpath(@__DIR__, "../.."))
const OUT = length(ARGS) >= 2 ? abspath(ARGS[2]) : joinpath(ROOT, "outputs", "figures", "manuscript_support_2026_10_13")
mkpath(OUT)
const S = joinpath(ROOT, "outputs", "statistical description of local support")
const T = joinpath(S, "opportunity_theory")
const A = joinpath(T, "neutral_allocation")
const ORDER = ["Quercus", "Nahuel", "Salix_Galpar", "Gottin_HP", "Gottin_PP", "Garraf_HP", "Garraf_PP", "Garraf_PP2", "Olot", "Montseny"]
const LABELS = replace.(ORDER, "_" => " ")
const PURPLE = "#332288"
const TEAL = "#228877"
const BLUE = "#4477AA"
const ROSE = "#CC6677"
const CHECKS = DataFrame(check=String[], passed=Bool[], detail=String[])
const SOURCES = DataFrame(path=String[], sha256=String[])
function check(name, value, detail="")
    push!(CHECKS, (name, value, detail))
    value || error("Validation failed: $name $detail")
end
function readcsv(path; kwargs...)
    push!(SOURCES, (relpath(path, ROOT), bytes2hex(sha256(read(path)))))
    CSV.read(path, DataFrame; kwargs...)
end
near(a,b) = isapprox(a,b; atol=1e-8, rtol=1e-8)
for folder in (T,A)
    manifest = readcsv(joinpath(folder,"input_versions.csv"))
    for row in eachrow(manifest)
        check("input fingerprint", bytes2hex(sha256(read(joinpath(ROOT,row.path)))) == row.sha256, row.path)
    end
end
kdata = readcsv(joinpath(S,"link_support.csv"))
sites = readcsv(joinpath(S,"site_universe.csv"); types=Dict(:site=>String))
qc = readcsv(joinpath(T,"quality_control.csv"))
env = readcsv(joinpath(T,"envelope_summary.csv"))
ec = readcsv(joinpath(T,"envelope_curves.csv"))
ns = readcsv(joinpath(A,"summary.csv"))
draws = readcsv(joinpath(A,"draws.csv"))
pm = readcsv(joinpath(T,"prediction_metrics.csv"))
pc = readcsv(joinpath(T,"prediction_curves.csv"))
ps = readcsv(joinpath(T,"performance_summary.csv"))
fits = readcsv(joinpath(T,"fitted_parameters.csv"))
splits = readcsv(joinpath(T,"site_split.csv"))
output_manifest = readcsv(joinpath(T,"output_versions.csv"))
for row in eachrow(output_manifest)
    path=joinpath(T,row.file)
    if any(SOURCES.path .== relpath(path,ROOT))
        check("output fingerprint",bytes2hex(sha256(read(path)))==row.sha256,row.file)
    end
end

function expected_gap(K,N)
    [mean(1 - (m<k ? 0.0 : prod((m-j)/(N-j) for j in 0:k-1)) for k in K) - (1-m/N) for m in 0:N]
end
function interp(x,y,z)
    [begin
        i=clamp(searchsortedlast(x,t),1,length(x)-1)
        y[i]+(t-x[i])/(x[i+1]-x[i])*(y[i+1]-y[i])
    end for t in z]
end
bins = DataFrame(dataset=String[], N=Int[], bin=String[], fraction=Float64[])
single = DataFrame(dataset=String[], empirical_minus_null_pp=Float64[], lower_pp=Float64[], upper_pp=Float64[])
grid=collect(0:.01:.8)
retentions=zeros(length(grid),10)
for (di,d) in enumerate(ORDER)
    kd=kdata[kdata.dataset.==d,:]; q=qc[qc.dataset.==d,:]
    N=only(unique(kd.N)); K=kd.K
    check("site count $d",N==count(==(d),sites.dataset))
    check("link support matches QC $d",Dict((r.consumer,r.resource)=>r.K for r in eachrow(kd))==Dict((r.consumer,r.resource)=>r.K for r in eachrow(q)))
    check("valid opportunities $d",all(1 .<=q.K .<=q.C .<=q.N))
    occ=readcsv(joinpath(ROOT,"outputs","support_prediction_zero","confirmed","inputs",d*"_occupancy.csv"); types=String)
    ints=readcsv(joinpath(ROOT,"data","processed",d*"_interactions.csv"); types=String)
    universe=Set(sites.site[sites.dataset.==d])
    check("occupancy sites $d",issubset(Set(occ.site),universe))
    occupancy=Dict(sp=>Set(occ.site[occ.species.==sp]) for sp in unique(occ.species))
    for row in eachrow(q)
        c=length(intersect(get(occupancy,row.consumer,Set{String}()),get(occupancy,row.resource,Set{String}())))
        check("confirmed C $d $(row.consumer) $(row.resource)",c==row.C)
        support=length(unique(ints.site[(ints.consumer.==row.consumer).&(ints.resource.==row.resource)]))
        check("interaction K $d $(row.consumer) $(row.resource)",support==row.K)
    end
    for (name,mask) in zip(["1","2","3–4","5–8","≥9"],[K.==1,K.==2,(3 .<=K .<=4),(5 .<=K .<=8),K.>=9])
        push!(bins,(d,N,name,mean(mask)))
    end
    e=sort(ec[(ec.dataset.==d).&(ec.envelope.=="A"),:],:m)
    check("integer grid $d",e.m==collect(0:N))
    check("exact empirical curve $d",all(near.(e.empirical,expected_gap(K,N))))
    for en in ["A","B"]
        z=ec[(ec.dataset.==d).&(ec.envelope.==en),:]
        h=only(eachrow(env[(env.dataset.==d).&(env.envelope.==en),:]))
        mask=z.removal_fraction.<=.8
        for (col,sumcol) in [(:empirical,:integrated_empirical_pp),(:minimum,:integrated_minimum_pp),(:maximum,:integrated_maximum_pp)]
            check("mean $col envelope $en $d",near(100mean(z[mask,col]),h[sumcol]))
        end
    end
    retentions[:,di]=interp(e.removal_fraction,100 .* (e.empirical.+1 .-e.removal_fraction),grid)
    n=only(eachrow(ns[ns.dataset.==d,:])); dr=draws[draws.dataset.==d,:]
    check("saved allocation mean $d",near(mean(dr.gap_pp),n.null_mean_gap_pp))
    check("saved allocation quantiles $d",all(near.(quantile(dr.gap_pp,[.025,.975]),[n.null_lower_gap_pp,n.null_upper_gap_pp])))
    check("singleton mean $d",near(mean(dr.singleton_fraction),n.null_singleton))
    check("empirical singleton $d",near(mean(K.==1),n.empirical_singleton))
    check("empirical gap null $d",near(n.empirical_gap_pp,100mean(e.empirical[e.removal_fraction.<=.8])))
    lo,hi=100 .* (quantile(dr.singleton_fraction,[.025,.975]).-n.null_singleton)
    push!(single,(d,100(n.empirical_singleton-n.null_singleton),lo,hi))
    for scope in ["LODO","Within"]
        counts=Int[]
        for model in ["Common","Pair"]
            fit=only(eachrow(fits[(fits.dataset.==d).&(fits.scope.==scope).&(fits.model.==model),:]))
            training=Set(split(fit.training_datasets,";"))
            check("training membership $d $scope $model",training==Set(scope=="LODO" ? filter(!=(d),ORDER) : [d]))
            p=sort(pc[(pc.dataset.==d).&(pc.scope.==scope).&(pc.model.==model),:],:m)
            met=only(eachrow(pm[(pm.dataset.==d).&(pm.scope.==scope).&(pm.model.==model),:]))
            check("prediction grid $d $scope $model",p.m==e.m)
            check("prediction reference $d $scope $model",all(near.(p.empirical,e.empirical)))
            err=100 .* (p.predicted.-p.empirical)[p.removal_fraction.<=.8]
            check("prediction MAE $d $scope $model",near(mean(abs.(err)),met.mae_pp))
            check("prediction max $d $scope $model",near(maximum(abs.(err)),met.maximum_error_pp))
            check("prediction bias $d $scope $model",near(mean(err),met.signed_error_pp))
            push!(counts,met.evaluated_cells)
            check("saved prediction draws $d $scope $model",met.prediction_draws==400)
        end
        check("comparable cells $d $scope",counts[1]==counts[2])
    end
    splitrows=splits[splits.dataset.==d,:]
    check("site split universe $d",sort(splitrows.site_index)==collect(1:N))
end
displaycurve=DataFrame(removal_fraction=grid,mean_regional_percent=vec(mean(retentions;dims=2)),minimum_percent=vec(minimum(retentions;dims=2)),maximum_percent=vec(maximum(retentions;dims=2)),local_percent=100 .* (1 .-grid))
for scope in ["LODO","Within"],model in ["Common","Pair"]
    value=mean(pm.mae_pp[(pm.scope.==scope).&(pm.model.==model)])
    saved=only(ps.mean_mae_pp[(ps.scope.==scope).&(ps.model.==model)])
    check("equal-dataset mean $scope $model",near(value,saved))
end
CSV.write(joinpath(OUT,"validation_checks.csv"),CHECKS)
CSV.write(joinpath(OUT,"source_versions.csv"),unique(SOURCES))
CSV.write(joinpath(OUT,"figure2A_support_bins.csv"),bins)
CSV.write(joinpath(OUT,"figure2B_display_interpolation.csv"),displaycurve)
CSV.write(joinpath(OUT,"figure2C_mean_gaps.csv"),env[env.envelope.=="A",[:dataset,:N,:integrated_empirical_pp]])
CSV.write(joinpath(OUT,"figure3B_singleton_departures.csv"),single)
CSV.write(joinpath(OUT,"figure4_plotted_metrics.csv"),pm[in.(pm.model,Ref(["Common","Pair"])),:])

set_theme!(Theme(font="Arial",fontsize=21,Axis=(xticklabelsize=18,yticklabelsize=18,xgridvisible=false,ygridvisible=false,topspinevisible=false,rightspinevisible=false),Legend=(framevisible=false,labelsize=18)))
letter(pos,s)=Label(pos,s;fontsize=28,font=:bold,halign=:left,tellwidth=false,padding=(0,0,0,5))
exportfig(f,name)=save(joinpath(OUT,name*".png"),f;px_per_unit=3)
function rowsaxis(pos; xlabel="", labels=LABELS)
    ax=Axis(pos;xlabel,yticks=(1:length(labels),labels),yreversed=true)
    ylims!(ax,.4,length(labels)+.6)
    ax.yreversed=true
    ax
end

# All example panels use these same matrices.
O=[1 1 1 0 0;1 1 1 1 1;1 0 1 1 1;0 1 0 1 1;1 0 0 1 1]
I=[1 0 0 0 0;0 1 1 0 1;1 0 0 0 1;0 1 0 0 1;0 0 0 0 1]
K=vec(sum(I;dims=2)); C=vec(sum(O;dims=2)); kept=copy(I);kept[:,1:2].=0
@assert all(I.<=O) && all(K.>0) && sum(kept)==5 && count(>(0),vec(sum(kept;dims=2)))==4
CSV.write(joinpath(OUT,"figure1_example.csv"),DataFrame(link=repeat(1:5,inner=5),site=repeat(1:5,5),opportunity=vec(permutedims(O)),interaction=vec(permutedims(I)),retained_interaction=vec(permutedims(kept))))
function matrixplot(pos,mat,counts,kind; removed=false)
    ax=Axis(pos;xticks=(1:5,string.(1:5)),yticks=(1:5,["Link $i" for i in 1:5]),xlabel=kind,yreversed=true)
    hidespines!(ax);ax.xticksvisible=false;ax.yticksvisible=false
    for i in 1:5,j in 1:5
        gone=removed && j<=2
        color=gone ? "#F2E4DC" : mat[i,j]==1 ? (kind=="Opportunities" ? "#CAD5DC" : "#B7D9CE") : "#F4F5F6"
        poly!(ax,Rect2f(j-.46,i-.46,.92,.92);color,strokewidth=0)
        text!(ax,j,i;text=gone ? "×" : string(mat[i,j]),align=(:center,:center),fontsize=20,color=gone ? "#AA6633" : "#223333")
    end
    for i in 1:5
        text!(ax,6,i;text=string(counts[i]),align=(:center,:center),font=:bold)
    end
    text!(ax,6,.1;text=kind=="Opportunities" ? "C" : removed ? "K′" : "K",align=(:center,:center),font=:bold)
    xlims!(ax,.45,6.5);ylims!(ax,-.25,5.6);ax.yreversed=true
    ax
end
f=Figure(size=(1520,890));left=f[1:2,1]=GridLayout();right=f[1:2,2]=GridLayout()
letter(left[1,1:2],"A")
matrixplot(left[2,1],O,C,"Opportunities")
matrixplot(left[2,2],I,K,"Interactions")
letter(left[3,1:2],"B")
matrixplot(left[4,1],kept,vec(sum(kept;dims=2)),"Sites 1–2 removed";removed=true)
ax=Axis(left[4,2]);hidedecorations!(ax);hidespines!(ax);xlims!(ax,0,1);ylims!(ax,0,1)
text!(ax,.05,.8;text="Local occurrences",fontsize=21)
text!(ax,.05,.61;text="9 → 5",fontsize=34,color=TEAL)
text!(ax,.05,.34;text="Regional links",fontsize=21)
text!(ax,.05,.15;text="5 → 4",fontsize=34,color=PURPLE)
letter(right[1,1],"C")
ax=Axis(right[2,1];xlabel="Sites removed (%)",ylabel="Retention (%)",xticks=0:20:100,yticks=0:20:100)
x=collect(0:5)./5;localret=1 .-x;reg=expected_gap(K,5).+localret
band!(ax,100x,100localret,100reg;color=(PURPLE,.17),label="Regional–local retention gap")
lines!(ax,100x,100reg;color=TEAL,linewidth=3,label="Regional links")
scatter!(ax,100x,100reg;color=TEAL,markersize=9)
lines!(ax,100x,100localret;color=:black,linestyle=:dash,linewidth=2,label="Local support")
vlines!(ax,[80];color=:gray,linestyle=:dot,linewidth=1.5)
lines!(ax,[0,80],[5,5];color=:gray,linewidth=2)
text!(ax,40,8;text="0–80% summary range",align=(:center,:bottom),fontsize=17,color="#555555")
xlims!(ax,0,100);ylims!(ax,0,102);axislegend(ax;position=:rt,labelsize=17)
colsize!(f.layout,1,Relative(.53));exportfig(f,"figure1_conceptual")

f=Figure(size=(1870,760));for (j,l) in enumerate(["A","B","C"]);letter(f[1,j],l);end
Ns=[only(unique(kdata.N[kdata.dataset.==d])) for d in ORDER]
ax=rowsaxis(f[2,1];xlabel="Regional links (%)",labels=["$(LABELS[i]) (N = $(Ns[i]))" for i in 1:10])
colors=["#E0E7F3","#AFBFDC","#7896C4","#476AA3","#203D70"]
for (i,d) in enumerate(ORDER)
    values=bins.fraction[bins.dataset.==d];offset=0.0
    for j in 1:5
        width=100values[j]
        poly!(ax,Rect2f(offset,i-.31,width,.62);color=colors[j],strokecolor=:white,strokewidth=.7)
        offset+=width
    end
end
xlims!(ax,0,100)
ax=Axis(f[2,2];xlabel="Sites removed (%)",ylabel="Retention (%)",xticks=0:20:80,yticks=0:20:100)
band!(ax,100grid,displaycurve.minimum_percent,displaycurve.maximum_percent;color=(TEAL,.18),label="Dataset min–max")
lines!(ax,100grid,displaycurve.mean_regional_percent;color=TEAL,linewidth=3,label="Mean regional retention")
lines!(ax,100grid,displaycurve.local_percent;color=:black,linestyle=:dash,linewidth=2,label="Local support")
xlims!(ax,0,80);ylims!(ax,0,102);axislegend(ax;position=:lb,labelsize=17)
ax=rowsaxis(f[2,3];xlabel="Mean gap (percentage points)")
vals=[only(env.integrated_empirical_pp[(env.dataset.==d).&(env.envelope.=="A")]) for d in ORDER]
for i in 1:10;lines!(ax,[0,vals[i]],[i,i];color=(PURPLE,.25),linewidth=2);end
scatter!(ax,vals,1:10;color=PURPLE,markersize=12);xlims!(ax,0,21)
Legend(f[3,1:3],[PolyElement(color=c) for c in colors],["K = 1","K = 2","K = 3–4","K = 5–8","K ≥ 9"];orientation=:horizontal)
colgap!(f.layout,35);exportfig(f,"figure2_support_retention")

f=Figure(size=(1510,760));letter(f[1,1],"A");letter(f[1,2],"B")
ax=rowsaxis(f[2,1];xlabel="Mean gap (percentage points)")
for (i,d) in enumerate(ORDER)
    for (en,col,lw) in [("A","#D9D9D9",18),("B","#858585",9)]
        row=only(eachrow(env[(env.dataset.==d).&(env.envelope.==en),:]))
        lines!(ax,[row.integrated_minimum_pp,row.integrated_maximum_pp],[i,i];color=col,linewidth=lw)
    end
    n=only(eachrow(ns[ns.dataset.==d,:]))
    lines!(ax,[n.null_lower_gap_pp,n.null_upper_gap_pp],[i,i];color=TEAL,linewidth=4)
    scatter!(ax,[n.null_mean_gap_pp],[i];color=:white,strokecolor=TEAL,strokewidth=2,markersize=11)
    scatter!(ax,[n.empirical_gap_pp],[i];color=PURPLE,marker=:diamond,markersize=15)
end
xlims!(ax,0,35)
ax=rowsaxis(f[2,2];xlabel="Excess singleton links (percentage points)")
vlines!(ax,[0];color="#777777",linestyle=:dash,linewidth=1)
for (i,d) in enumerate(ORDER)
    r=only(eachrow(single[single.dataset.==d,:]))
    lines!(ax,[r.lower_pp,r.upper_pp],[i,i];color=TEAL,linewidth=4)
    scatter!(ax,[0],[i];color=:white,strokecolor=TEAL,strokewidth=2,markersize=11)
    scatter!(ax,[r.empirical_minus_null_pp],[i];color=PURPLE,marker=:diamond,markersize=15)
end
Legend(f[3,1:2],[LineElement(color="#D9D9D9",linewidth=14),LineElement(color="#858585",linewidth=8),[LineElement(color=TEAL,linewidth=3),MarkerElement(color=:white,strokecolor=TEAL,strokewidth=2,marker=:circle,markersize=10)],MarkerElement(color=PURPLE,marker=:diamond,markersize=13)],["N, L, M range","Opportunity-constrained range","Allocation 95% interval / mean","Empirical"];orientation=:horizontal,labelsize=17)
colgap!(f.layout,40);exportfig(f,"figure3_constraints_allocation")

f=Figure(size=(1470,790))
for (j,scope) in enumerate(["LODO","Within"])
    letter(f[1,j],j==1 ? "A" : "B")
    task=scope=="LODO" ? "Entire-dataset holdout" : "Within-system holdout"
    local ax=rowsaxis(f[2,j];xlabel=task*": gap MAE (pp)",labels=vcat(LABELS,["Equal-dataset mean"]))
    values=Dict(model=>[only(pm.mae_pp[(pm.dataset.==d).&(pm.scope.==scope).&(pm.model.==model)]) for d in ORDER] for model in ["Common","Pair"])
    for model in ["Common","Pair"];push!(values[model],mean(values[model]));end
    hlines!(ax,[10.5];color="#DDDDDD",linewidth=1)
    for i in 1:11
        lines!(ax,[values["Common"][i],values["Pair"][i]],[i,i];color="#BBBBBB",linewidth=i==11 ? 3 : 1.7)
    end
    for (model,col,marker) in [("Common",BLUE,:circle),("Pair",ROSE,:diamond)]
        scatter!(ax,values[model],1:11;color=col,marker,markersize=vcat(fill(12,10),18))
    end
    xlims!(ax,0,j==1 ? 12.5 : 2.7)
end
Legend(f[3,1:2],[MarkerElement(color=BLUE,marker=:circle,markersize=13),MarkerElement(color=ROSE,marker=:diamond,markersize=14)],["Common probability","Persistent pair"];orientation=:horizontal)
colgap!(f.layout,45);exportfig(f,"figure4_prediction")
println("Saved four figures; ",nrow(CHECKS)," checks passed. Output: ",OUT)

module GapReconstruction
using CSV, DataFrames, Statistics, SHA, Printf, CairoMakie
include(joinpath(@__DIR__,"..","src","config.jl"))
const SOURCE=joinpath(OUT,"statistical description of local support")
const FIT=joinpath(SOURCE,"two_parameter_fit")
const DEST=joinpath(SOURCE,"gap_reconstruction")
const PALETTE=["#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499","#661100"]
mkpath(DEST)
versions=CSV.read(joinpath(FIT,"input_versions.csv"),DataFrame)
descriptive=CSV.read(joinpath(FIT,"descriptive_input_versions.csv"),DataFrame)
for row in eachrow(vcat(select(versions,:path,:sha256),descriptive))
    @assert bytes2hex(sha256(read(joinpath(ROOT,row.path))))==row.sha256 "Input version changed"
end
CSV.write(joinpath(DEST,"input_versions.csv"),versions)
paths=[joinpath(SOURCE,"link_support.csv"),joinpath(SOURCE,"site_universe.csv"),joinpath(FIT,"fit_summary.csv"),joinpath(FIT,"fitted_distributions.csv")]
CSV.write(joinpath(DEST,"reconstruction_input_versions.csv"),DataFrame(path=relpath.(paths,ROOT),sha256=[bytes2hex(sha256(read(p))) for p in paths]))
links=CSV.read(paths[1],DataFrame);sites=CSV.read(paths[2],DataFrame;types=Dict(:site=>String))
parameters=CSV.read(paths[3],DataFrame);fitted=CSV.read(paths[4],DataFrame)
curves=NamedTuple[];summary=NamedTuple[]
for ds in DATA
    v=links[links.dataset.==ds,:];N=only(unique(v.N));K=Int.(v.K);L=length(K)
    @assert nrow(sites[sites.dataset.==ds,:])==N
    ints=CSV.read(joinpath(INPUT,ds*"_interactions.csv"),DataFrame;types=String)
    raw=combine(groupby(ints,[:consumer,:resource]),:site=>(x->length(unique(x)))=>:K)
    joined=innerjoin(raw,v;on=[:consumer,:resource],makeunique=true)
    @assert nrow(joined)==L==nrow(raw) && all(joined.K.==joined.K_1)
    par=only(eachrow(parameters[(parameters.dataset.==ds).&(parameters.model.=="two_parameter"),:]))
    fit=sort(fitted[(fitted.dataset.==ds).&(fitted.model.=="two_parameter"),:],:K)
    @assert fit.K==collect(1:N)
    p=fit.fitted_frequency;emp=[count(==(k),K)/L for k in 1:N]
    @assert isapprox(sum(p),1;atol=1e-12) && all(isapprox.(fit.observed_frequency,emp;atol=1e-12))
    # Reconstruct probabilities from saved parameters only to verify consistency.
    tail=par.r.^(0:N-2);tail./=sum(tail)
    @assert all(isapprox.(p,vcat(par.s,(1-par.s).*tail);atol=1e-12))
    for m in 0:N
        survival=ones(N);logext=0.
        for k in 1:m
            logext+=log((m-k+1)/(N-k+1))
            survival[k]=-expm1(logext)
        end
        re=sum(emp.*survival);rf=sum(p.*survival);local_retention=1-m/N
        @assert re>=local_retention-1e-12 && rf>=local_retention-1e-12
        push!(curves,(dataset=ds,N=N,removed_sites=m,removal_fraction=m/N,in_main_range=5m<=4N,
            empirical_regional_retention=re,fitted_regional_retention=rf,local_support_retention=local_retention,
            empirical_gap=re-local_retention,fitted_gap=rf-local_retention,error_pp=100(rf-re)))
    end
    d=DataFrame(filter(x->x.dataset==ds,curves));main=d[d.in_main_range,:];err=main.error_pp
    @assert isapprox(first(d.empirical_regional_retention),1) && last(d.empirical_regional_retention)==0
    @assert isapprox(first(d.fitted_regional_retention),1) && last(d.fitted_regional_retention)==0
    @assert all(diff(d.empirical_regional_retention).<=1e-12) && all(diff(d.fitted_regional_retention).<=1e-12)
    @assert isapprox(d.empirical_regional_retention[2],1-emp[1]/N;atol=1e-12)
    @assert isapprox(d.empirical_regional_retention[end-1],mean(K)/N;atol=1e-12)
    worst=argmax(abs.(err));peak=argmax(main.empirical_gap)
    push!(summary,(dataset=ds,N=N,regional_links=L,comparison_levels=nrow(main),maximum_removal_fraction=last(main.removal_fraction),
        mean_absolute_error_pp=mean(abs.(err)),maximum_absolute_error_pp=maximum(abs.(err)),signed_mean_error_pp=mean(err),
        worst_removed_sites=main.removed_sites[worst],worst_removal_fraction=main.removal_fraction[worst],signed_error_at_worst_pp=err[worst],
        error_at_last_main_level_pp=last(err),empirical_peak_gap_pp=100main.empirical_gap[peak],empirical_peak_removal_fraction=main.removal_fraction[peak],
        s=par.s,r=par.r))
end
d=DataFrame(curves);s=DataFrame(summary)
CSV.write(joinpath(DEST,"curves.csv"),d);CSV.write(joinpath(DEST,"summary.csv"),s)
set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=15,yticklabelsize=15)))
main=d[d.in_main_range,:];ymax=5ceil(maximum(vcat(main.empirical_gap,main.fitted_gap))*100/5)
errmax=max(1.,ceil(maximum(abs.(main.error_pp))))
for kind in ["gap","error"]
    f=Figure(size=(1850,850))
    for (j,ds) in enumerate(DATA)
        v=main[main.dataset.==ds,:]
        ax=Axis(f[cld(j,5),mod1(j,5)],title=DISPLAY[ds],titlesize=19,xlabel="Sites removed (%)",
            ylabel=mod1(j,5)==1 ? (kind=="gap" ? "Regional − local retention (pp)" : "Fitted − empirical gap (pp)") : "",xticks=0:20:80)
        if kind=="gap"
            lines!(ax,100v.removal_fraction,100v.empirical_gap;color=PALETTE[j],linewidth=3)
            lines!(ax,100v.removal_fraction,100v.fitted_gap;color=:gray25,linestyle=:dash,linewidth=2)
            ylims!(ax,0,ymax)
        else
            hlines!(ax,[0],color=:gray65,linestyle=:dash)
            lines!(ax,100v.removal_fraction,v.error_pp;color=PALETTE[j],linewidth=2.5)
            scatter!(ax,100v.removal_fraction,v.error_pp;color=PALETTE[j],markersize=4)
            ylims!(ax,-errmax,errmax)
        end
        xlims!(ax,0,80)
    end
    if kind=="gap"
        Legend(f[3,1:5],[LineElement(color=:black,linewidth=3),LineElement(color=:gray25,linestyle=:dash)],
            ["Empirical support","Fitted singleton + geometric support"];orientation=:horizontal,framevisible=false)
    end
    save(joinpath(DEST,kind=="gap" ? "01_gap_reconstruction.png" : "02_gap_errors.png"),f;px_per_unit=2)
end
println(s)
end

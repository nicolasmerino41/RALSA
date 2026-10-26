const SCENARIOS=["Observed_M","Known_Proportional","All_Proportional","Known_Intercept","All_Intercept"]
const DISPLAY=["Observed_M","Known_Proportional","All_Proportional"]
const COLOR=Dict("Observed_M"=>"#CC6677","Known_Proportional"=>"#228833","All_Proportional"=>"#4477AA")
const LABEL=Dict("Observed_M"=>"Observed M","Known_Proportional"=>"M from known-link CO","All_Proportional"=>"M from all-pair CO")
function compare()
    allroot=joinpath(SOURCE,"cooccurrence_estimated_M");oldroot=joinpath(SOURCE,"conditioned_M")
    used=CSV.read(joinpath(DEST,"source_versions.csv"),DataFrame)
    for file in [joinpath(allroot,"source_versions.csv"),joinpath(oldroot,"input_versions.csv")]
        for r in eachrow(CSV.read(file,DataFrame));@assert bytes2hex(sha256(read(joinpath(ROOT,r.path))))==r.sha256 "Changed comparison input: $(r.path)";end
    end
    function load(path)
        push!(used,(relpath(path,ROOT),bytes2hex(sha256(read(path)))))
        CSV.read(path,DataFrame)
    end
    newf=load(joinpath(DEST,"support_frequencies.csv"));allf=load(joinpath(allroot,"support_frequencies.csv"));oldf=load(joinpath(oldroot,"support_frequencies.csv"))
    newr=load(joinpath(DEST,"regression_folds.csv"));allr=load(joinpath(allroot,"regression_folds.csv"));totals=load(joinpath(DEST,"dataset_totals.csv"));oldcurves=load(joinpath(oldroot,"curves.csv"))
    oldtotals=load(joinpath(allroot,"dataset_totals.csv"));@assert totals==oldtotals
    curves=NamedTuple[];distributions=NamedTuple[];metrics=NamedTuple[];decomp=NamedTuple[]
    for ds in DATA,model in ["Pair","Common"]
        t=only(eachrow(totals[totals.dataset.==ds,:]));ref=oldf[(oldf.dataset.==ds).&(oldf.model.==model),:]
        @assert ref.K==collect(1:t.N)
        empirical=ref.empirical_frequency;oracle=ref.predicted_frequency;H=gapmatrix(t.N);eg=100H*empirical;og=100H*oracle;idx=1:fld(4t.N,5)+1
        oldc=oldcurves[(oldcurves.dataset.==ds).&(oldcurves.model.==model),:]
        @assert maximum(abs.(oldc.predicted_gap_pp.-og))<1e-8 && maximum(abs.(oldc.empirical_gap_pp.-eg))<1e-8
        for scenario in SCENARIOS
            isold=scenario=="Observed_M";known=startswith(scenario,"Known");kind=endswith(scenario,"Intercept") ? "Intercept" : "Proportional"
            if isold
                p=oracle;raw=Float64(t.M);usedM=t.M;flag="Observed reference"
            else
                source=known ? newf : allf;reg=known ? newr : allr
                v=source[(source.dataset.==ds).&(source.model.==model).&(source.scenario.==kind),:];@assert v.K==collect(1:t.N)
                p=v.predicted_frequency;r=only(eachrow(reg[(reg.dataset.==ds).&(reg.scenario.==kind),:]));raw=r.predicted_M_raw;usedM=r.predicted_M_used;flag=r.feasibility
                @assert Set(split(r.training_datasets,";"))==Set(filter(!=(ds),DATA))
            end
            @assert isapprox(sum(p),1;atol=1e-8) && isapprox(dot(1:t.N,p),usedM/t.L;atol=1e-8)
            gap=100H*p;error=gap[idx].-eg[idx];cdfdiff=cumsum(p.-empirical)
            tv=.5sum(abs.(p.-empirical));alloc=.5sum(abs.(oracle.-empirical));shift=.5sum(abs.(p.-oracle))
            @assert maximum(abs.((p.-empirical).-((oracle.-empirical).+(p.-oracle))))<1e-12
            @assert abs(sum(cdfdiff)-(t.M-usedM)/t.L)<1e-7
            push!(metrics,(dataset=ds,model=model,scenario=scenario,N=t.N,L=t.L,observed_M=t.M,predicted_M_raw=raw,predicted_M_used=usedM,feasibility=flag,relative_M_error_percent=100(raw/t.M-1),absolute_relative_M_error_percent=100abs(raw/t.M-1),mean_support_error=(usedM-t.M)/t.L,TV_percent=100tv,KS_pp=100maximum(abs.(cdfdiff)),Wasserstein_sites=sum(abs.(cdfdiff)),singleton_error_pp=100(p[1]-empirical[1]),observed_p90=findfirst(>=(.9),cumsum(empirical)),predicted_p90=findfirst(>=(.9),cumsum(p)),allocation_TV_percent=100alloc,amount_shift_TV_percent=100shift,gap_mae_pp=mean(abs.(error)),gap_signed_error_pp=mean(error),gap_max_error_pp=maximum(abs.(error)),mean_empirical_gap_pp=mean(eg[idx]),relative_gap_error_percent=100mean(abs.(error))/mean(eg[idx]),allocation_gap_mae_pp=mean(abs.(og[idx].-eg[idx])),amount_shift_gap_mae_pp=mean(abs.(gap[idx].-og[idx]))))
            empcc=reverse(cumsum(reverse(empirical)));pcc=reverse(cumsum(reverse(p)))
            for k in 1:t.N
                push!(distributions,(dataset=ds,model=model,scenario=scenario,K=k,empirical_frequency=empirical[k],predicted_frequency=p[k],empirical_ccdf=empcc[k],predicted_ccdf=pcc[k]))
                push!(decomp,(dataset=ds,model=model,scenario=scenario,quantity="frequency",position=k,total_error=p[k]-empirical[k],allocation_error=oracle[k]-empirical[k],amount_effect=p[k]-oracle[k]))
            end
            for m in 0:t.N
                push!(curves,(dataset=ds,model=model,scenario=scenario,m=m,removal_fraction=m/t.N,empirical_gap_pp=eg[m+1],predicted_gap_pp=gap[m+1]))
                push!(decomp,(dataset=ds,model=model,scenario=scenario,quantity="gap_pp",position=m,total_error=gap[m+1]-eg[m+1],allocation_error=og[m+1]-eg[m+1],amount_effect=gap[m+1]-og[m+1]))
            end
        end
    end
    tab=DataFrame(metrics);agg=NamedTuple[]
    for model in ["Pair","Common"],scenario in SCENARIOS
        v=tab[(tab.model.==model).&(tab.scenario.==scenario),:];baseline=tab[(tab.model.==model).&(tab.scenario.=="All_Proportional"),:]
        push!(agg,(model=model,scenario=scenario,mean_absolute_relative_M_error_percent=mean(v.absolute_relative_M_error_percent),mean_TV_percent=mean(v.TV_percent),median_TV_percent=median(v.TV_percent),mean_KS_pp=mean(v.KS_pp),mean_Wasserstein_sites=mean(v.Wasserstein_sites),mean_absolute_singleton_error_pp=mean(abs.(v.singleton_error_pp)),mean_gap_mae_pp=mean(v.gap_mae_pp),median_gap_mae_pp=median(v.gap_mae_pp),mean_gap_signed_error_pp=mean(v.gap_signed_error_pp),mean_relative_gap_error_percent=mean(v.relative_gap_error_percent),TV_improved_vs_all=count(v.TV_percent.<baseline.TV_percent),gap_improved_vs_all=count(v.gap_mae_pp.<baseline.gap_mae_pp)))
    end
    CSV.write(joinpath(DEST,"distribution_metrics.csv"),tab);CSV.write(joinpath(DEST,"comparison_summary.csv"),DataFrame(agg));CSV.write(joinpath(DEST,"comparison_distributions.csv"),DataFrame(distributions));CSV.write(joinpath(DEST,"comparison_curves.csv"),DataFrame(curves));CSV.write(joinpath(DEST,"signed_error_decomposition.csv"),DataFrame(decomp));CSV.write(joinpath(DEST,"comparison_source_versions.csv"),used)
end
function plots()
    t=CSV.read(joinpath(DEST,"distribution_metrics.csv"),DataFrame);d=CSV.read(joinpath(DEST,"comparison_distributions.csv"),DataFrame);c=CSV.read(joinpath(DEST,"comparison_curves.csv"),DataFrame)
    set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
    for (kind,file) in [("frequency","01_support_frequencies"),("ccdf","02_support_ccdf"),("gap","03_gap_curves")]
        f=Figure(size=(1800,880))
        for (i,ds) in enumerate(DATA)
            z=only(eachrow(t[(t.dataset.==ds).&(t.model.=="Pair").&(t.scenario.=="Observed_M"),:]))
            ax=Axis(f[cld(i,5),mod1(i,5)];title=replace(ds,"_"=>" "),xlabel=kind=="gap" ? "Sites removed (%)" : "Supporting sites, K",ylabel=mod1(i,5)==1 ? (kind=="gap" ? "Gap (pp)" : kind=="ccdf" ? "Proportion with support ≥ K" : "Proportion of links") : "",yscale=kind=="ccdf" ? log10 : identity,xscale=kind=="ccdf" ? log10 : identity)
            maxk=1
            for scenario in DISPLAY
                if kind=="gap"
                    v=c[(c.dataset.==ds).&(c.model.=="Pair").&(c.scenario.==scenario).&(c.removal_fraction.<=.8),:]
                    lines!(ax,100v.removal_fraction,v.predicted_gap_pp;color=COLOR[scenario],linewidth=2.5)
                    scenario=="Observed_M" && lines!(ax,100v.removal_fraction,v.empirical_gap_pp;color=:black,linestyle=:dash,linewidth=2.5)
                else
                    v=d[(d.dataset.==ds).&(d.model.=="Pair").&(d.scenario.==scenario),:]
                    field=kind=="ccdf" ? :predicted_ccdf : :predicted_frequency;empfield=kind=="ccdf" ? :empirical_ccdf : :empirical_frequency
                    mask=v[!,field].>0;lines!(ax,v.K[mask],v[mask,field];color=COLOR[scenario],linewidth=2)
                    maxk=max(maxk,maximum(v.K[v.empirical_frequency.>0]),something(findfirst(>=(.995),cumsum(v.predicted_frequency)),z.N))
                    if scenario == "Observed_M"
                        keep = v[!, empfield] .> 0

                        if kind == "ccdf"
                            # Straight lines instead of stairs
                            lines!(
                                ax, v.K[keep], v[keep, empfield];
                                color=:black,
                                linewidth=1.8
                            )
                        else
                            # Connect empirical frequency points with straight lines
                            lines!(
                                ax, v.K, v[!, empfield];
                                color=:black,
                                linewidth=1.8
                            )

                            # Retain empirical points
                            scatter!(
                                ax, v.K[keep], v[keep, empfield];
                                color=:black,
                                markersize=6
                            )
                        end
                    end
                end
            end
            if kind=="gap";xlims!(ax,0,80);ylims!(ax,0,50)
            elseif kind=="frequency"
                xlims!(ax,.5,min(z.N,maxk+1)+.5);ylims!(ax,0,.85)
                ticks=maxk<=14 ? filter(x->x<=maxk+1,[1,2,4,6,8,10,12,14]) : filter(x->x<=maxk+1,[1,5,10,15,20,25,30,40,50,100,200])
                ax.xticks=(ticks,string.(ticks))
            else;xlims!(ax,.95,min(z.N,maxk+1)+.5);ylims!(ax,.5/z.L,1.1);ticks=filter(x->x<=maxk+1,[1,2,5,10,20,50,100,200]);ax.xticks=(ticks,string.(ticks));end
        end
        empirical_element=kind=="frequency" ? MarkerElement(color=:black,marker=:circle) : LineElement(color=:black,linestyle=kind=="gap" ? :dash : :solid)
        Legend(f[3,1:5],vcat([LineElement(color=COLOR[s]) for s in DISPLAY],[empirical_element]),vcat([LABEL[s] for s in DISPLAY],["Empirical"]);orientation=:horizontal,framevisible=false)
        save(joinpath(DEST,file*".png"),f;px_per_unit=2)
    end
    f=Figure(size=(1760,800))
    for (j,field,label) in [(1,:relative_M_error_percent,"Error in estimated M (%)"),(2,:TV_percent,"Distribution difference: TV (%)"),(3,:gap_mae_pp,"Gap MAE (pp)")]
        ax=Axis(f[1,j];xlabel=label,yticks=(1:10,replace.(DATA,"_"=>" ")));j==1 && vlines!(ax,[0];color=:gray,linestyle=:dash)
        for (i,ds) in enumerate(DATA),(si,scenario) in enumerate(DISPLAY)
            j==1 && scenario=="Observed_M" && continue
            z=only(eachrow(t[(t.dataset.==ds).&(t.model.=="Pair").&(t.scenario.==scenario),:]));scatter!(ax,[z[field]],[i+(si-2)*.17];color=COLOR[scenario],markersize=10)
        end
        ylims!(ax,.5,10.5);ax.yreversed=true
    end
    Legend(f[2,1:3],[MarkerElement(color=COLOR[s],marker=:circle) for s in DISPLAY],[LABEL[s] for s in DISPLAY];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"04_amount_and_allocation_errors.png"),f;px_per_unit=2)
end
function report()
    a=CSV.read(joinpath(DEST,"comparison_summary.csv"),DataFrame);t=CSV.read(joinpath(DEST,"distribution_metrics.csv"),DataFrame);reg=CSV.read(joinpath(DEST,"full_data_regressions.csv"),DataFrame);r=CSV.read(joinpath(DEST,"regression_folds.csv"),DataFrame)
    d=CSV.read(joinpath(DEST,"comparison_distributions.csv"),DataFrame);c=CSV.read(joinpath(DEST,"comparison_curves.csv"),DataFrame);de=CSV.read(joinpath(DEST,"signed_error_decomposition.csv"),DataFrame);checks=NamedTuple[]
    @assert maximum(abs.(de.total_error.-de.allocation_error.-de.amount_effect))<1e-10
    for z in eachrow(t)
        f=d[(d.dataset.==z.dataset).&(d.model.==z.model).&(d.scenario.==z.scenario),:];v=c[(c.dataset.==z.dataset).&(c.model.==z.model).&(c.scenario.==z.scenario),:]
        @assert v.m==collect(0:z.N) && f.K==collect(1:z.N)
        @assert maximum(abs.(100gapmatrix(z.N)*f.predicted_frequency.-v.predicted_gap_pp))<1e-8
        @assert isapprox(50sum(abs.(f.predicted_frequency.-f.empirical_frequency)),z.TV_percent;atol=1e-9)
        @assert isapprox(mean(abs.((v.predicted_gap_pp.-v.empirical_gap_pp)[v.removal_fraction.<=.8])),z.gap_mae_pp;atol=1e-9)
        @assert abs(dot(f.K,f.predicted_frequency)-z.predicted_M_used/z.L)<1e-8
        push!(checks,(dataset=z.dataset,model=z.model,scenario=z.scenario,passed=true))
    end
    CSV.write(joinpath(DEST,"comparison_checks.csv"),DataFrame(checks))
    fmt(x)=@sprintf("%.3f",x)
    open(joinpath(DEST,"FINDINGS.md"),"w") do io
        println(io,"# Known-link opportunities: amount versus allocation\n\nAll ten datasets are retained. Both total-support regressions exclude the target, as do the unchanged persistent-pair parameters. Known regional interaction identities and occurrence-derived opportunities remain supplied.\n")
        for z in eachrow(reg);println(io,"All-ten descriptive $(z.model) fit: M=$(fmt(z.intercept)) + $(fmt(z.slope)) T_known; slope $(fmt(z.slope_percent))%, Pearson r=$(fmt(z.pearson_r)). These all-ten coefficients do not generate held-out predictions.\n");end
        println(io,"| Pair model input | Mean absolute M error (%) | Mean distribution TV (%) | Mean KS (pp) | Mean gap MAE (pp) |\n|---|---:|---:|---:|---:|")
        for z in eachrow(a[a.model.=="Pair",:]);println(io,"| $(z.scenario) | $(fmt(z.mean_absolute_relative_M_error_percent)) | $(fmt(z.mean_TV_percent)) | $(fmt(z.mean_KS_pp)) | $(fmt(z.mean_gap_mae_pp)) |");end
        known=only(eachrow(a[(a.model.=="Pair").&(a.scenario.=="Known_Proportional"),:]));all=only(eachrow(a[(a.model.=="Pair").&(a.scenario.=="All_Proportional"),:]));obs=only(eachrow(a[(a.model.=="Pair").&(a.scenario.=="Observed_M"),:]))
        println(io,"\nKnown-link CO changes mean Pair gap MAE from $(fmt(all.mean_gap_mae_pp)) to $(fmt(known.mean_gap_mae_pp)) pp, improving $(known.gap_improved_vs_all)/10 datasets. Distribution TV changes from $(fmt(all.mean_TV_percent))% to $(fmt(known.mean_TV_percent))%, improving $(known.TV_improved_vs_all)/10. With observed M, the allocation-only reference still has TV=$(fmt(obs.mean_TV_percent))%, despite gap MAE=$(fmt(obs.mean_gap_mae_pp)) pp. A good gap curve therefore does not imply an equally close support-frequency distribution.\n")
        println(io,"The known-link predictor substantially reduces total-support error (mean absolute relative error 28.8%, versus 58.2% for all-pair CO). Salix Galpar and Garraf PP now have almost correct amounts and approach the observed-M reference. Major remaining amount errors occur in Gottin HP (+75.6%), Gottin PP (+98.4%) and Garraf PP2 (-51.1%). Those systems also show large distribution and gap errors. This points to poor amount estimation as an important remaining limitation there. Known-link prediction is not uniformly superior: Nahuel, Olot and Montseny have larger gap errors than under all-pair prediction.\n\nEven with observed M, distribution TV ranges from 4.3% to 15.6%; the support-frequency distribution is not recovered exactly. Some small-frequency tail differences remain visible in CCDFs. TV includes empirical sampling fluctuations, so this is descriptive mismatch rather than a test rejecting the allocation model. With known-link-estimated M, Pair still improves mean gap MAE over Common (3.252 versus 5.056 pp), and mean distribution TV (14.601% versus 19.080%). Allocation heterogeneity therefore retains descriptive value alongside the unresolved amount prediction.\n")
        println(io,"| Dataset | M error: known CO (%) | TV: observed M (%) | TV: known CO (%) | TV: all CO (%) | Gap MAE: known CO (pp) | Gap MAE: all CO (pp) |\n|---|---:|---:|---:|---:|---:|---:|")
        for ds in DATA
            v=t[(t.dataset.==ds).&(t.model.=="Pair"),:];k=only(eachrow(v[v.scenario.=="Known_Proportional",:]));o=only(eachrow(v[v.scenario.=="Observed_M",:]));b=only(eachrow(v[v.scenario.=="All_Proportional",:]))
            println(io,"| $(replace(ds,"_"=>" ")) | $(fmt(k.relative_M_error_percent)) | $(fmt(o.TV_percent)) | $(fmt(k.TV_percent)) | $(fmt(b.TV_percent)) | $(fmt(k.gap_mae_pp)) | $(fmt(b.gap_mae_pp)) |")
        end
        println(io,"\n## Separation\n\nObserved-M predictions isolate the remaining allocation mismatch within this model: their mean support equals observed M/L exactly. Estimated-M predictions add a change in amount while holding pair parameters, opportunities, known links and N fixed. The exact signed decomposition is predicted minus empirical = (observed-M model minus empirical) + (estimated-M model minus observed-M model). It is saved for every frequency bin and removal level. Absolute distances are not additive: the components can reinforce or cancel. An estimated-M curve can occasionally improve through compensation, not necessarily better allocation. This is a model-conditional diagnostic, not a causal or variance decomposition.\n")
        println(io,"TV is half the sum of absolute frequency differences: the fraction of probability mass that must move among K bins. KS is maximum CDF separation; Wasserstein is summed absolute CDF separation in supporting-site units. Singleton error and empirical/model 90th percentiles are also exported. Metrics use expected model distributions, without tests of sampling significance.\n")
        flagged=r[r.feasibility.!="Feasible",:];println(io,"Known-link regressions have $(nrow(flagged)) infeasible rounded predictions. All raw estimates are retained; predictions are projected to [L,sum(C)] only where necessary, as in the earlier all-pair analysis.")
        for z in eachrow(flagged);println(io,"- $(z.dataset) $(z.scenario): raw $(fmt(z.predicted_M_raw)), used $(z.predicted_M_used); $(z.feasibility).");end
        println(io,"\nOLS with an intercept is retained as a sensitivity check, not selected based on held-out performance. Gap errors average integer removal levels through 80%, including zero. Across-dataset summaries give equal weight to each dataset. Shared sites, collection dependence, interaction-derived occupancy, unresolved corrections, known-link conditioning and unpropagated parameter/regression uncertainty remain; see README. No datasets were removed, and no allocation parameters were refitted.")
    end
    CSV.write(joinpath(DEST,"code_versions.csv"),DataFrame(file=["predict_M.jl","compare.jl"],sha256=[bytes2hex(sha256(read(joinpath(@__DIR__,f)))) for f in ["predict_M.jl","compare.jl"]]))
end

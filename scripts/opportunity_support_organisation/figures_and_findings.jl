const PALETTE=Dict("Opportunity"=>"#4477AA","Support"=>"#111111","Common"=>"#999999","Pair"=>"#CC6677")
const LABELS=Dict("Opportunity"=>"Co-occurrence C","Support"=>"Observed support K","Common"=>"Common probability, fixed M","Pair"=>"Persistent pair, fixed M")
function legend!(fig,row,cols,quantities)
    Legend(fig[row,cols],[LineElement(color=PALETTE[q],linestyle=q=="Support" ? :dash : :solid) for q in quantities],[LABELS[q] for q in quantities];orientation=:horizontal,framevisible=false,labelsize=17)
end
function plots()
    d=CSV.read(joinpath(DEST,"distributions.csv"),DataFrame);a=CSV.read(joinpath(DEST,"allocation_summary.csv"),DataFrame);c=CSV.read(joinpath(DEST,"conditional_C.csv"),DataFrame);l=CSV.read(joinpath(DEST,"concentration_curves.csv"),DataFrame);s=CSV.read(joinpath(DEST,"site_subsampling_summary.csv"),DataFrame);e=CSV.read(joinpath(DEST,"empirical_summary.csv"),DataFrame)
    set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,xgridvisible=false,ygridvisible=false,xticklabelsize=16,yticklabelsize=16)))
    for kind in ["ccdf","frequency"]
        f=Figure(size=(1800,880))
        for (i,ds) in enumerate(DATA)
            ax=Axis(f[cld(i,5),mod1(i,5)];title=replace(ds,"_"=>" "),xlabel="Number of sites",ylabel=mod1(i,5)==1 ? (kind=="ccdf" ? "Proportion with ≥ x sites" : "Proportion of links") : "",xscale=kind=="ccdf" ? log10 : identity,yscale=kind=="ccdf" ? log10 : identity)
            L=only(e.L[(e.dataset.==ds).&(e.quantity.=="Support")]);maxk=maximum(d.K[(d.dataset.==ds).&(d.quantity.=="Opportunity").&(d.frequency.>0)])
            for q in ["Opportunity","Common","Pair","Support"]
                v=d[(d.dataset.==ds).&(d.quantity.==q),:];field=Symbol(kind);keep=kind=="frequency" ? trues(nrow(v)) : v[!,field].>0
                if kind=="ccdf" && q in ["Opportunity","Support"]
                    stairs!(ax,v.K[keep],v[keep,field];color=PALETTE[q],linewidth=2,step=:post,linestyle=q=="Support" ? :dash : :solid)
                else
                    lines!(ax,v.K[keep],v[keep,field];color=PALETTE[q],linewidth=2,linestyle=q=="Support" ? :dash : :solid)
                end
            end
            if kind=="ccdf"
                xlims!(ax,.95,maxk+.5);ylims!(ax,.5/L,1.1);ticks=filter(x->x<=maxk,[1,2,5,10,20,50,100,200]);ax.xticks=(ticks,string.(ticks))
            else;xlims!(ax,.5,maxk+.5);ylims!(ax,0,.85);end
        end
        legend!(f,3,1:5,["Opportunity","Support","Common","Pair"])
        save(joinpath(DEST,kind=="ccdf" ? "01_opportunities_and_support.png" : "08_support_frequencies.png"),f;px_per_unit=2)
    end
    f=Figure(size=(1800,820))
    fields=[("singleton","Singleton links (%)"),("gini","Support Gini (%)"),("top10","Support in top 10% of links (%)"),("within_C_variance","Within-C realisation spread (pp)")]
    for (j,(field,label)) in enumerate(fields)
        ax=Axis(f[1,j];xlabel=label,yticks=(1:10,replace.(DATA,"_"=>" ")))
        transform(x)=field=="within_C_variance" ? 100sqrt(max(0,x)) : 100x
        for (i,ds) in enumerate(DATA),model in ["Common","Pair"]
            v=a[(a.dataset.==ds).&(a.model.==model).&(a.metric.==field),:];isempty(v) && continue;r=only(eachrow(v));y=i+(model=="Common" ? -.15 : .15)
            lines!(ax,transform.([r.lower,r.upper]),[y,y];color=PALETTE[model],linewidth=3);scatter!(ax,[transform(r.predicted_mean)],[y];color=PALETTE[model],markersize=8)
            model=="Pair" && scatter!(ax,[transform(r.observed)],[i];color=:black,marker=:diamond,markersize=11)
        end
        ylims!(ax,.5,10.5);ax.yreversed=true
    end
    Legend(f[2,1:4],[LineElement(color=PALETTE["Common"]),LineElement(color=PALETTE["Pair"]),MarkerElement(color=:black,marker=:diamond)],["Common: allocation mean and 95% interval","Pair: allocation mean and 95% interval","Empirical"];orientation=:horizontal,framevisible=false,labelsize=16)
    save(joinpath(DEST,"02_allocation_features.png"),f;px_per_unit=2)
    for kind in ["mean","sd"]
        f=Figure(size=(1800,880))
        for (i,ds) in enumerate(DATA)
            ax=Axis(f[cld(i,5),mod1(i,5)];title=replace(ds,"_"=>" "),xlabel="Co-occurring sites, C",ylabel=mod1(i,5)==1 ? (kind=="mean" ? "Mean support K at this C" : "SD of support K at this C") : "")
            for model in ["Common","Pair"]
                v=c[(c.dataset.==ds).&(c.model.==model).&c.display,:];sort!(v,:C);isempty(v) && continue
                band!(ax,v.C,v[!,Symbol("lower_"*kind*"_K")],v[!,Symbol("upper_"*kind*"_K")];color=(PALETTE[model],.15))
                lines!(ax,v.C,v[!,Symbol("predicted_"*kind*"_K")];color=PALETTE[model],linewidth=2)
                model=="Pair" && scatter!(ax,v.C,v[!,Symbol("observed_"*kind*"_K")];color=:black,markersize=7)
            end
            ylims!(ax,0,nothing)
        end
        Legend(f[3,1:5],[LineElement(color=PALETTE["Common"]),LineElement(color=PALETTE["Pair"]),MarkerElement(color=:black,marker=:circle)],["Common: allocation mean and 95% interval","Pair: allocation mean and 95% interval","Empirical (≥3 links per C)"];orientation=:horizontal,framevisible=false,labelsize=17)
        save(joinpath(DEST,kind=="mean" ? "03_support_at_equal_opportunity.png" : "04_variation_at_equal_opportunity.png"),f;px_per_unit=2)
    end
    f=Figure(size=(1800,880))
    for (i,ds) in enumerate(DATA)
        ax=Axis(f[cld(i,5),mod1(i,5)];title=replace(ds,"_"=>" "),xlabel="Cumulative share of links",ylabel=mod1(i,5)==1 ? "Cumulative share of occurrences" : "")
        lines!(ax,[0,1],[0,1];color=:gray,linestyle=:dot)
        for q in ["Opportunity","Common","Pair","Support"]
            v=l[(l.dataset.==ds).&(l.quantity.==q),:]
            q in ["Common","Pair"] && band!(ax,v.link_fraction,v.lower_share,v.upper_share;color=(PALETTE[q],.12))
            lines!(ax,v.link_fraction,v.mean_share;color=PALETTE[q],linewidth=2,linestyle=q=="Support" ? :dash : :solid)
        end
        xlims!(ax,0,1);ylims!(ax,0,1)
    end
    legend!(f,3,1:5,["Opportunity","Support","Common","Pair"]);save(joinpath(DEST,"05_concentration.png"),f;px_per_unit=2)
    f=Figure(size=(1800,820))
    for (j,field,label) in [(1,"singleton","Singleton links (%)"),(2,"gini","Support Gini (%)"),(3,"top10","Support in top 10% of links (%)"),(4,"lost_link_fraction","Original links unobserved at 8 sites (%)")]
        ax=Axis(f[1,j];xlabel=label,yticks=(1:10,replace.(DATA,"_"=>" ")))
        for (i,ds) in enumerate(DATA)
            v=only(eachrow(s[(s.dataset.==ds).&(s.sampled_sites.==8).&(s.metric.==field),:]));lines!(ax,100 .* [v.lower,v.upper],[i,i];color="#4477AA",linewidth=3);scatter!(ax,[100v.mean],[i];color="#4477AA",markersize=9)
            if j<4;full=only(eachrow(e[(e.dataset.==ds).&(e.quantity.=="Support"),:]));scatter!(ax,[100full[Symbol(field)]],[i];color=:black,marker=:diamond,markersize=10);end
        end
        ylims!(ax,.5,10.5);ax.yreversed=true
    end
    Legend(f[2,1:4],[LineElement(color="#4477AA"),MarkerElement(color=:black,marker=:diamond)],["Eight sites: mean and central 95% subsampling interval","Full observed site set"];orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"06_sampling_sensitivity.png"),f;px_per_unit=2)
    v=CSV.read(joinpath(DEST,"eight_site_distributions.csv"),DataFrame);f=Figure(size=(1800,880))
    for (i,ds) in enumerate(DATA)
        ax=Axis(f[cld(i,5),mod1(i,5)];title=replace(ds,"_"=>" "),xlabel="Number of sites (out of 8)",ylabel=mod1(i,5)==1 ? "Mean proportion with ≥ x sites" : "",yscale=log10,xticks=[1,2,4,6,8])
        g=v[v.dataset.==ds,:]
        for (field,q) in [(:support_frequency,"Support"),(:opportunity_frequency,"Opportunity")]
            cc=reverse(cumsum(reverse(g[!,field])));keep=cc.>0;lines!(ax,g.K[keep],cc[keep];color=PALETTE[q],linewidth=2.5,linestyle=q=="Support" ? :dash : :solid)
        end
        xlims!(ax,1,8);ylims!(ax,1e-4,1.1)
    end
    legend!(f,3,1:5,["Opportunity","Support"]);save(joinpath(DEST,"07_eight_site_comparison.png"),f;px_per_unit=2)
end
function report()
    a=CSV.read(joinpath(DEST,"allocation_summary.csv"),DataFrame);s=CSV.read(joinpath(DEST,"site_subsampling_summary.csv"),DataFrame);e=CSV.read(joinpath(DEST,"empirical_summary.csv"),DataFrame)
    dist=CSV.read(joinpath(DEST,"distributions.csv"),DataFrame);checks=NamedTuple[]
    for v in groupby(dist,[:dataset,:quantity])
        @assert isapprox(sum(v.frequency),1;atol=1e-8)
        @assert maximum(abs.(reverse(cumsum(reverse(v.frequency))).-v.ccdf))<1e-8
        ds=v.dataset[1];q=v.quantity[1];ref=q in ["Common","Pair"] ? "Support" : q
        r=only(eachrow(e[(e.dataset.==ds).&(e.quantity.==ref),:]));@assert abs(dot(v.K,v.frequency)-r.total/r.L)<1e-8
        if q in ["Common","Pair"]
            z=only(eachrow(a[(a.dataset.==ds).&(a.model.==q).&(a.metric.=="singleton"),:]))
            @assert abs(z.predicted_mean-v.frequency[1])<6z.mean_mcse+1e-10
        end
        push!(checks,(dataset=ds,quantity=q,passed=true))
    end
    rr=CSV.read(joinpath(DEST,"site_subsampling_draws.csv"),DataFrame)
    for ds in DATA
        x=only(eachrow(e[(e.dataset.==ds).&(e.quantity.=="Support"),:]));v=rr[(rr.dataset.==ds).&(rr.sampled_sites.==x.N),:]
        for field in [:singleton,:gini,:top10,:mean_support];@assert maximum(abs.(v[!,field].-x[field]))<1e-10;end
        @assert all(v.retained_links.==x.L)
    end
    CSV.write(joinpath(DEST,"result_checks.csv"),DataFrame(checks))
    fmt(x)=@sprintf("%.3f",x)
    open(joinpath(DEST,"FINDINGS.md"),"w") do io
        println(io,"# Opportunity, realisation and local support\n\nThis analysis holds known links, site count, opportunities and observed total support fixed to compare allocation. Pair parameters are unchanged nine-dataset estimates. Common is the uniform eligible positive-row allocation conditional on M. Results concern whether opportunities alone plus a shared probability suffice, and which features persistent pair heterogeneity improves. They do not identify a biological mechanism.\n")
        println(io,"The strongest recurring pattern is a departure from common-probability allocation: all ten empirical datasets exceed its central 95% intervals for singleton abundance, Gini and support in the highest-supported 10% of links. Opportunity heterogeneity and known total support therefore do not reproduce these features under this null. Persistent-pair heterogeneity moves concentration toward the observations in all ten datasets and singleton abundance in nine, but does not fully reproduce the distributions.\n\nAt equal opportunities, the picture is more specific: nine datasets exceed the Common interval for within-C realisation variance, yet the transferred Pair model overshoots the empirical variance (empirical below its interval) in seven. Pair is closer than Common on this scalar in only three datasets. Thus heterogeneous realisation is useful for the marginal distribution, but the existing two-state specification does not correctly distribute all that heterogeneity among opportunity classes. This is an important limitation that gap accuracy alone concealed.\n\nSite count has a substantial effect on apparent shape: reducing Salix Galpar from 341 to eight sites raises the average singleton fraction from 45.5% to 93.8%, while about 92.7% of its original links disappear from the subsample. Other systems remain more broadly supported at eight sites. These results support recurring concentration relative to the Common null, not a sampling-invariant universal support distribution. The null comparison itself has not yet been repeated on subsamples.\n")
        println(io,"## Which features need more than common realisation?\n\n| Feature | Empirical above Common interval | Below Common interval | Inside Common interval | Inside Pair interval | Pair closer than Common |\n|---|---:|---:|---:|---:|---:|")
        for metric in ["singleton","gini","top10","within_C_variance"]
            x=a[(a.model.=="Common").&(a.metric.==metric),:];y=a[(a.model.=="Pair").&(a.metric.==metric),:];@assert x.dataset==y.dataset
            println(io,"| $metric | $(count(x.observed.>x.upper)) | $(count(x.observed.<x.lower)) | $(count(.!x.outside_interval)) | $(count(.!y.outside_interval)) | $(count(abs.(y.observed.-y.predicted_mean).<abs.(x.observed.-x.predicted_mean))) / $(nrow(x)) |")
        end
        println(io,"\nIntervals are central 95% conditional allocation intervals from 499 independent exact draws, not confidence intervals, removal-sequence variation or multiple-testing-corrected significance tests. Empirical positions are descriptive ranks. Mean Monte Carlo SE is exported.\n\nOpportunity/support CCDFs compare the same known links. Common preserves each link's C and total M: its discrepancies cannot be attributed merely to differing opportunities or total support. Pair is a specific transferred two-state model; improvements support persistent heterogeneity as a useful statistical description, not proof of a two-class ecological mechanism. Observation differences could also contribute.\n")
        println(io,"## Variation at equal opportunity\n\nEach conditional point groups links with exactly the same C; figures show only groups with at least three links. Both empirical and model summaries use identical groups. The within-C scalar is the link-weighted mean squared deviation of K/C from its within-class mean, across eligible classes. The figure displays its square root in percentage points. This prevents variation in C alone from creating the measured within-class spread. Sparse groups and omitted classes are retained in conditional_C.csv, with link counts.\n")
        println(io,"## Eight-site sampling check\n\n| Dataset | Full singleton % | Eight-site mean singleton % | Full Gini % | Eight-site mean Gini % | Mean original links lost % |\n|---|---:|---:|---:|---:|---:|")
        for ds in DATA
            r=only(eachrow(e[(e.dataset.==ds).&(e.quantity.=="Support"),:]));v=s[(s.dataset.==ds).&(s.sampled_sites.==8),:]
            getmean(field)=only(v.mean[v.metric.==field])
            println(io,"| $(replace(ds,"_"=>" ")) | $(fmt(100r.singleton)) | $(fmt(100getmean("singleton"))) | $(fmt(100r.gini)) | $(fmt(100getmean("gini"))) | $(fmt(100getmean("lost_link_fraction"))) |")
        end
        println(io,"\nEach replicate samples sites uniformly without replacement, preserving their interaction and occupancy records. K and C are recomputed. Positive-support summaries use links observed in that subsample; opportunity distributions use those exact same surviving links. Original known links disappearing from the sample are explicitly counted, including those with remaining co-occurrence. This separates apparent positive-distribution shape from loss of link identities. Three hundred permutations per dataset provide intervals at 8 sites and a sequence of other sample sizes; each dataset's full-site result reproduces the input. Nahuel already has eight sites and has no eight-site subsampling variation.\n\nMatching site count tests one source of sampling dependence. It does not equalise within-site effort, detection, spatial extent, spatial clustering or seasonal coverage. Eight sites provide only eight possible support values and cannot establish a tail family. Pair/Common are evaluated on full datasets; they are not refitted or reconditioned in the subsampling exercise. Comparing subsampled positive links changes the cohort, so persistence of a qualitative shape is evidence of descriptive recurrence only.\n")
        println(io,"## Files and limitations\n\nAll distributions and concentration curves, exact expected model frequencies, conditional group summaries, allocation draws and site-subsampling draws are exported. Gini and top-10% support share use fractional-rank interpolation; model concentration curves average curves calculated separately for each allocation, not curves of expected individual support. No power-law fitting, new predictive model or ecological mechanism test is added. Shared sites, related collections, interaction-derived occupancy and unresolved data flags are retained in README and metadata. All source versions are hashed and verified.")
    end
    CSV.write(joinpath(DEST,"code_versions.csv"),DataFrame(file=["run_analysis.jl","figures_and_findings.jl"],sha256=[bytes2hex(sha256(read(joinpath(@__DIR__,f)))) for f in ["run_analysis.jl","figures_and_findings.jl"]]))
end

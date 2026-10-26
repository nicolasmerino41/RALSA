using CSV, DataFrames, Statistics, Random, LinearAlgebra, Printf, CairoMakie
include(joinpath(@__DIR__, "..", "src", "data_loading.jl"))
include(joinpath(ROOT, "src", "site_removal.jl"))

const DEST = joinpath(OUT, "sampling_sensitivity")
const DRAWS = 300
const SIZES = [5, 8, 14, 20, 30, 50, 100, 200]
const H = collect(0.0:0.1:0.8)
mkpath(DEST)
avgarea(y) = sum(diff(H) .* (y[1:end-1] .+ y[2:end]) ./ 2) / last(H)
finiteavg(x) = isempty(filter(isfinite,x)) ? NaN : mean(filter(isfinite,x))

# Mix neighbouring removal counts to match the same expected fraction at small N.
function lossgrid(N)
    [begin
        m = h*N; a = floor(Int,m); b = ceil(Int,m)
        (1-(m-a))*loss(k,N,a) + (m-a)*loss(k,N,b)
    end for h in H, k in 1:N]
end

function group_ids(deg, nodes)
    groups = Dict{Tuple{String,String},Vector{Int}}()
    for guild in ["consumer","resource"]
        ids = findall(i -> nodes[i][1]==guild && deg[i]>0, eachindex(nodes))
        vals = sort(unique(deg[ids]))
        length(vals)<2 && continue
        cut = vals[fld(length(vals),2)]
        groups[(guild,"Specialists")] = filter(i -> deg[i]<=cut, ids)
        groups[(guild,"Generalists")] = filter(i -> deg[i]>cut, ids)
    end
    groups
end

function run_sampling()
    rows, grows, checks = NamedTuple[], NamedTuple[], NamedTuple[]
    original = table("conversion_and_pink_area")
    for (di,ds) in enumerate(DATA)
        d = load_data(ds); N=d.N; L=length(d.links)
        si=Dict(s=>i for (i,s) in enumerate(d.sites))
        li=Dict(p=>i for (i,p) in enumerate(d.links))
        co_links=[Int[] for _ in 1:N]; opportunities=zeros(Int,N)
        for r in eachrow(d.co)
            i=si[r.site]; opportunities[i]+=1
            j=get(li,(r.consumer,r.resource),0)
            j>0 && push!(co_links[i],j)
        end
        nodes=vcat([("consumer",s) for s in sort(unique(d.ints.consumer))],
                   [("resource",s) for s in sort(unique(d.ints.resource))])
        ni=Dict(s=>i for (i,s) in enumerate(nodes))
        ends=[(ni[("consumer",p[1])],ni[("resource",p[2])]) for p in d.links]
        full_degree=zeros(Int,length(nodes))
        for (a,b) in ends; full_degree[a]+=1;full_degree[b]+=1;end
        full_groups=group_ids(full_degree,nodes)
        node_sites=[Set{Int}() for _ in nodes]
        sizes=sort(unique(vcat(filter(<=(N),SIZES),N)))
        grids=Dict(n=>lossgrid(n) for n in sizes)
        pink0=[mean(loss(k,N,removed_count(N,h)) for k in d.K)-
               mean(loss(k,N,removed_count(N,h)) for k in d.n) for h in H]
        @assert isapprox(avgarea(pink0),only(original.mean_pink[original.dataset.==ds]);atol=1e-12)
        push!(checks,(dataset=ds,full_sites=N,full_links=L,original_pink_verified=true))
        rng=MersenneTwister(90100+di)
        for rep in 1:DRAWS
            order=randperm(rng,N); K=zeros(Int,L); n=zeros(Int,L); co_total=0; previous=0
            for s in node_sites;empty!(s);end
            for size in sizes
                for pos in previous+1:size
                    site=order[pos];co_total+=opportunities[site]
                    for j in d.site_links[site]
                        K[j]+=1
                        a,b=ends[j];push!(node_sites[a],site);push!(node_sites[b],site)
                    end
                    for j in co_links[site];n[j]+=1;end
                end
                previous=size
                size==N && rep>1 && continue
                observed=findall(>(0),K); isempty(observed) && error("Empty network: $ds $size")
                q=grids[size]
                regional=[mean(1-q[t,K[j]] for j in observed) for t in eachindex(H)]
                pink=[mean(q[t,K[j]]-q[t,n[j]] for j in observed) for t in eachindex(H)]
                @assert all(pink.>=-1e-12) && all(regional .>= 1 .- H .- 1e-12)
                push!(rows,(dataset=ds,full_sites=N,sampled_sites=size,replicate=rep,
                    links=length(observed),mean_K=mean(K[observed]),prevalence=mean(K[observed])/size,
                    p_emp=sum(K)/co_total,mean_gap=avgarea(regional .- (1 .- H)),mean_pink=avgarea(pink),
                    gap80=regional[end]-.2,pink80=pink[end]))
                degree=zeros(Int,length(nodes));support=zeros(length(nodes));lost=zeros(length(nodes))
                for j in observed
                    for a in ends[j];degree[a]+=1;support[a]+=K[j];lost[a]+=q[end,K[j]];end
                end
                active=[isempty(s) ? 0.0 : 1-q[end,length(s)] for s in node_sites]
                perpartner=[degree[i]>0 ? support[i]/degree[i] : NaN for i in eachindex(nodes)]
                for (definition,groups) in [("refitted",group_ids(degree,nodes)),("fixed",full_groups)]
                    vals=Dict{String,NTuple{3,Float64}}()
                    for group in GROUPS
                        gs=Float64[];ga=Float64[];gl=Float64[]
                        for guild in ["consumer","resource"]
                            ids=filter(i->degree[i]>0,get(groups,(guild,group),Int[]))
                            isempty(ids) && continue
                            push!(gs,mean(perpartner[ids]));push!(ga,mean(active[ids]));push!(gl,mean(lost[ids]))
                        end
                        vals[group]=length(gs)==2 ? (mean(gs),mean(ga),mean(gl)) : (NaN,NaN,NaN)
                    end
                    a,b=vals["Generalists"],vals["Specialists"]
                    push!(grows,(dataset=ds,sampled_sites=size,replicate=rep,definition=definition,
                        support_difference=a[1]-b[1],active_difference80=a[2]-b[2],loss_difference80=a[3]-b[3]))
                end
            end
        end
        println("Sampling complete: ",ds)
    end
    r,g=DataFrame(rows),DataFrame(grows)
    CSV.write(joinpath(DEST,"replicates.csv"),r)
    CSV.write(joinpath(DEST,"group_replicates.csv"),g)
    CSV.write(joinpath(DEST,"validation.csv"),DataFrame(checks))
end

function summarise()
    r=CSV.read(joinpath(DEST,"replicates.csv"),DataFrame)
    g=CSV.read(joinpath(DEST,"group_replicates.csv"),DataFrame)
    metrics=[:links,:mean_K,:prevalence,:p_emp,:mean_gap,:mean_pink,:gap80,:pink80]
    s=combine(groupby(r,[:dataset,:full_sites,:sampled_sites]),
        [v=>mean=>v for v in metrics]...,
        [v=>(x->quantile(x,.025))=>Symbol(v,"_low") for v in [:mean_gap,:mean_pink]]...,
        [v=>(x->quantile(x,.975))=>Symbol(v,"_high") for v in [:mean_gap,:mean_pink]]...)
    full=s[s.sampled_sites.==s.full_sites,:]
    correlations=NamedTuple[]
    for size in [5,8,14,20,30]
        d=s[s.sampled_sites.==size,:];nrow(d)<3 && continue
        b=full[indexin(d.dataset,full.dataset),:]
        for metric in [:mean_gap,:mean_pink]
            rank(x)=invperm(sortperm(x))
            push!(correlations,(sampled_sites=size,datasets=nrow(d),metric=String(metric),
                full_mean=mean(b[!,metric]),standardised_mean=mean(d[!,metric]),
                mean_change=mean(d[!,metric].-b[!,metric]),
                rank_correlation=cor(rank(d[!,metric]),rank(b[!,metric])),
                between_dataset_sd=std(d[!,metric]),
                p_emp_correlation=cor(d.p_emp,d[!,metric])))
        end
    end
    CSV.write(joinpath(DEST,"summary.csv"),s)
    CSV.write(joinpath(DEST,"comparisons.csv"),DataFrame(correlations))
    gs=combine(groupby(g,[:dataset,:sampled_sites,:definition]),
        [v=>finiteavg=>v for v in [:support_difference,:active_difference80,:loss_difference80]]...,
        :support_difference=>(x->count(isfinite,x))=>:valid_replicates,
        :support_difference=>(x->finiteavg([isfinite(v) ? Float64(v>0) : NaN for v in x]))=>:support_positive_fraction)
    CSV.write(joinpath(DEST,"group_summary.csv"),gs)
    # Distribution of cross-dataset correlations across independent rarefaction draws.
    cr=NamedTuple[]
    for size in [5,8,14,20], rep in 1:DRAWS
        vals=NamedTuple[]
        for ds in DATA
            d=r[(r.dataset.==ds).&(r.sampled_sites.==size),:]
            isempty(d) && continue
            v=d[min(rep,nrow(d)),:]
            push!(vals,(p=v.p_emp,gap=v.mean_gap,pink=v.mean_pink))
        end
        length(vals)<3 && continue
        push!(cr,(sampled_sites=size,replicate=rep,datasets=length(vals),
            p_pink=cor(getproperty.(vals,:p),getproperty.(vals,:pink))))
    end
    CSV.write(joinpath(DEST,"correlation_replicates.csv"),DataFrame(cr))
    s,gs,full
end

function plots(s,gs,full)
    set_theme!(Theme(font="Arial",fontsize=19,Axis=(topspinevisible=false,rightspinevisible=false,
        xgridvisible=false,ygridvisible=false)))
    colors=["#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77","#CC6677","#882255","#AA4499","#666666"]
    fig=Figure(size=(1650,1100))
    for (col,metric,label) in [(1,:mean_gap,"Mean regional–local gap"),(2,:mean_pink,"Mean co-occurrence-only fraction")]
        ax=Axis(fig[1,col],xlabel="Number of sampled sites",ylabel=label,xscale=log10,
            xticks=([5,8,14,30,100,341],string.([5,8,14,30,100,341])))
        for (ds,color) in zip(DATA,colors)
            d=sort(s[s.dataset.==ds,:],:sampled_sites)
            band!(ax,d.sampled_sites,d[!,Symbol(metric,"_low")],d[!,Symbol(metric,"_high")],color=(color,.12))
            lines!(ax,d.sampled_sites,d[!,metric],color=color,linewidth=2)
            scatter!(ax,d.sampled_sites,d[!,metric],color=color,markersize=8)
        end
        ax2=Axis(fig[2,col],xlabel="Dataset",ylabel=label,xticks=(1:10,[DISPLAY[x] for x in DATA]),xticklabelrotation=pi/4)
        for (offset,size,color) in [(-.20,0,"#333333"),(0.,5,"#4477AA"),(.20,8,"#CC6677")]
            for (i,ds) in enumerate(DATA)
                d=size==0 ? full[full.dataset.==ds,:] : s[(s.dataset.==ds).&(s.sampled_sites.==size),:]
                rangebars!(ax2,[i+offset],d[!,Symbol(metric,"_low")],d[!,Symbol(metric,"_high")],color=color)
                scatter!(ax2,[i+offset],d[!,metric],color=color,markersize=10)
            end
        end
    end
    Legend(fig[3,1:2],[MarkerElement(color=c,marker=:circle) for c in ["#333333","#4477AA","#CC6677"]],
        ["Full sample","5 sites","8 sites"],orientation=:horizontal,framevisible=false)
    Legend(fig[4,1:2],[LineElement(color=c) for c in colors],[DISPLAY[x] for x in DATA],nbanks=2,orientation=:horizontal,framevisible=false)
    colgap!(fig.layout,50)
    save(joinpath(DEST,"sampling_effects.png"),fig;px_per_unit=2)
    f=Figure(size=(1550,650))
    for (j,metric,label) in [(1,:support_difference,"Generalist − specialist support"),
                           (2,:active_difference80,"Difference in fraction retaining interactions"),
                           (3,:loss_difference80,"Difference in absolute partners lost")]
        ax=Axis(f[1,j],xlabel="Number of sampled sites",ylabel=label,xscale=log10,
            xticks=([5,8,30,100,341],string.([5,8,30,100,341])))
        hlines!(ax,[0],color=:gray60,linestyle=:dash)
        for (ds,color) in zip(DATA,colors)
            d=sort(gs[(gs.dataset.==ds).&(gs.definition.=="refitted"),:],:sampled_sites)
            lines!(ax,d.sampled_sites,d[!,metric],color=color,linewidth=2)
        end
    end
    Legend(f[2,1:3],[LineElement(color=c) for c in colors],[DISPLAY[x] for x in DATA],nbanks=2,orientation=:horizontal,framevisible=false)
    save(joinpath(DEST,"sampling_species_groups.png"),f;px_per_unit=2)
end

!("--plot-only" in ARGS) && run_sampling()
plots(summarise()...)
println("Saved sampling sensitivity to ",DEST)

function generalists_and_specialists()
    out = OUT
    baseline = NamedTuple[]
    repout = NamedTuple[]
    groupinfo = NamedTuple[]
    for (di, ds) in enumerate(DATA)
        d = load_data(ds)
        ints, sites, links, N = d.ints, d.sites, d.links, d.N
        L = length(links)
        nodes = vcat([("Consumer", s) for s in sort(unique(ints.consumer))], [("Resource", s) for s in sort(unique(ints.resource))])
        ni = Dict(s => i for (i, s) in enumerate(nodes))
        ends = [(ni[("Consumer", p[1])], ni[("Resource", p[2])]) for p in links]
        site_links, K = d.site_links, d.K
        degree = zeros(Int, length(nodes))
        support = zeros(length(nodes))
        for (j, (a, b)) in enumerate(ends)
            degree[a] += 1
            degree[b] += 1
            support[a] += K[j]
            support[b] += K[j]
        end
        support ./= degree
        memberships = Dict{ Tuple{String,String}, Vector{Int} }()
        for guild in ["Consumer", "Resource"]
            ids = findall(x -> x[1] == guild, nodes)
            values = sort(unique(degree[ids]))
            @assert length(values) > 1
            # Split distinct degree values, not species counts.
            cut = values[fld(length(values), 2)]
            for (group, predicate) in [("Specialists", k -> k <= cut), ("Generalists", k -> k > cut)]
                idx = filter(i -> predicate(degree[i]), ids)
                @assert !isempty(idx)
                memberships[(guild, group)] = idx
                push!(baseline, (dataset=ds, guild=guild, group=group, support=mean(support[idx])))
                push!(groupinfo, (dataset=ds, guild=guild, group=group, n_species=length(idx), minimum_degree=minimum(degree[idx]), maximum_degree=maximum(degree[idx])))
            end
        end
        rng = MersenneTwister(SEED + di)
        for rep in 1:NREPS
            order = randperm(rng, N)
            remaining = copy(K)
            last = 0
            retained = copy(degree)
            for f in LEVELS
                m = removed_count(N, f)
                for pos in last+1:m
                    for j in site_links[order[pos]]
                        remaining[j] -= 1
                        if remaining[j] == 0
                            a, b = ends[j]
                            retained[a] -= 1
                            retained[b] -= 1
                        end
                    end
                end
                last = m
                for ((guild, group), ids) in memberships
                    lost = mean(degree[ids] .- retained[ids])
                    active = mean(retained[ids] .> 0)
                    @assert lost >= 0
                    @assert 0 <= active <= 1
                    if f == 0
                        @assert lost == 0
                        @assert active == 1
                    end
                    push!(repout, (dataset=ds, guild=guild, group=group, removal=f, replicate=rep, links_lost=lost, active_fraction=active))
                end
            end
        end
        println("Group analysis: ", ds)
    end
    a = DataFrame(baseline)
    r = DataFrame(repout)
    ag = combine(groupby(a, [:dataset, :group]), :support => mean => :support)
    rg = combine(groupby(r, [:dataset, :group, :removal, :replicate]), :links_lost => mean => :links_lost, :active_fraction => mean => :active_fraction)
    rs = combine(groupby(rg, [:dataset, :group, :removal]), :links_lost => mean => :links_lost, :active_fraction => mean => :active_fraction)
    for (name, t) in [("support_by_guild", a), ("support_merged", ag), ("group_definitions", DataFrame(groupinfo)), ("group_replicates", r), ("group_curves", rs)]
        CSV.write(joinpath(out, "tables", name * ".csv"), t)
    end
    ag, rs
end

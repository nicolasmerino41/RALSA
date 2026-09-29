function degree_distribution()
    d = load_data(EXAMPLE)
    consumers = sort(unique(d.ints.consumer))
    ci = Dict(s => i for (i,s) in enumerate(consumers))
    link_consumer = [ci[p[1]] for p in d.links]
    initial = [count(==(i), link_consumer) for i in eachindex(consumers)]
    samples = Dict(f => Vector{Vector{Float64}}() for f in [0.,.4,.8])
    rng = MersenneTwister(SEED + findfirst(==(EXAMPLE), DATA))
    for rep in 1:NREPS
        order = randperm(rng,d.N)
        remaining, last_removed = copy(d.K), 0
        for f in REMOVAL
            m = removed_count(d.N,f)
            for pos in last_removed+1:m, j in d.site_links[order[pos]]
                remaining[j] -= 1
            end
            last_removed = m
            deg = zeros(Int,length(consumers))
            for j in eachindex(remaining)
                remaining[j] > 0 && (deg[link_consumer[j]] += 1)
            end
            active = deg .> 0
            if haskey(samples,f) && any(active)
                push!(samples[f], [mean(deg[active] .>= k) for k in 1:maximum(initial)])
            end
            f == 0 && @assert deg == initial
        end
    end
    rows = NamedTuple[]
    for f in [0.,.4,.8]
        mat = reduce(hcat,samples[f])
        for k in axes(mat,1)
            v = mat[k,:]
            push!(rows,(dataset=EXAMPLE,removal=f,degree=k,probability=median(v),
                        low=quantile(v,.025),high=quantile(v,.975)))
        end
    end
    save_table("degree_example",DataFrame(rows))
end

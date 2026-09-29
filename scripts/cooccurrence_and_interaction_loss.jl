function cooccurrence_and_interaction_loss()
    states, summaries = NamedTuple[], NamedTuple[]
    for ds in DATA
        d = load_data(ds)
        pink = Float64[]
        for f in REMOVAL
            m = removed_count(d.N, f)
            grey = mean(loss(k,d.N,m) for k in d.n)
            extinct = mean(loss(k,d.N,m) for k in d.K)
            push!(pink, extinct-grey)
            @assert extinct-grey >= -1e-12
            push!(states, (dataset=ds, removal=f, green=1-extinct, pink=extinct-grey, grey=grey))
        end
        area = sum(diff(REMOVAL) .* (pink[1:end-1] .+ pink[2:end]) ./ 2)
        push!(summaries, (dataset=ds, sites=d.N, links=length(d.links),
            p_emp=nrow(d.ints)/nrow(d.co), f_regional=length(d.links)/nrow(unique(d.co[:,[:consumer,:resource]])),
            pink_area=area, mean_pink=area/last(REMOVAL), mean_cooccurrence_support=mean(d.n)))
    end
    c = innerjoin(table("retention"), DataFrame(states); on=[:dataset,:removal])
    sort!(c, [:dataset,:removal])
    @assert all(isapprox.(c.green .+ c.pink .+ c.grey, 1))
    save_table("retention_and_states", c)
    save_table("conversion_and_pink_area", DataFrame(summaries))
    pooled = combine(groupby(c,:removal), [v => mean => v for v in
                     [:species,:regional,:local_support,:green,:pink,:grey]]...)
    save_table("pooled_retention", sort!(pooled,:removal))
end

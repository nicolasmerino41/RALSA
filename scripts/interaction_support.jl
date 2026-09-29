function interaction_support()
    support, retention = NamedTuple[], NamedTuple[]
    for ds in DATA
        d = load_data(ds)
        for (j, pair) in enumerate(d.links)
            push!(support, (dataset=ds, consumer=pair[1], resource=pair[2], support=d.K[j],
                            cooccurrence_support=d.n[j], sites=d.N))
        end
        for f in REMOVAL
            m = removed_count(d.N, f)
            push!(retention, (dataset=ds, removal=f, actual_removal=m/d.N,
                species=mean(1-loss(k,d.N,m) for k in d.species_support),
                regional=mean(1-loss(k,d.N,m) for k in d.K), local_support=1-m/d.N))
        end
    end
    save_table("link_support", DataFrame(support))
    save_table("retention", DataFrame(retention))
end

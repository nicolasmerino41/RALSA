function prepare_data()
    mkpath(INPUT)
    checks = NamedTuple[]
    expected = Dict(row.file => row.sha256 for row in
                    eachrow(CSV.read(joinpath(ROOT, "data", "checksums.csv"), DataFrame)))
    for ds in DATA
        inputs = Dict{String,DataFrame}()
        for kind in ["interactions", "cooccurrences", "occupancy"]
            file = ds * "_" * kind * ".csv"
            path = joinpath(RAW, file)
            @assert bytes2hex(sha256(read(path))) == expected[file] "Input checksum changed: $file"
            df = CSV.read(path, DataFrame; types=String)
            @assert !any(ismissing, Matrix(df)) && all(x -> !isempty(strip(x)), Matrix(df))
            @assert nrow(unique(df)) == nrow(df) "Duplicate records in $file"
            inputs[kind] = df
            CSV.write(joinpath(INPUT, file), df)
        end
        ints, co, occ = inputs["interactions"], inputs["cooccurrences"], inputs["occupancy"]
        @assert issubset(Set(Tuple.(eachrow(ints))), Set(Tuple.(eachrow(co))))
        occupied = Set((r.site, r.species, r.trophic_level) for r in eachrow(occ))
        for r in eachrow(ints)
            @assert (r.site, r.consumer, "consumer") in occupied
            @assert (r.site, r.resource, "resource") in occupied
        end
        d = load_data(ds)
        @assert all(1 .<= d.K .<= d.n .<= d.N)
        if ds == "Gottin_HP"
            @assert length(d.links) == 108 && !("V1" in ints.consumer)
        end
        push!(checks, (dataset=ds, sites=d.N, links=length(d.links),
                       consumers=length(unique(ints.consumer)), checks_passed=true))
    end
    save_table("validation", DataFrame(checks))
end

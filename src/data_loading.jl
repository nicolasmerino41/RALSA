using CSV, DataFrames, Statistics, Random, LinearAlgebra, Printf, SHA
include("config.jl")

table(name) = CSV.read(joinpath(OUT, "tables", name * ".csv"), DataFrame)
function save_table(name, df)
    mkpath(joinpath(OUT, "tables"))
    CSV.write(joinpath(OUT, "tables", name * ".csv"), df)
end

function load_data(ds)
    read_input(kind) = CSV.read(joinpath(INPUT, ds * "_" * kind * ".csv"), DataFrame; types=String)
    ints, co, occ = read_input("interactions"), read_input("cooccurrences"), read_input("occupancy")
    sites = sort(unique(vcat(ints.site, co.site, occ.site)))
    links = unique(collect(zip(ints.consumer, ints.resource)))
    li = Dict(p => i for (i, p) in enumerate(links))
    si = Dict(s => i for (i, s) in enumerate(sites))
    site_links = [Int[] for _ in sites]
    K = zeros(Int, length(links))
    for r in eachrow(ints)
        j = li[(r.consumer, r.resource)]
        push!(site_links[si[r.site]], j)
        K[j] += 1
    end
    counts = Dict((r.consumer, r.resource) => r.count for r in
                  eachrow(combine(groupby(co, [:consumer, :resource]), nrow => :count)))
    n = [counts[p] for p in links]
    species_support = combine(groupby(occ, [:trophic_level, :species]), nrow => :support).support
    (; ints, co, occ, sites, links, site_links, K, n, species_support, N=length(sites))
end

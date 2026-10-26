module EmpiricalSupportScatter

using CSV, DataFrames, Statistics, SHA, CairoMakie

const ROOT = isempty(ARGS) ? normpath(joinpath(@__DIR__, "../..")) : abspath(ARGS[1])
const DEST = length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(ROOT, "outputs/statistical description of local support/empirical_support_vs_cooccurrence")
const INPUT = joinpath(ROOT, "outputs/support_prediction_zero/confirmed")

const DATA = [
    "Quercus", "Nahuel", "Salix_Galpar", "Gottin_HP", "Gottin_PP",
    "Garraf_HP", "Garraf_PP", "Garraf_PP2", "Olot", "Montseny"
]

const COLORS = [
    "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
    "#DDCC77", "#CC6677", "#882255", "#AA4499", "#661100"
]

const CASES = [
    (true,  :K,        "01_absolute_K_including_zeros", "All co-occurring pairs"),
    (false, :K,        "02_absolute_K_positive_only",   "Observed regional links"),
    (true,  :fraction, "03_fraction_including_zeros",   "All co-occurring pairs"),
    (false, :fraction, "04_fraction_positive_only",    "Observed regional links")
]

mkpath(DEST)

versions = CSV.read(joinpath(INPUT, "input_versions.csv"), DataFrame)
provenance = DataFrame(path=String[], sha256=String[])

for r in eachrow(versions)
    path = joinpath(INPUT, "inputs", r.dataset * "_" * r.kind * ".csv")
    digest = bytes2hex(sha256(read(path)))
    @assert digest == r.sha256
    push!(provenance, (relpath(path, ROOT), digest))
end

pairs = CSV.read(joinpath(INPUT, "pairs.csv"), DataFrame)

push!(
    provenance,
    (
        relpath(joinpath(INPUT, "pairs.csv"), ROOT),
        bytes2hex(sha256(read(joinpath(INPUT, "pairs.csv"))))
    )
)

qcpath = joinpath(
    ROOT,
    "outputs/statistical description of local support/opportunity_theory/quality_control.csv"
)

qc = CSV.read(qcpath, DataFrame)

push!(
    provenance,
    (relpath(qcpath, ROOT), bytes2hex(sha256(read(qcpath))))
)

checks = NamedTuple[]

for ds in DATA

    co = CSV.read(
        joinpath(INPUT, "inputs", ds * "_cooccurrences.csv"),
        DataFrame;
        types=String
    )

    ints = CSV.read(
        joinpath(INPUT, "inputs", ds * "_interactions.csv"),
        DataFrame;
        types=String
    )

    n = Dict(
        (g.consumer[1], g.resource[1]) => length(unique(g.site))
        for g in groupby(co, [:consumer, :resource])
    )

    k = Dict(
        (g.consumer[1], g.resource[1]) => length(unique(g.site))
        for g in groupby(ints, [:consumer, :resource])
    )

    v = pairs[pairs.dataset .== ds, :]

    @assert length(n) == nrow(v)

    for r in eachrow(v)
        @assert r.n == n[(r.consumer, r.resource)] &&
                r.K == get(k, (r.consumer, r.resource), 0)
    end

    q = qc[qc.dataset .== ds, :]

    @assert Dict(
        (r.consumer, r.resource) => (r.K, r.C, r.N)
        for r in eachrow(q)
    ) == Dict(
        (r.consumer, r.resource) => (r.K, r.n, r.N)
        for r in eachrow(v) if r.K > 0
    )

    @assert all(0 .<= v.K .<= v.n .<= v.N) && all(v.n .>= 1)

    push!(
        checks,
        (
            dataset=ds,
            N=only(unique(v.N)),
            cooccurring_pairs=nrow(v),
            positive_links=count(v.K .> 0),
            zero_support_pairs=count(v.K .== 0),
            max_n=maximum(v.n),
            max_K=maximum(v.K),
            records_verified=true
        )
    )
end

pairs.fraction = pairs.K ./ pairs.n

CSV.write(joinpath(DEST, "empirical_pairs.csv"), pairs)
CSV.write(joinpath(DEST, "input_versions.csv"), provenance)
CSV.write(joinpath(DEST, "dataset_counts_and_checks.csv"), DataFrame(checks))

cells = combine(
    groupby(pairs, [:dataset, :n, :K]),
    nrow => :pair_count
)

cells.fraction = cells.K ./ cells.n

CSV.write(joinpath(DEST, "plotted_cells.csv"), cells)

# Retain empirical means in CSV, but do not plot them.
means = NamedTuple[]

for (zeros, field, _, _) in CASES, ds in DATA
    v = pairs[(pairs.dataset .== ds) .& (zeros .| (pairs.K .> 0)), :]

    for g in groupby(v, :n)
        push!(
            means,
            (
                dataset=ds,
                include_zero_K=zeros,
                quantity=String(field),
                n=g.n[1],
                pairs=nrow(g),
                mean=mean(g[!, field])
            )
        )
    end
end

meanframe = DataFrame(means)
CSV.write(joinpath(DEST, "empirical_means.csv"), meanframe)

set_theme!(
    Theme(
        font="Arial",
        fontsize=20,
        Axis=(
            topspinevisible=false,
            rightspinevisible=false,
            xgridvisible=false,
            ygridvisible=false,
            xticklabelsize=17,
            yticklabelsize=17
        )
    )
)

# Scatter points only: no mean lines or fitted curves.
function drawpanel(ax, ds, zeros, field, col; combined=false)

    v = cells[
        (cells.dataset .== ds) .&
        (zeros .| (cells.K .> 0)),
        :
    ]

    # Coincident pairs share one marker.
    # Multiplicities remain available in plotted_cells.csv.
    scatter!(
        ax,
        v.n,
        v[!, field];
        color=(col, combined ? 0.20 : 0.40),
        markersize=combined ? 5 : 7
    )
end

# Individual dataset figures
for (zeros, field, name, label) in CASES

    fig = Figure(size=(1750, 850))

    for (i, ds) in enumerate(DATA)

        v = pairs[pairs.dataset .== ds, :]

        ax = Axis(
            fig[cld(i, 5), mod1(i, 5)];
            title=replace(ds, "_" => " "),
            xlabel="Co-occurring sites, n",
            ylabel=mod1(i, 5) == 1 ?
                (field == :K ? "Local support, K" : "Realisation fraction, K/n") :
                ""
        )

        drawpanel(ax, ds, zeros, field, COLORS[i])

        xlims!(ax, 0.5, maximum(v.n) + 0.5)

        ylims!(
            ax,
            field == :K ?
                (-0.5, maximum(v.K) + 0.8) :
                (-0.04, 1.04)
        )

        field == :fraction && (ax.yticks = 0:0.25:1)
    end

    save(
        joinpath(DEST, name * ".png"),
        fig;
        px_per_unit=2
    )
end

# Combined figure: four cases, ten datasets
fig = Figure(size=(1500, 1120))

xt = filter(
    x -> x <= maximum(pairs.n),
    [1, 2, 5, 10, 20, 50, 100, 200, 500]
)

for (zeros, field, _, label) in CASES

    row = field == :K ? 1 : 2
    col = zeros ? 1 : 2

    ax = Axis(
        fig[row, col];
        title=label,
        xlabel="Co-occurring sites, n (log scale)",
        ylabel=field == :K ?
            "Local support, K" :
            "Realisation fraction, K/n",
        xscale=log10,
        xticks=(xt, string.(xt))
    )

    for (i, ds) in enumerate(DATA)
        drawpanel(
            ax,
            ds,
            zeros,
            field,
            COLORS[i];
            combined=true
        )
    end

    xlims!(ax, 0.85, maximum(pairs.n) * 1.08)

    ylims!(
        ax,
        field == :K ?
            (-1.0, maximum(pairs.K) + 2.0) :
            (-0.04, 1.04)
    )

    field == :fraction && (ax.yticks = 0:0.25:1)
end

Legend(
    fig[3, 1:2],
    [MarkerElement(color=c, marker=:circle, markersize=12) for c in COLORS],
    replace.(DATA, "_" => " ");
    orientation=:horizontal,
    nbanks=2,
    framevisible=false,
    labelsize=17
)

save(
    joinpath(DEST, "05_all_datasets_four_cases.png"),
    fig;
    px_per_unit=2
)

CSV.write(
    joinpath(DEST, "code_versions.csv"),
    DataFrame(
        file=["plot_empirical.jl"],
        sha256=[bytes2hex(sha256(read(@__FILE__)))]
    )
)

println("Saved five empirical figures without lines: ", DEST)

end
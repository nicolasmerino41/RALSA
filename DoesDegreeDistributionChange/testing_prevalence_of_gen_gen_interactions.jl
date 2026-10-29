
using CSV, DataFrames, Statistics, CairoMakie, Printf

# ============================================================
# CONSUMER–RESOURCE DEGREE MIXING
# Regional links only; no local support analysis
# ============================================================

function find_root(dir)
    while true
        if isfile(joinpath(dir, "outputs",
                "support_prediction_zero", "confirmed", "pairs.csv"))
            return dir
        end
        parent = dirname(dir)
        parent == dir && error("Cannot find pairs.csv")
        dir = parent
    end
end

ROOT = find_root(@__DIR__)
INPUT = joinpath(ROOT, "outputs",
    "support_prediction_zero", "confirmed", "pairs.csv")
OUT = joinpath(ROOT, "outputs", "consumer_resource_mixing")
mkpath(OUT)

KINDS = ["S-S", "S-G", "G-S", "G-G"]
COLORS = ["#20363B", "#248524", "#9B7152", "#6378A7"]

# ============================================================
# LOAD DATA
# ============================================================

pairs = CSV.read(INPUT, DataFrame)
dropmissing!(pairs, [:dataset, :consumer, :resource, :K])
filter!(r -> r.K > 0, pairs)

pairs.dataset = string.(pairs.dataset)

@assert nrow(unique(pairs,
    [:dataset, :consumer, :resource])) == nrow(pairs)

results = DataFrame(
    dataset=String[],
    type=String[],
    links=Int[],
    proportion=Float64[],
    expected_proportion=Float64[],
    difference_pp=Float64[]
)

mixing = DataFrame(
    dataset=String[],
    total_links=Int[],
    phi=Float64[]
)

# ============================================================
# ANALYSIS
# ============================================================

for dataset in sort(unique(pairs.dataset))

    d = copy(pairs[pairs.dataset .== dataset, :])

    # Regional degrees
    dc = combine(groupby(d, :consumer), nrow => :degree_c)
    dr = combine(groupby(d, :resource), nrow => :degree_r)

    d = leftjoin(d, dc, on=:consumer)
    d = leftjoin(d, dr, on=:resource)

    # Specialist = one partner; Generalist = 2+ partners
    d.C = ifelse.(d.degree_c .== 1, "S", "G")
    d.R = ifelse.(d.degree_r .== 1, "S", "G")
    d.type = string.(d.C, "-", d.R)

    N = nrow(d)
    counts = [count(==(k), d.type) for k in KINDS]
    proportions = counts ./ N

    # Expected composition under independent endpoint mixing
    pC = (counts[1] + counts[2]) / N
    pR = (counts[1] + counts[3]) / N

    expected = [
        pC * pR,
        pC * (1 - pR),
        (1 - pC) * pR,
        (1 - pC) * (1 - pR)
    ]

    # Phi: assortativity of the two degree categories
    a, b, c, e = counts

    denominator = sqrt(
        Float64(a+b) * (c+e) * (a+c) * (b+e)
    )

    phi = denominator > 0 ?
        (a*e - b*c) / denominator : NaN

    push!(mixing, (dataset, N, phi))

    for j in 1:4
        push!(results, (
            dataset,
            KINDS[j],
            counts[j],
            proportions[j],
            expected[j],
            100 * (proportions[j] - expected[j])
        ))
    end
end

# Sort datasets by assortativity
sort!(mixing, :phi)

CSV.write(joinpath(OUT, "interaction_frequencies.csv"), results)
CSV.write(joinpath(OUT, "degree_assortativity.csv"), mixing)

# ============================================================
# COMPARATIVE FIGURE
# ============================================================

datasets = mixing.dataset
nd = length(datasets)

fig = Figure(
    size=(1200, max(450, 65*nd + 150)),
    fontsize=15,
    backgroundcolor=:white
)

ax1 = Axis(fig[1, 1],
    title="A   Composition of regional links",
    xlabel="Percentage of links",
    yticks=(1:nd, datasets),
    xticks=(0:0.25:1, ["0%", "25%", "50%", "75%", "100%"]),
    yreversed=true,
    xgridvisible=false,
    ygridvisible=false
)

ax2 = Axis(fig[1, 2],
    title="B   Degree assortativity",
    xlabel="Assortativity (ϕ)",
    yticks=(1:nd, fill("", nd)),
    xticks=-1:0.5:1,
    yreversed=true,
    xgridvisible=false,
    ygridvisible=false
)

# Panel A: stacked proportions
for (i, dataset) in enumerate(datasets)

    s = results[results.dataset .== dataset, :]
    left = 0.0

    for (j, kind) in enumerate(KINDS)

        p = s.proportion[s.type .== kind][1]
        p <= 0 && continue

        poly!(ax1,
            Rect2f(left, i - 0.32, p, 0.64),
            color=COLORS[j],
            strokecolor=:white,
            strokewidth=1
        )

        if p >= 0.075
            text!(ax1, @sprintf("%.0f%%", 100p),
                position=(left + p/2, i),
                align=(:center, :center),
                color=:white,
                fontsize=12
            )
        end

        left += p
    end
end

# Panel B: assortativity
vlines!(ax2, [0.0],
    color=:gray50,
    linestyle=:dash,
    linewidth=1.5
)

for (i, row) in enumerate(eachrow(mixing))

    isfinite(row.phi) || continue

    color = row.phi >= 0 ? "#248524" : "#B46A53"

    lines!(ax2, [0.0, row.phi], [i, i],
        color=color, linewidth=2.5)

    scatter!(ax2, [row.phi], [i],
        color=color, markersize=15)
end

xlims!(ax1, 0, 1)
xlims!(ax2, -1, 1)
ylims!(ax1, 0.5, nd + 0.5)
ylims!(ax2, 0.5, nd + 0.5)

Legend(fig[0, 1:2],
    [PolyElement(color=c) for c in COLORS],
    KINDS,
    orientation=:horizontal,
    framevisible=false,
    nbanks=4
)

colsize!(fig.layout, 1, Relative(0.65))
colgap!(fig.layout, 30)

save(joinpath(OUT, "degree_mixing_comparison.png"),
    fig, px_per_unit=2)

println("\nINTERACTION COMPOSITION")
show(stdout, MIME("text/plain"), results)
println("\n\nDEGREE ASSORTATIVITY")
show(stdout, MIME("text/plain"), mixing)
println("\n\nSaved to: ", OUT)

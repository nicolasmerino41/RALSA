using CSV, DataFrames, Statistics, CairoMakie, Printf

# ============================================================
# CONSUMER–RESOURCE DEGREE COUPLING IN LOCAL SUPPORT
#
# Input: outputs/support_prediction_zero/confirmed/pairs.csv
# Columns: dataset, consumer, resource, K, n
#
# Output: outputs/consumer_resource_coupling/
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

const ROOT = find_root(@__DIR__)
const INPUT = joinpath(ROOT, "outputs",
    "support_prediction_zero", "confirmed", "pairs.csv")
const OUT = joinpath(ROOT, "outputs", "consumer_resource_coupling")

mkpath(OUT)

pairs = CSV.read(INPUT, DataFrame)
dropmissing!(pairs, [:dataset, :consumer, :resource, :K, :n])

pairs.dataset = string.(pairs.dataset)
pairs.K = Int.(pairs.K)
pairs.n = Int.(pairs.n)

@assert all((pairs.K .>= 0) .& (pairs.K .<= pairs.n))
@assert nrow(unique(pairs,
    [:dataset, :consumer, :resource])) == nrow(pairs)

# Only observed regional links determine regional degree
filter!(r -> r.K > 0, pairs)

# ============================================================
# HELPERS
# ============================================================

safe_mean(x) = isempty(x) ? NaN : mean(x)

function cell_value(s, c, r, field)
    rows = s[(s.consumer_group .== c) .&
             (s.resource_group .== r), :]
    return getproperty(rows, field)[1]
end

function draw_panel!(fig, j, summary, field, title;
                     adjusted=false)

    Z = fill(NaN, 2, 2)
    counts = zeros(Int, 2, 2)

    for row in eachrow(summary)
        x = row.resource_group == "S" ? 1 : 2
        y = row.consumer_group == "S" ? 1 : 2

        Z[x, y] = Float64(getproperty(row, field))
        counts[x, y] = adjusted ? row.adjusted_links : row.links
    end

    vals = filter(isfinite, vec(Z))

    if adjusted
        span = isempty(vals) ? 0.2 :
            max(0.2, maximum(abs.(vals .- 1.0)))

        colorrange = (1.0 - span, 1.0 + span)
        palette = cgrad(["#7FB6D4", "#FFFFFF", "#E5A16C"])
    else
        upper = isempty(vals) ? 1.0 : max(1.0, maximum(vals))
        colorrange = (0.0, upper)
        palette = cgrad(["#F4F8FC", "#B8D9EA", "#70ABD0"])
    end

    ax = Axis(fig[1, j],
        title=title,
        xlabel="Resource degree",
        ylabel=j == 1 ? "Consumer degree" : "",
        xticks=([1, 2], ["Lower", "Higher"]),
        yticks=([1, 2], ["Lower", "Higher"]),
        titlesize=17,
        xlabelsize=14,
        ylabelsize=14
    )

    hm = heatmap!(ax, 1:2, 1:2, Z;
        colormap=palette,
        colorrange=colorrange
    )

    for x in 1:2, y in 1:2
        value = Z[x, y]

        txt = isfinite(value) ?
            @sprintf("%.2f\n(%d)", value, counts[x, y]) :
            "—"

        text!(ax, txt;
            position=(x, y),
            align=(:center, :center),
            fontsize=16,
            color=:black
        )
    end

    limits!(ax, 0.5, 2.5, 0.5, 2.5)
    ax.aspect = DataAspect()

    Colorbar(fig[2, j], hm; vertical=false)
end

# ============================================================
# ANALYSIS
# ============================================================

summary = DataFrame(
    dataset=String[],
    consumer_group=String[],
    resource_group=String[],
    links=Int[],
    mean_K=Float64[],
    mean_n=Float64[],
    mean_q=Float64[],
    relative_K=Float64[],
    adjusted_links=Int[]
)

contrasts = DataFrame(
    dataset=String[],
    specialist_consumer_partner_effect=Float64[],
    specialist_resource_partner_effect=Float64[],
    GG_interaction_contrast=Float64[]
)

for dataset in sort(unique(pairs.dataset))

    d = copy(pairs[pairs.dataset .== dataset, :])
    isempty(d) && continue

    # --------------------------------------------------------
    # 1. Regional degree of both endpoints
    # --------------------------------------------------------

    dc = combine(groupby(d, :consumer), nrow => :degree_c)
    dr = combine(groupby(d, :resource), nrow => :degree_r)

    d = leftjoin(d, dc, on=:consumer)
    d = leftjoin(d, dr, on=:resource)

    # Split lower/higher using median of distinct degrees
    cut_c = median(unique(d.degree_c))
    cut_r = median(unique(d.degree_r))

    d.consumer_group = ifelse.(d.degree_c .<= cut_c, "S", "G")
    d.resource_group = ifelse.(d.degree_r .<= cut_r, "S", "G")

    d.pair = string.(d.consumer_group, d.resource_group)

    # --------------------------------------------------------
    # 2. Opportunity-adjusted support
    #
    # Compare K with E[K | n] across regional links.
    # Use only n strata represented by at least 2 pair types.
    # --------------------------------------------------------

    baseline = combine(groupby(d, :n),
        :K => mean => :expected_K,
        :pair => (x -> length(unique(x))) => :n_groups
    )

    d = leftjoin(d, baseline, on=:n)

    d.adjusted_K = [
        row.n_groups >= 2 ?
            row.K / row.expected_K : NaN
        for row in eachrow(d)
    ]

    # --------------------------------------------------------
    # 3. Summarise four degree combinations
    # --------------------------------------------------------

    local_summary = summary[1:0, :]

    for c in ("S", "G"), r in ("S", "G")

        subset = d[(d.consumer_group .== c) .&
                   (d.resource_group .== r), :]

        adjusted_values = filter(isfinite, subset.adjusted_K)

        push!(local_summary, (
            dataset,
            c,
            r,
            nrow(subset),
            safe_mean(subset.K),
            safe_mean(subset.n),
            safe_mean(subset.K ./ subset.n),
            safe_mean(adjusted_values),
            length(adjusted_values)
        ))
    end

    append!(summary, local_summary)

    # --------------------------------------------------------
    # 4. Simple contrasts in adjusted support
    # --------------------------------------------------------

    SS = cell_value(local_summary, "S", "S", :relative_K)
    SG = cell_value(local_summary, "S", "G", :relative_K)
    GS = cell_value(local_summary, "G", "S", :relative_K)
    GG = cell_value(local_summary, "G", "G", :relative_K)

    # Does a generalist resource increase support for
    # specialist consumers?
    effect_SC = SG - SS

    # Does a generalist consumer increase support for
    # specialist resources?
    effect_SR = GS - SS

    # Non-additivity of the two degree effects
    interaction = GG - GS - SG + SS

    push!(contrasts, (
        dataset, effect_SC, effect_SR, interaction
    ))

    # --------------------------------------------------------
    # 5. Figure
    # --------------------------------------------------------

    fig = Figure(size=(1200, 460), fontsize=15)

    Label(fig[0, 1:3], dataset;
        fontsize=23, font=:bold)

    draw_panel!(fig, 1, local_summary,
        :mean_K, "A   Mean support (K)")

    draw_panel!(fig, 2, local_summary,
        :mean_n, "B   Mean co-occurrence (n)")

    draw_panel!(fig, 3, local_summary,
        :relative_K, "C   Support relative to E[K | n]";
        adjusted=true)

    Label(fig[3, 1:3],
        "Numbers in parentheses = links contributing to each estimate. " *
        "S/G = lower/higher regional degree.";
        fontsize=12)

    filename = replace(dataset, r"[^A-Za-z0-9_-]" => "_")
    save(joinpath(OUT, "$(filename)_coupling.png"), fig)

    println("\n", dataset,
        " | Links: ", nrow(d),
        " | Consumer cutoff: ", cut_c,
        " | Resource cutoff: ", cut_r)

    @printf("  Specialist consumer, generalist resource: %+.3f\n",
        effect_SC)
    @printf("  Specialist resource, generalist consumer: %+.3f\n",
        effect_SR)
    @printf("  G-G non-additivity:                     %+.3f\n",
        interaction)
end

# ============================================================
# EXPORT
# ============================================================

CSV.write(joinpath(OUT, "coupling_summary.csv"), summary)
CSV.write(joinpath(OUT, "coupling_contrasts.csv"), contrasts)

println("\nSaved all results to: ", OUT)

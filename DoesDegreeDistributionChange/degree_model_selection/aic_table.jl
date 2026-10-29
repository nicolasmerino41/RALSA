
using CSV, DataFrames, Statistics, CairoMakie, Printf, SHA

folder = length(ARGS) > 0 ? ARGS[1] : @__DIR__
out = joinpath(folder, "outputs")
input = joinpath(out, "fit_diagnostics.csv")

t = CSV.read(input, DataFrame)

models = [
    "Power law",
    "Truncated power law",
    "Log-normal",
    "Exponential"
]

labels = [
    "Power law",
    "Truncated\npower law",
    "Log-normal",
    "Exponential"
]

fractions = [0.0, 0.4, 0.8]
datasets = unique(t.dataset)

rows = NamedTuple[]

# ------------------------------------------------------------
# Summarise AIC and AIC weights
# ------------------------------------------------------------

for ds in datasets, f in fractions, model in models

    v = t[
        (t.dataset .== ds) .&
        (t.removal .== f) .&
        (t.model .== model),
        :
    ]

    @assert nrow(v) == (f == 0 ? 1 : 100)

    push!(rows, (
        community = v.community[1],
        dataset = ds,
        removal = f,
        model = model,
        AIC = mean(v.AIC),
        AIC_sd = nrow(v) > 1 ? std(v.AIC) : 0.0,
        mean_AIC_weight = mean(v.AIC_weight),
        sequences = nrow(v)
    ))
end

s = DataFrame(rows)

# ------------------------------------------------------------
# Identify winning model for each dataset/removal
# Criterion: highest mean AIC weight
# Mark exact or numerically equivalent ties.
# ------------------------------------------------------------

s.winner = falses(nrow(s))

for ds in datasets, f in fractions

    idx = findall(
        (s.dataset .== ds) .&
        (s.removal .== f)
    )

    weights = s.mean_AIC_weight[idx]
    best = maximum(weights)

    for i in idx
        s.winner[i] = isapprox(
            s.mean_AIC_weight[i],
            best;
            atol = 1e-10,
            rtol = 0
        )
    end
end

CSV.write(
    joinpath(out, "03_all_models_AIC.csv"),
    s
)

# ------------------------------------------------------------
# Figure
# ------------------------------------------------------------

set_theme!(Theme(font = "Arial", fontsize = 17))

fig = Figure(size = (1900, 750))

winner_color = "#0072B2"
normal_color = "#303030"

Label(
    fig[0, 1:13],
    "AIC by model and removal level",
    font = :bold,
    fontsize = 24
)

Label(
    fig[1:2, 1],
    "Community",
    font = :bold,
    fontsize = 19,
    halign = :left
)

# Column headings
for (j, f) in enumerate(fractions)

    cols = (2 + 4 * (j - 1)):(5 + 4 * (j - 1))

    Label(
        fig[1, cols],
        "$(Int(100f))% removal",
        font = :bold,
        fontsize = 21
    )

    for mi in 1:4
        Label(
            fig[2, cols[mi]],
            labels[mi],
            fontsize = 16
        )
    end
end

# Values and highlighted winners
for (i, ds) in enumerate(datasets)

    Label(
        fig[i + 2, 1],
        "C$i  " * replace(ds, "_" => " "),
        fontsize = 17,
        halign = :left
    )

    for (j, f) in enumerate(fractions), mi in 1:4

        r = only(eachrow(s[
            (s.dataset .== ds) .&
            (s.removal .== f) .&
            (s.model .== models[mi]),
            :
        ]))

        is_winner = r.winner

        value = @sprintf("%.3f", r.AIC)

        if is_winner
            value *= " ★"
        end

        Label(
            fig[i + 2, 1 + 4 * (j - 1) + mi],
            value,
            fontsize = 17,
            font = is_winner ? :bold : :regular,
            color = is_winner ? winner_color : normal_color
        )
    end
end

Label(
    fig[13, 1:13],
    "★ Winner: highest mean AIC weight. 0%: baseline AIC. 40% and 80%: mean AIC across 100 removals.",
    fontsize = 16,
    color = winner_color
)

colgap!(fig.layout, 18)
rowgap!(fig.layout, 16)
resize_to_layout!(fig)

save(
    joinpath(out, "03_all_models_AIC.png"),
    fig;
    px_per_unit = 2
)

# ------------------------------------------------------------
# Markdown report
# ------------------------------------------------------------

md = [
    "# All four models: AIC",
    "",
    "Saved fits only; no models refitted.",
    "",
    "**★ Winner = model with the highest mean AIC weight.**",
    "",
    "At 0%, values are the baseline AIC. At 40% and 80%, values are arithmetic means over the 100 separately fitted removal sequences.",
    "",
    "Winning models are highlighted in bold. Compare models within each community and treatment, not absolute AIC across communities or removal levels.",
    "",
    "The existing best-model selection remains based on mean per-sequence AIC weights. Ranking mean AIC need not give exactly the same ordering as mean AIC weights. No existing selections, figures or KS results have been changed."
]

for f in fractions

    append!(md, [
        "",
        "## $(Int(100f))% removal",
        "",
        "| Community | Power law | Truncated power law | Log-normal | Exponential |",
        "|---|---:|---:|---:|---:|"
    ])

    for (i, ds) in enumerate(datasets)

        values = String[]

        for model in models

            r = only(eachrow(s[
                (s.dataset .== ds) .&
                (s.removal .== f) .&
                (s.model .== model),
                :
            ]))

            value = @sprintf("%.3f", r.AIC)

            if r.winner
                value = "**$(value) ★**"
            end

            push!(values, value)
        end

        push!(
            md,
            "| C$i — $(replace(ds, "_" => " ")) | " *
            join(values, " | ") * " |"
        )
    end
end

push!(
    md,
    "\nSource: fit_diagnostics.csv\n\nSHA256: " *
    bytes2hex(sha256(read(input)))
)

write(
    joinpath(out, "03_all_models_AIC.md"),
    join(md, "\n")
)

println("Saved AIC table, figure and report.")
println("\nWinning distributions:")

show(
    stdout,
    MIME("text/plain"),
    s[s.winner, [
        :dataset,
        :removal,
        :model,
        :AIC,
        :mean_AIC_weight
    ]]
)

println()

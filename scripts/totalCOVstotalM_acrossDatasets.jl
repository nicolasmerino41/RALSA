using CSV, DataFrames, CairoMakie, Statistics, Printf

# Paths
const ROOT = normpath(joinpath(@__DIR__, ".."))
const INPUT = joinpath(ROOT, "outputs/support_prediction_zero/confirmed/pairs.csv")
const DEST = joinpath(ROOT, "outputs/statistical description of local support/total_support_vs_cooccurrence")

mkpath(DEST)

# Dataset order and colours
const DATA = [
    "Quercus", "Nahuel", "Salix_Galpar", "Gottin_HP",
    "Garraf_HP", "Garraf_PP", "Olot", "Montseny"
]

const COLORS = [
    "#332288", "#88CCEE", "#44AA99", "#117733",
    "#999933", "#DDCC77", "#CC6677", "#882255"
]

# Load data
pairs = CSV.read(INPUT, DataFrame)

# Exclude datasets
filter!(row -> !(row.dataset in ["Garraf_PP2", "Gottin_PP"]), pairs)

# Aggregate by dataset
summary = combine(
    groupby(filter(:K => >(0), pairs), :dataset),
    :n => sum => :total_cooccurrences,
    :K => sum => :total_local_support
)

# Preserve dataset order
summary.order = [findfirst(==(d), DATA) for d in summary.dataset]
sort!(summary, :order)

CSV.write(
    joinpath(DEST, "total_support_vs_cooccurrence.csv"),
    select(summary, Not(:order))
)

# Linear regression and Pearson correlation
x = Float64.(summary.total_cooccurrences)
y = Float64.(summary.total_local_support)

pearson_r = cor(x, y)

slope = sum((x .- mean(x)) .* (y .- mean(y))) /
        sum((x .- mean(x)).^2)

intercept = mean(y) - slope * mean(x)

println(@sprintf("Pearson r = %.4f", pearson_r))
println(@sprintf("Linear regression: y = %.4f + %.4f x", intercept, slope))

# Plot
set_theme!(Theme(
    font="Arial",
    fontsize=20,
    Axis=(
        topspinevisible=false,
        rightspinevisible=false,
        xgridvisible=false,
        ygridvisible=false
    )
))

fig = Figure(size=(1100, 800))

ax = Axis(
    fig[1, 1],
    xlabel="Total co-occurrences",
    ylabel="Total local support",
    xlabelsize=22,
    ylabelsize=22
)

# Regression line
xfit = range(minimum(x), maximum(x), length=200)
yfit = intercept .+ slope .* xfit

lines!(
    ax, xfit, yfit,
    color=:black,
    linewidth=2,
    linestyle=:dash
)

# Dataset points
for (i, ds) in enumerate(DATA)

    v = summary[summary.dataset .== ds, :]

    xi = v.total_cooccurrences[1]
    yi = v.total_local_support[1]

    scatter!(
        ax, [xi], [yi],
        color=COLORS[i],
        markersize=18
    )

    text!(
        ax, xi, yi,
        text=replace(ds, "_" => " "),
        align=(:left, :bottom),
        offset=(7, 5),
        fontsize=15,
        color=:black
    )
end

# Pearson annotation
text!(
    ax, 0.05, 0.95,
    text=@sprintf("Pearson r = %.3f", pearson_r),
    space=:relative,
    align=(:left, :top),
    fontsize=22,
    color=:black
)

# Axis limits
xlims!(ax, 0, maximum(x) * 1.20)
ylims!(ax, 0, maximum(y) * 1.15)

save(
    joinpath(DEST, "total_support_vs_cooccurrence.png"),
    fig;
    px_per_unit=2
)

println("Figure saved in: ", DEST)
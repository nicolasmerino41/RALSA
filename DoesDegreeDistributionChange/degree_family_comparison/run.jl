
# degree_family_comparison/run.jl
# Figures 06–09:
# Original / Observed / Shuffled / Shuffled CO
#
# CO-constrained null:
# Randomly assign surviving regional links only to pairs
# with >=1 co-occurrence site remaining after removal.
# Preserve the observed number of surviving links per replicate.
#
# This is a binary-retention shuffle, not a full shuffle of
# the site-by-link interaction matrix.

using CSV, DataFrames, CairoMakie, Random, Statistics

const ROOT = length(ARGS) >= 1 ?
    abspath(ARGS[1]) :
    abspath(joinpath(@__DIR__, "..", ".."))

const OUT = length(ARGS) >= 2 ?
    abspath(ARGS[2]) :
    abspath(joinpath(@__DIR__, "..", "outputs"))

const DATA = [
    "Quercus", "Nahuel", "Salix_Galpar",
    "Gottin_HP", "Gottin_PP",
    "Garraf_HP", "Garraf_PP", "Garraf_PP2",
    "Olot", "Montseny"
]

const FIGURES = [
    ("consumer", 0.4, "06_consumer_40pct_degree_ccdf.png"),
    ("consumer", 0.8, "07_consumer_80pct_degree_ccdf.png"),
    ("resource", 0.4, "08_resource_40pct_degree_ccdf.png"),
    ("resource", 0.8, "09_resource_80pct_degree_ccdf.png")
]

const SCENARIOS = [
    "Original", "Observed", "Shuffled", "Shuffled_CO"
]

const COLORS = Dict(
    "Original" => :black,
    "Observed" => "#0072B2",
    "Shuffled" => "#D55E00",
    "Shuffled_CO" => "#009E73"
)

# ------------------------------------------------------------
# CO-constrained support retention
# ------------------------------------------------------------

function co_histograms()

    checks = CSV.read(
        joinpath(OUT, "checks.csv"), DataFrame
    )

    B = only(unique(Int.(checks.replicates)))

    original_sim = CSV.read(
        joinpath(OUT, "removal_replicates.csv"),
        DataFrame
    )

    observed_counts = Dict{Tuple{String,Int,Int},Int}()

    for r in eachrow(original_sim)
        if r.role == "consumer" && r.scenario == "Observed"
            observed_counts[(
                String(r.dataset),
                round(Int, 100 * r.target_removal),
                Int(r.replicate)
            )] = Int(r.regional_links)
        end
    end

    qc = CSV.read(
        joinpath(
            ROOT, "outputs", "support_prediction_zero",
            "confirmed", "pairs.csv"
        ),
        DataFrame
    )

    qc = qc[qc.K .> 0, :]

    # Reproduce original site-removal sequences.
    site_rng = MersenneTwister(20261028)

    # Independent RNG for the new null.
    co_rng = MersenneTwister(20261029)

    result = NamedTuple[]
    diagnostics = NamedTuple[]

    for ds in DATA

        v = qc[qc.dataset .== ds, :]

        L = nrow(v)
        N = only(unique(Int.(v.N)))

        @assert L > 0

        input_dir = joinpath(
            ROOT, "outputs", "support_prediction_zero",
            "confirmed", "inputs"
        )

        occ = CSV.read(
            joinpath(input_dir, ds * "_occupancy.csv"),
            DataFrame; types=String
        )

        ints = CSV.read(
            joinpath(
                ROOT, "data", "processed",
                ds * "_interactions.csv"
            ),
            DataFrame; types=String
        )

        for column in (:site, :species, :trophic_level)
            column in propertynames(occ) ||
                error("$ds occupancy lacks $column")
        end

        sites = sort(unique(occ.site))

        @assert length(sites) == N

        si = Dict(
            s => i for (i, s) in enumerate(sites)
        )

        links = Dict(
            (r.consumer, r.resource) => i
            for (i, r) in enumerate(eachrow(v))
        )

        # ----------------------------------------------------
        # Reconstruct species occupancy
        # ----------------------------------------------------

        consumers_at = [
            Set{String}() for _ in 1:N
        ]

        resources_at = [
            Set{String}() for _ in 1:N
        ]

        for r in eachrow(occ)

            s = si[r.site]

            if r.trophic_level == "consumer"
                push!(consumers_at[s], r.species)

            elseif r.trophic_level == "resource"
                push!(resources_at[s], r.species)

            else
                error(
                    "Unknown trophic level $(r.trophic_level)"
                )
            end
        end

        # ----------------------------------------------------
        # Reconstruct CO opportunities per site and link
        # ----------------------------------------------------

        co_at = [Int[] for _ in 1:N]

        reconstructed_n = zeros(Int, L)

        for s in 1:N, l in 1:L

            if (v.consumer[l] in consumers_at[s]) &&
               (v.resource[l] in resources_at[s])

                push!(co_at[s], l)
                reconstructed_n[l] += 1
            end
        end

        @assert reconstructed_n == Int.(v.n) (
            "CO counts disagree for $ds"
        )

        # ----------------------------------------------------
        # Reconstruct original interaction support
        # ----------------------------------------------------

        observed_at = [Int[] for _ in 1:N]
        reconstructed_K = zeros(Int, L)

        for r in eachrow(ints)

            s = si[r.site]

            l = get(
                links,
                (r.consumer, r.resource),
                0
            )

            l > 0 || error("Unknown link in $ds")

            @assert r.consumer in consumers_at[s]
            @assert r.resource in resources_at[s]

            push!(observed_at[s], l)
            reconstructed_K[l] += 1
        end

        @assert reconstructed_K == Int.(v.K) (
            "Support counts disagree for $ds"
        )

        @assert all(
            length(unique(a)) == length(a)
            for a in observed_at
        )

        # ----------------------------------------------------
        # Species indices and original degrees
        # ----------------------------------------------------

        name_lists = [
            sort(unique(v.consumer)),
            sort(unique(v.resource))
        ]

        name_to_id = [
            Dict(name => i for (i, name) in
                 enumerate(name_lists[h]))
            for h in 1:2
        ]

        link_ids = [
            [
                name_to_id[h][name]
                for name in (
                    h == 1 ? v.consumer : v.resource
                )
            ]
            for h in 1:2
        ]

        original_degree = [
            [
                count(==(i), link_ids[h])
                for i in eachindex(name_lists[h])
            ]
            for h in 1:2
        ]

        histogram = Dict(
            (j, h) =>
                zeros(Int, maximum(original_degree[h]) + 1)
            for j in 1:2 for h in 1:2
        )

        fractions = [0.4, 0.8]

        mvals = [
            round(Int, f * N) for f in fractions
        ]

        stats = [
            (
                eligible=0.0,
                retained=0.0,
                changed=0.0,
                alternatives=0
            )
            for _ in 1:2
        ]

        # ----------------------------------------------------
        # Replicate site removal and CO-constrained shuffling
        # ----------------------------------------------------

        for b in 1:B

            order = randperm(site_rng, N)

            # Consume the original free-shuffle RNG draw
            # to preserve the exact subsequent site orders.
            randperm(site_rng, L)

            remaining_K = copy(reconstructed_K)
            remaining_n = copy(reconstructed_n)

            previous_m = 0

            for (j, m) in enumerate(mvals)

                # Remove additional sites.
                for pos in (previous_m + 1):m

                    s = order[pos]

                    for l in observed_at[s]
                        remaining_K[l] -= 1
                    end

                    for l in co_at[s]
                        remaining_n[l] -= 1
                    end
                end

                previous_m = m

                @assert all(
                    (remaining_K .>= 0) .&
                    (remaining_K .<= remaining_n)
                )

                kept_observed = remaining_K .> 0

                # Number of regional links surviving
                # in the observed support assignment.
                R = count(kept_observed)

                # Only links with remaining CO are eligible.
                feasible_links = findall(
                    >(0), remaining_n
                )

                @assert length(feasible_links) >= R

                # Validate against original simulation.
                key = (
                    ds,
                    round(Int, 100 * fractions[j]),
                    b
                )

                haskey(observed_counts, key) ||
                    error("Missing replicate $key")

                @assert observed_counts[key] == R (
                    "Site removal differs from original $key"
                )

                # --------------------------------------------
                # CO-constrained null
                #
                # Uniformly choose R distinct links among
                # the links with remaining co-occurrence.
                # --------------------------------------------

                shuffle!(co_rng, feasible_links)

                kept_co = falses(L)

                kept_co[
                    feasible_links[1:R]
                ] .= true

                @assert count(kept_co) == R

                @assert all(
                    !kept_co[l] || remaining_n[l] > 0
                    for l in 1:L
                )

                # --------------------------------------------
                # Consumer and resource degree distributions
                # --------------------------------------------

                for h in 1:2

                    d = zeros(
                        Int,
                        length(name_lists[h])
                    )

                    for l in findall(kept_co)
                        d[link_ids[h][l]] += 1
                    end

                    @assert all(d .<= original_degree[h])
                    @assert sum(d) == R

                    for k in d
                        histogram[(j, h)][k + 1] += 1
                    end
                end

                old = stats[j]

                stats[j] = (
                    eligible=old.eligible +
                        length(feasible_links),

                    retained=old.retained + R,

                    changed=old.changed +
                        count(kept_co .!= kept_observed) / L,

                    alternatives=old.alternatives +
                        (0 < R < length(feasible_links))
                )
            end
        end

        # ----------------------------------------------------
        # Export degree histograms
        # ----------------------------------------------------

        for (j, fraction) in enumerate(fractions)

            for (h, role) in enumerate(
                ("consumer", "resource")
            )

                counts = histogram[(j, h)]

                @assert sum(counts) ==
                    B * length(name_lists[h])

                for k in 0:(length(counts) - 1)

                    push!(
                        result,
                        (
                            dataset=ds,
                            role=role,
                            target_removal=fraction,
                            scenario="Shuffled_CO",
                            degree=k,
                            mean_species_count=
                                counts[k + 1] / B
                        )
                    )
                end
            end

            z = stats[j]

            push!(
                diagnostics,
                (
                    dataset=ds,
                    target_removal=fraction,
                    replicates=B,
                    mean_feasible_links=z.eligible / B,
                    mean_retained_links=z.retained / B,
                    mean_fraction_link_status_changed=
                        z.changed / B,
                    fraction_replicates_with_alternatives=
                        z.alternatives / B
                )
            )
        end

        println(
            "CO-constrained shuffle: $ds completed"
        )

        flush(stdout)
    end

    CSV.write(
        joinpath(
            OUT, "co_constrained_degree_histograms.csv"
        ),
        DataFrame(result)
    )

    CSV.write(
        joinpath(
            OUT, "co_constrained_shuffle_checks.csv"
        ),
        DataFrame(diagnostics)
    )

    return DataFrame(result)
end

# ------------------------------------------------------------
# Original degrees
# ------------------------------------------------------------

function initial_degrees(sp, ds, role)

    v = sp[
        (sp.dataset .== ds) .&
        (sp.role .== role) .&
        (abs.(sp.target_removal .- 0.4) .< 1e-8),
        :
    ]

    isempty(v) && error(
        "Missing initial degrees for $ds, $role"
    )

    @assert length(unique(v.species)) == nrow(v)

    return Int.(v.initial_degree)
end

# ------------------------------------------------------------
# CCDF from integer-degree histograms
# ------------------------------------------------------------

function curve_from_counts(counts)

    active = sum(counts[2:end])

    active > 0 ||
        return Int[], Float64[], Float64[]

    k = collect(1:(length(counts) - 1))

    masses = Float64.(counts[2:end])

    y = reverse(
        cumsum(reverse(masses))
    ) ./ active

    @assert all(diff(y) .<= 1e-9)

    return k, y, masses
end

# ------------------------------------------------------------
# Straight declining segments instead of staircases
# ------------------------------------------------------------

function curve_points(k, y, masses, floor_y)

    visible = (
        (masses .> 1e-12) .&
        (y .>= floor_y - 1e-12)
    )

    return (
        k[visible],
        y[visible],
        visible
    )
end

function plot_line!(ax, x, y, scenario)

    style = scenario in (
        "Shuffled", "Shuffled_CO"
    ) ? :dash : :solid

    if length(x) >= 2

        lines!(
            ax, x, y;
            color=COLORS[scenario],
            linewidth=2.6,
            linestyle=style
        )

    elseif length(x) == 1

        scatter!(
            ax, x, y;
            color=COLORS[scenario],
            markersize=8
        )
    end
end

# ------------------------------------------------------------
# Figures 06–09
# ------------------------------------------------------------

function main()

    sp = CSV.read(
        joinpath(OUT, "species_expected_retention.csv"),
        DataFrame
    )

    old_hist = CSV.read(
        joinpath(OUT, "raw_degree_histograms.csv"),
        DataFrame
    )

    # Calculate additional null.
    co_hist = co_histograms()

    hist = vcat(
        old_hist, co_hist;
        cols=:union
    )

    set_theme!(
        Theme(
            font="Arial",
            fontsize=16,
            Axis=(
                topspinevisible=false,
                rightspinevisible=false,
                xgridvisible=false,
                ygridvisible=false,
                titlesize=17,
                xticklabelsize=12,
                yticklabelsize=12,
                xlabelsize=15,
                ylabelsize=15
            )
        )
    )

    exported = NamedTuple[]

    for (role, fraction, filename) in FIGURES

        fig = Figure(size=(1850, 920))

        for (i, ds) in enumerate(DATA)

            d0 = initial_degrees(sp, ds, role)

            S = length(d0)
            dmax = maximum(d0)

            floor_y = 1 / S

            ax = Axis(
                fig[cld(i, 5), mod1(i, 5)],
                title=replace(ds, "_" => " "),
                xlabel=i > 5 ? "Regional degree k" : "",
                ylabel=mod1(i, 5) == 1 ?
                    "P(D ≥ k | D > 0)" : "",
                xscale=log10,
                yscale=log10
            )

            # --------------------------------------------
            # All four scenarios
            # --------------------------------------------

            for scenario in SCENARIOS

                counts = if scenario == "Original"

                    [
                        count(==(k), d0)
                        for k in 0:dmax
                    ]

                else

                    subset = hist[
                        (hist.dataset .== ds) .&
                        (hist.role .== role) .&
                        (abs.(
                            hist.target_removal .- fraction
                        ) .< 1e-8) .&
                        (hist.scenario .== scenario),
                        :
                    ]

                    @assert nrow(subset) == dmax + 1

                    sort!(subset, :degree)

                    @assert Int.(subset.degree) ==
                        collect(0:dmax)

                    Float64.(subset.mean_species_count)
                end

                @assert isapprox(
                    sum(counts), S;
                    atol=1e-5
                )

                k, y, mass = curve_from_counts(counts)

                isempty(k) && continue

                xx, yy, shown = curve_points(
                    k, y, mass, floor_y
                )

                plot_line!(
                    ax, xx, yy, scenario
                )

                for t in eachindex(k)

                    push!(
                        exported,
                        (
                            dataset=ds,
                            role=role,
                            target_removal=fraction,
                            scenario=scenario,
                            degree=k[t],
                            ccdf=y[t],
                            displayed=shown[t],
                            original_max_degree=dmax,
                            display_floor=floor_y
                        )
                    )
                end
            end

            # --------------------------------------------
            # Identical degree domain within each panel
            # --------------------------------------------

            xlims!(
                ax,
                1,
                max(1.35, 1.10 * dmax)
            )

            ylims!(
                ax,
                max(1e-9, 0.83 * floor_y),
                1.13
            )

            ticks = filter(
                k -> k <= dmax,
                [
                    1, 2, 5, 10, 20, 50,
                    100, 200, 500, 1000
                ]
            )

            ax.xticks = (
                ticks,
                string.(ticks)
            )

            ax.yticks = [
                p for p in (
                    0.001, 0.01, 0.1, 1.0
                )
                if p >= 0.83 * floor_y
            ]
        end

        # ------------------------------------------------
        # Legend
        # ------------------------------------------------

        Legend(
            fig[3, 1:5],
            [
                LineElement(
                    color=COLORS[s],
                    linewidth=3,
                    linestyle=s in (
                        "Shuffled", "Shuffled_CO"
                    ) ? :dash : :solid
                )
                for s in SCENARIOS
            ],
            [
                "Original",
                "Observed support",
                "Shuffled (free)",
                "Shuffled (CO-feasible)"
            ],
            orientation=:horizontal,
            framevisible=false
        )

        Label(
            fig[0, 1:5],
            "$(uppercasefirst(role))s · " *
            "$(round(Int, 100 * fraction))% of sites removed",
            fontsize=22
        )

        Label(
            fig[4, 1:5],
            "Active species only · CO-feasible: ≥1 remaining co-occurrence site · Display floor 1/S₀",
            fontsize=13,
            color=:gray40
        )

        save(
            joinpath(OUT, filename),
            fig;
            px_per_unit=2
        )

        println("Saved $filename")
    end

    CSV.write(
        joinpath(
            OUT, "raw_degree_ccdf_figures_06_09.csv"
        ),
        DataFrame(exported)
    )
end

main()

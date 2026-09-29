include(joinpath(@__DIR__, "..", "src", "data_loading.jl"))
include(joinpath(@__DIR__, "..", "src", "site_removal.jl"))
for script in ["prepare_data", "interaction_support", "cooccurrence_and_interaction_loss",
               "generalists_and_specialists", "degree_distribution"]
    include(joinpath(@__DIR__, script * ".jl"))
end
if !("--figures-only" in ARGS)
    prepare_data()
    interaction_support()
    cooccurrence_and_interaction_loss()
    generalists_and_specialists()
    degree_distribution()
end
if !("--analysis-only" in ARGS)
    include(joinpath(ROOT, "src", "figure_style.jl"))
    for number in 2:4
        include(joinpath(ROOT, "figures", "figure_$(number).jl"))
    end
    figure_2()
    figure_3()
    figure_4()
end
println("RALSA outputs: ", OUT)

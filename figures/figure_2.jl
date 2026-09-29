function figure_2()
        c, sup, pooled = table("retention_and_states"), table("link_support"), table("pooled_retention")
    f2 = Figure(size=(1800, 700), figure_padding=(80, 35, 35, 35))
    a = Axis(f2[1, 1], xlabel= "Sites removed (%)", ylabel= "Retained (%)", xticks= percent_ticks, yticks= percent_ticks)
    for (col, color) in [(:species, SPECIES), (:regional, INTERACTION), (:local_support, SUPPORT)]
        for ds in DATA
            v = c[c.dataset .== ds, :]
            lines!(a, v.removal, v[!, col], color=(color, 0.17), linewidth=1.3)
        end
        lines!(a, pooled.removal, pooled[!, col], color=color, linewidth=3.5)
    end
    xlims!(a, 0, 0.8)
    ylims!(a, 0, 1.02)
    axislegend(a, [LineElement(color=SPECIES, linewidth=3), LineElement(color=INTERACTION, linewidth=3), LineElement(color=SUPPORT, linewidth=3)], ["Species", "Regional interactions", "Local interaction occurrences"]; position=:lb, labelsize=18, patchsize=(25, 12), framevisible=false)
    panel!(f2[1, 1, TopLeft()], "a")
    b = Axis(f2[1, 2], ylabel= "Supporting sites per interaction", yscale= log10, xticks=(1:10, [DISPLAY[s] for s in DATA]), xticklabelrotation= pi / 4, yticks=([1, 2, 5, 10, 20, 50, 100, 200], string.([1, 2, 5, 10, 20, 50, 100, 200])), xticklabelsize=17)
    rng = MersenneTwister(43)
    for (i, ds) in enumerate(DATA)
        v = sup.support[sup.dataset .== ds]
        scatter!(b, i .+ 0.52 .* (rand(rng, length(v)) .- 0.5), v; color=(SUPPORT, 0.23), markersize=4.5)
        boxplot!(b, fill(i, length(v)), v; color=(SUPPORT, 0.25), strokecolor= SUPPORT, mediancolor= :black, whiskercolor= SUPPORT, show_outliers= false, width= 0.52, whiskerwidth= 0.5)
    end
    xlims!(b, 0.4, 10.6)
    ylims!(b, 0.85, maximum(sup.support) * 1.3)
    panel!(f2[1, 2, TopLeft()], "b")

    colgap!(f2.layout,75)
    save_figure(f2,"Figure2_interaction_support")
end

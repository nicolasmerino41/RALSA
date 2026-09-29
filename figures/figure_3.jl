function figure_3()
        pooled, conv = table("pooled_retention"), table("conversion_and_pink_area")
        f2 = Figure(size=(2000,650),figure_padding=(80,45,35,35))
    d = Axis(f2[1, 1], xlabel= "Sites removed (%)", ylabel= "Original interactions (%)", xticks= percent_ticks, yticks= percent_ticks)
    x = pooled.removal
    band!(d, x, zeros(length(x)), pooled.grey, color= CO_LOST)
    band!(d, x, pooled.grey, pooled.grey + pooled.pink, color= CO_ONLY)
    band!(d, x, pooled.grey + pooled.pink, ones(length(x)), color= INTERACTION)
    xlims!(d, 0, 0.8)
    ylims!(d, 0, 1)
    panel!(f2[1, 1, TopLeft()], "a")
    dgrid = GridLayout()
    f2[1, 2] = dgrid
    e1 = Axis(dgrid[1, 1], xlabel= "Empirical local interaction realisation (%)", ylabel= "Mean co-occurrence-only fraction (%)", xticks=([0.1, 0.2, 0.3, 0.4], ["10", "20", "30", "40"]), yticks=([0, 0.05, 0.10, 0.15, 0.20], ["0", "5", "10", "15", "20"]))
    xx1 = range(minimum(conv.p_emp), maximum(conv.p_emp), length=100)
    beta1 = hcat(ones(nrow(conv)), conv.p_emp) \ conv.mean_pink
    lines!(e1, xx1, beta1[1] .+ beta1[2] .* xx1; color=:gray40, linewidth=1.6)
    scatter!(e1, conv.p_emp, conv.mean_pink; color= CO_ONLY, markersize= 13, strokecolor= :white, strokewidth= 1)
    offsets_d1 = Dict("Quercus" => (-8, -15), "Nahuel" => (8, -12), "Salix_Galpar" => (8, 7), "Gottin_HP" => (8, 9), "Gottin_PP" => (8, 8), "Garraf_HP" => (-10, -12), "Garraf_PP" => (12, -28), "Garraf_PP2" => (-8, 8), "Olot" => (-8, 15), "Montseny" => (-8, 9))
    for r in eachrow(conv)
        dx, dy = offsets_d1[r.dataset]
        text!(e1, r.p_emp, r.mean_pink; text= DISPLAY[r.dataset], offset=(dx, dy), align=(dx < 0 ? :right : :left, dy < 0 ? :top : :bottom), fontsize=14)
    end
    xlims!(e1, 0.075, 0.40)
    ylims!(e1, 0, 0.19)
    text!(e1, 0.97, 0.97; space=:relative, text=@sprintf("r = %.2f", cor(conv.p_emp, conv.mean_pink)), align=(:right, :top), fontsize=18)
    e2 = Axis(dgrid[1, 2], xlabel= "Mean co-occurring sites per interaction", ylabel= "")
    xlims!(e2, 2.1, 8.5)
    xx2 = range(minimum(conv.mean_cooccurrence_support), maximum(conv.mean_cooccurrence_support), length=100)
    beta2 = hcat(ones(nrow(conv)), conv.mean_cooccurrence_support) \ conv.mean_pink
    lines!(e2, xx2, beta2[1] .+ beta2[2] .* xx2; color=:gray40, linewidth=1.6)
    scatter!(e2, conv.mean_cooccurrence_support, conv.mean_pink; color= CO_ONLY, markersize= 13, strokecolor= :white, strokewidth= 1)
    offsets_d2 = Dict("Quercus" => (8, -12), "Nahuel" => (8, 8), "Salix_Galpar" => (8, 8), "Gottin_HP" => (8, 8), "Gottin_PP" => (8, 8), "Garraf_HP" => (8, -12), "Garraf_PP" => (8, -12), "Garraf_PP2" => (8, 8), "Olot" => (8, 8), "Montseny" => (8, 8))
    for r in eachrow(conv)
        dx, dy = offsets_d2[r.dataset]
        text!(e2, r.mean_cooccurrence_support, r.mean_pink; text= DISPLAY[r.dataset], offset=(dx, dy), align=(dx < 0 ? :right : :left, dy < 0 ? :top : :bottom), fontsize=14)
    end
    ylims!(e2, 0, 0.19)
    text!(e2, 0.97, 0.97; space=:relative, text=@sprintf("r = %.2f", cor(conv.mean_cooccurrence_support, conv.mean_pink)), align=(:right, :top), fontsize=18)
    hideydecorations!(e2; grid=false)
    panel!(f2[1, 2, TopLeft()], "b")
    Legend(f2[2, 1:2], [PolyElement(color=INTERACTION), PolyElement(color=CO_ONLY), PolyElement(color=CO_LOST)], ["Interaction retained", "Co-occurrence only", "Co-occurrence lost"], orientation=:horizontal)
    colsize!(f2.layout, 1, Relative(.40))
    colsize!(f2.layout, 2, Relative(.60))
    rowgap!(f2.layout, 40)
    colgap!(f2.layout, 75)
    colgap!(dgrid, 50)

    save_figure(f2,"Figure3_cooccurrence_and_interaction_loss")
end

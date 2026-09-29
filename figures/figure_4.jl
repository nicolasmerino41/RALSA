function figure_4()
    ag, rs, dd = table("support_merged"), table("group_curves"), table("degree_example")
    fig = Figure(size=(2050, 850), figure_padding=(75,35,35,35))
    left = fig[1,1] = GridLayout()
    middle = fig[1,2] = GridLayout()
    right = fig[1,3] = GridLayout()
    colsize!(fig.layout, 1, Relative(.25))
    colsize!(fig.layout, 2, Relative(.40))
    colsize!(fig.layout, 3, Relative(.35))
    colgap!(fig.layout, 1, 55)
    colgap!(fig.layout, 2, 28)
    a = Axis(left[1,1], ylabel= "Mean supporting sites per realised partner", xticks=([1,2], GROUPS), limits=(.65, 2.35, 0, nothing))
    for ds in DATA
        d = ag[ag.dataset .== ds, :]
        v = [only(d.support[d.group .== g]) for g in GROUPS]
        lines!(a, [1,2], v, color=(:gray50,.28), linewidth=1.5)
        for i in 1:2
            scatter!(a, [i], [v[i]], color=([SPECIALIST, GENERALIST][i], .5), markersize=10)
        end
    end
    means = [mean(ag.support[ag.group .== g]) for g in GROUPS]
    lines!(a, [1,2], means, color=:gray25, linewidth=2.5)
    for i in 1:2
        v = ag.support[ag.group .== GROUPS[i]]
        q = quantile(v, [.25,.75])
        color = [SPECIALIST, GENERALIST][i]
        rangebars!(a, [i], [q[1]], [q[2]], color=color, linewidth=4, whiskerwidth=16)
        scatter!(a, [i], [means[i]], color=color, markersize=19, strokecolor=:white, strokewidth=1.5)
    end
    ticks = ([0,.2,.4,.6,.8], ["0","20","40","60","80"])
    b1 = Axis(middle[1,1], ylabel="Links lost per original species", xticks=ticks)
    b2 = Axis(middle[2,1], xlabel="Sites removed (%)", ylabel="Species retaining interactions (%)", xticks=ticks, yticks=([0,.25,.5,.75,1], ["0","25","50","75","100"]))
    hidexdecorations!(b1; grid=false)
    linkxaxes!(b1, b2)
    for (group, color) in zip(GROUPS, [SPECIALIST, GENERALIST])
        d = rs[rs.group .== group, :]
        for ds in DATA
            v = d[d.dataset .== ds, :]
            sort!(v, :removal)
            lines!(b1, v.removal, v.links_lost, color=(color,.18), linewidth=1.3)
            lines!(b2, v.removal, v.active_fraction, color=(color,.18), linewidth=1.3)
        end
        p = combine(groupby(d, :removal), :links_lost => mean => :links_lost, :active_fraction => mean => :active_fraction)
        sort!(p, :removal)
        lines!(b1, p.removal, p.links_lost, color=color, linewidth=3.7)
        lines!(b2, p.removal, p.active_fraction, color=color, linewidth=3.7)
    end
    xlims!(b1, 0, .8)
    xlims!(b2, 0, .8)
    ylims!(b1, 0, nothing)
    ylims!(b2, 0, 1.02)
    rowgap!(middle, 25)
    c = Axis(right[1,1], xlabel="Degree", ylabel="Cumulative probability", xscale=log10, yscale=log10, xlabelsize=19, ylabelsize=19, xticklabelsize=17, yticklabelsize=17, xticks=([1,2,5,10,20,50], ["1","2","5","10","20","50"]), yticks=([.01,.1,1], ["0.01","0.1","1"]))
    colors = DEGREE_COLORS
    for (j, f) in enumerate([0.,.4,.8])
        v = dd[dd.removal .== f, :]
        sort!(v, :degree)
        y = v.probability
        keep = [i for i in eachindex(y) if (y[i] > 0 && (i == length(y) || y[i] > y[i+1] + 1e-12))]
        lines!(c, v.degree[keep], y[keep], color=colors[j], linewidth=2.5)
    end
    axislegend(c, [LineElement(color=x, linewidth=3) for x in colors], ["0% removed", "40% removed", "80% removed"]; position=:lb, orientation=:horizontal, framevisible=false, labelsize=16, patchsize=(20,12))
    for (slot, letter) in [(left[1, 1, TopLeft()], "a"), (middle[1, 1, TopLeft()], "b"), (right[1, 1, TopLeft()], "c")]
        Label(slot, letter, font=:bold, fontsize=28, halign=:left, valign=:top, padding=(-42,0,10,0), tellwidth=false, tellheight=false)
    end
    Legend(fig[2,1:3], [LineElement(color=SPECIALIST, linewidth=4), LineElement(color=GENERALIST, linewidth=4)], GROUPS, orientation=:horizontal, framevisible=false, labelsize=21)
    rowsize!(fig.layout, 2, Auto(0.08))
    rowgap!(fig.layout, 15)

    save_figure(fig,"Figure4_generalists_and_specialists")
end

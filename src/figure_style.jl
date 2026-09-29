using CairoMakie
const SUPPORT = "#332288"
const INTERACTION = "#44AA99"
const SPECIES = "#AA4499"
const CO_ONLY = "#CC6677"
const CO_LOST = "#BBBBBB"
const SPECIALIST = "#4477AA"
const GENERALIST = "#CC6677"
const DEGREE_COLORS = ["#88CCEE", "#4477AA", "#332288"]
set_theme!(Theme(font="Arial", fontsize=22, Axis=(backgroundcolor=:white, xgridvisible=false, ygridvisible=false, topspinevisible=false, rightspinevisible=false, spinewidth=1.1, xtickwidth=1.1, ytickwidth=1.1, xticklabelsize=19, yticklabelsize=19, xlabelsize=22, ylabelsize=22), Legend=(framevisible=false, labelsize=19, patchsize=(24, 15))))
percent_ticks = ([0, 0.2, 0.4, 0.6, 0.8, 1.0], string.([0, 20, 40, 60, 80, 100]))
function panel!(slot, letter)
    Label(slot, letter; font=:bold, fontsize=28, halign=:left, valign=:top, tellwidth=false, tellheight=false, padding=(-45, 0, 12, 0))
end
function save_figure(fig, name)
    mkpath(joinpath(OUT, "figures", "main"))
    for ext in ["png"]
        path = joinpath(OUT, "figures", "main", name * "." * ext)
        ext == "png" ? save(path, fig; px_per_unit=2) : save(path, fig)
    end
end

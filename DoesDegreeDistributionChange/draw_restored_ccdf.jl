function draw_restored_ccdf(output)
    root=dirname(dirname(output))
    script=joinpath(dirname(output),"degree_family_comparison","run.jl")
    run(`$(Base.julia_cmd()) --startup-file=no --project=$root $script $root`)
end
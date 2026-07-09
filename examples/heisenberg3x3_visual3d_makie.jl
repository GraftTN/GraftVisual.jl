# CairoMakie version of the 3D rotating visual check. Scene construction comes
# from src/GraftVisual.jl; this script only does the Makie drawing. Run with:
#
#   julia --project=. examples/heisenberg3x3_visual3d_makie.jl
#
# Outputs:
#   - heisenberg3x3_visual3d_makie.png / .svg / .pdf   static snapshot
#   - heisenberg3x3_visual3d_makie.gif                 camera orbit animation

using GraftVisual
using CairoMakie
const GV = GraftVisual

root = normpath(joinpath(@__DIR__, ".."))
artifact_dir = joinpath(root, "artifacts")
mkpath(artifact_dir)
out(ext) = joinpath(artifact_dir, "heisenberg3x3_visual3d_makie.$ext")

function draw_scene!(ax, sc)
    segs = Dict(kind => Pair{Point3f,Point3f}[] for kind in (:bond, :nn, :grid, :leg))
    arrows_from = Point3f[]
    arrows_dir = Vec3f[]
    for (a, b, kind) in sc.edges
        pa, pb = Point3f(sc.pos[a]...), Point3f(sc.pos[b]...)
        if kind == :tree
            q = GV.ctrlpoint(sc.pos[a], sc.pos[b])
            pts = [Point3f(GV.bez(sc.pos[a], q, sc.pos[b], t / 24)...) for t in 0:24]
            lines!(ax, pts; color=:black, linewidth=1.8)
        elseif kind == :arrow
            push!(arrows_from, pa)
            push!(arrows_dir, Vec3f(pb - pa))
        else
            push!(segs[kind], pa => pb)
        end
    end
    isempty(segs[:grid]) ||
        linesegments!(ax, segs[:grid]; color=GV.GRID_GRAY, linewidth=1.0, linestyle=:dash)
    isempty(segs[:bond]) ||
        linesegments!(ax, segs[:bond]; color=:black, linewidth=1.4, linestyle=:dash)
    isempty(segs[:nn]) ||
        linesegments!(ax, segs[:nn]; color=GV.NN_RED, linewidth=2.6)
    isempty(segs[:leg]) ||
        linesegments!(ax, segs[:leg]; color=:black, linewidth=1.6)
    isempty(arrows_from) ||
        arrows3d!(ax, arrows_from, arrows_dir; color=:black,
            shaftradius=0.012, tipradius=0.04, tiplength=0.12)

    ks = sort([k for k in keys(sc.pos) if sc.style[k].r > 0]; by=string)
    pts = [Point3f(sc.pos[k]...) for k in ks]
    halos = [sc.style[k].halo for k in ks]
    any(halos) && scatter!(ax, pts[halos];
        color=:white, markersize=[2 * sc.style[k].r + 6 for k in ks if sc.style[k].halo])
    scatter!(ax, pts;
        color=[sc.style[k].fill for k in ks],
        markersize=[2 * sc.style[k].r for k in ks],
        strokecolor=[sc.style[k].strokewidth > 0 ? sc.style[k].stroke : "#00000000" for k in ks],
        strokewidth=[max(sc.style[k].strokewidth, 0.01) for k in ks])

    labeled = sort([k for k in keys(sc.pos) if !isempty(sc.style[k].label)]; by=string)
    isempty(labeled) || text!(ax, [Point3f(sc.pos[k]...) for k in labeled];
        text=[sc.style[k].label for k in labeled],
        offset=(0, -16), align=(:center, :top), fontsize=12,
        color=[sc.style[k].labelcolor for k in labeled], font=:bold_italic)
    sublabeled = sort([k for k in keys(sc.pos) if !isempty(sc.style[k].sublabel)]; by=string)
    isempty(sublabeled) || text!(ax, [Point3f(sc.pos[k]...) for k in sublabeled];
        text=[sc.style[k].sublabel for k in sublabeled],
        offset=(0, -30), align=(:center, :top), fontsize=9, color=:black)
end

const AZ0, EL0 = 1.65π, 0.42

fig = Figure(size=(1500, 560), backgroundcolor=:white)
Label(fig[0, 1:3], "3x3 Heisenberg visual check - 3D view"; fontsize=24, font=:bold)
Label(fig[1, 1:3],
    "model graph: SAFIRE-style lattice, red = nearest-neighbor bonds; " *
    "tree panels: PRX-style tree tensors, blue = low level, green = high, dark green = root";
    fontsize=13, color="#475569", tellwidth=false)

panels = [("model graph", model_scene(3, 3)),
          ("balanced binary leaf tree", binary_tree_scene(3, 3; leaftensors=true)),
          ("comb / T3NS tree", comb_scene(3, 3))]

axes = map(enumerate(panels)) do (i, (title, sc))
    ax = Axis3(fig[2, i]; aspect=:data, azimuth=AZ0, elevation=EL0,
        viewmode=:fitzoom, protrusions=0)
    hidedecorations!(ax)
    hidespines!(ax)
    draw_scene!(ax, sc)
    Label(fig[3, i], title; fontsize=17, font=:bold, tellwidth=false)
    ax
end
rowgap!(fig.layout, 8)

save(out("png"), fig; px_per_unit=2)
save(out("svg"), fig)
save(out("pdf"), fig)

# :fitzoom refits the zoom every frame (visible breathing while orbiting), so
# switch to the rotation-stable :fit before recording
for ax in axes
    ax.viewmode = :fit
end
record(fig, out("gif"), range(AZ0, AZ0 + 2π; length=73)[1:72]; framerate=18) do az
    for ax in axes
        ax.azimuth[] = az
    end
end

foreach(ext -> println(out(ext)), ("png", "svg", "pdf", "gif"))

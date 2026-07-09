# 3D rotating-view visual check for the 3x3 open-boundary Heisenberg model.
# All rendering machinery lives in src/GraftVisual.jl. Run with:
#
#   julia --project=. examples/heisenberg3x3_visual3d.jl
#
# Outputs:
#   - heisenberg3x3_visual3d.svg / .png   static snapshot at a fixed camera angle
#   - heisenberg3x3_visual3d.gif          camera orbiting the scene (needs magick)
#   - heisenberg3x3_visual3d.html         self-contained interactive viewer
#                                         (auto-rotate toggle, drag to orbit)

using GraftVisual

root = normpath(joinpath(@__DIR__, ".."))
artifact_dir = joinpath(root, "artifacts")
mkpath(artifact_dir)
out(ext) = joinpath(artifact_dir, "heisenberg3x3_visual3d.$ext")

panels = [
    Panel("model graph", model_scene(3, 3)),
    Panel("balanced binary leaf tree", [
        "full tree: explicit leaf tensors" => binary_tree_scene(3, 3; leaftensors=true),
        "compact: tree attaches to sites" => binary_tree_scene(3, 3; leaftensors=false, dz=0.62),
    ]),
    Panel("comb / T3NS tree", comb_scene(3, 3)),
]

title = "3x3 Heisenberg visual check - 3D view"
subtitle = "model graph: SAFIRE-style lattice, red = nearest-neighbor bonds; " *
           "tree panels: PRX-style tree tensors over the lattice, blue = low level, green = high, dark green = root"

write_snapshot_svg(out("svg"), panels; title, subtitle)
Sys.which("magick") !== nothing && run(`magick $(out("svg")) $(out("png"))`)
write_orbit_gif(out("gif"), panels; title, subtitle)
write_viewer_html(out("html"), panels; title, subtitle=subtitle * ". Drag to orbit.")

for ext in ("svg", "png", "gif", "html")
    isfile(out(ext)) && println(out(ext))
end

# GraftVisual: visual checks for tree tensor networks over 2D lattice models.
#
# Style: PRX Quantum "Hybrid Tree Tensor Networks" figures for the 3D tree
# scenes (lattice flat in z = 0, colored tensor circles, drooping black
# curves) and SAFIRE's afqmctools `plot_lattice` for the model-graph scene
# (light gray sites with black edge + white halo, royal-blue bold-italic
# labels with coordinates, red nearest-neighbor bonds, dashed grid, lattice
# vectors).
#
# Rendering targets:
#   - write_viewer_html : self-contained interactive canvas viewer
#                         (auto-rotate toggle, drag to orbit, scene variants)
#   - write_snapshot_svg: static SVG at a fixed camera angle
#   - write_orbit_gif   : rotating GIF (needs ImageMagick's `magick`)
#
# The SVG path is restricted to primitives that ImageMagick's built-in SVG
# rasterizer (MSVG, no librsvg) renders faithfully: black <line>, colored
# <circle>/<polygon> fills, <text>. Curves are sampled into short segments;
# colored or light strokes are drawn as filled quads.

module GraftVisual

export NodeStyle, Scene, Panel,
    square_lattice_scene, model_scene, binary_tree_scene, comb_scene,
    write_viewer_html, write_snapshot_svg, write_orbit_gif

const ROYAL_BLUE = "#001452"   # SAFIRE site-label color
const NN_RED = "#d62728"       # SAFIRE nearest-neighbor bond color
const GRID_GRAY = "#c9ced6"

Base.@kwdef struct NodeStyle
    r::Float64 = 11.0
    fill::String = "#4f81e8"
    stroke::String = "white"
    strokewidth::Float64 = 1.6
    halo::Bool = false            # SAFIRE-style white halo ring
    label::String = ""            # bold italic, below the node
    sublabel::String = ""         # monospace, below the label
    labelcolor::String = ROYAL_BLUE
end

struct Scene
    pos::Dict{Symbol,NTuple{3,Float64}}
    style::Dict{Symbol,NodeStyle}
    # (a, b, kind); kinds: :bond (dashed black), :nn (solid red), :tree
    # (curved black), :leg (straight black), :grid (dashed light gray),
    # :arrow (black, arrowhead at b)
    edges::Vector{Tuple{Symbol,Symbol,Symbol}}
end
Scene() = Scene(Dict{Symbol,NTuple{3,Float64}}(), Dict{Symbol,NodeStyle}(),
    Tuple{Symbol,Symbol,Symbol}[])

struct Panel
    title::String
    variants::Vector{Pair{String,Scene}}
end
Panel(title::AbstractString, scene::Scene) = Panel(title, ["default" => scene])

# --- scene builders -----------------------------------------------------------

site(i) = Symbol(:s, i)
sitexy(i, Lx) = ((i - 1) % Lx, (i - 1) ÷ Lx)
function sitepos(i, Lx, Ly)
    x, y = sitexy(i, Lx)
    (x - (Lx - 1) / 2, y - (Ly - 1) / 2, 0.0)
end

function lattice_bonds(Lx, Ly)
    bonds = Tuple{Int,Int}[]
    for y in 0:Ly-1, x in 1:Lx
        i = x + Lx * y
        x < Lx && push!(bonds, (i, i + 1))
        y < Ly - 1 && push!(bonds, (i, i + Lx))
    end
    return bonds
end

"Open-boundary square lattice as small black dots joined by dashed bonds."
function square_lattice_scene(Lx=3, Ly=3)
    sc = Scene()
    for i in 1:Lx*Ly
        sc.pos[site(i)] = sitepos(i, Lx, Ly)
        sc.style[site(i)] = NodeStyle(r=5.0, fill="#111111", strokewidth=0.0)
    end
    for (i, j) in lattice_bonds(Lx, Ly)
        push!(sc.edges, (site(i), site(j), :bond))
    end
    return sc
end

"""
SAFIRE-style model graph: light gray sites with black edge and white halo,
royal-blue bold-italic labels with monospace coordinates below, red
nearest-neighbor bonds, a dashed background grid, and lattice vectors a1/a2.
"""
function model_scene(Lx=3, Ly=3)
    sc = Scene()
    ghost = Ref(0)
    function node!(p, st)
        ghost[] += 1
        g = Symbol(:g, ghost[])
        sc.pos[g] = p
        sc.style[g] = st
        return g
    end
    hidden = NodeStyle(r=0.0)

    # dashed grid through every lattice row/column, extended past the sites
    xs = [sitepos(i, Lx, Ly)[1] for i in 1:Lx]
    ys = [sitepos(1 + Lx * k, Lx, Ly)[2] for k in 0:Ly-1]
    lo = (minimum(xs) - 0.6, minimum(ys) - 0.6)
    hi = (maximum(xs) + 0.6, maximum(ys) + 0.6)
    for x in xs
        push!(sc.edges, (node!((x, lo[2], 0.0), hidden), node!((x, hi[2], 0.0), hidden), :grid))
    end
    for y in ys
        push!(sc.edges, (node!((lo[1], y, 0.0), hidden), node!((hi[1], y, 0.0), hidden), :grid))
    end

    for (i, j) in lattice_bonds(Lx, Ly)
        push!(sc.edges, (site(i), site(j), :nn))
    end
    for i in 1:Lx*Ly
        x, y = sitexy(i, Lx)
        sc.pos[site(i)] = sitepos(i, Lx, Ly)
        sc.style[site(i)] = NodeStyle(r=10.0, fill="#d9d9d9", stroke="#000000",
            strokewidth=1.2, halo=true, label="s$i", sublabel="($x,$y)")
    end

    # lattice vectors from an offset corner, SAFIRE-style
    o = (lo[1] - 0.35, lo[2] - 0.35, 0.0)
    push!(sc.edges, (node!(o, hidden), node!((o[1] + 1, o[2], 0.0), hidden), :arrow))
    push!(sc.edges, (node!(o, hidden), node!((o[1], o[2] + 1, 0.0), hidden), :arrow))
    node!((o[1] + 0.55, o[2] - 0.22, 0.0), NodeStyle(r=0.0, label="a1", labelcolor="#000000"))
    node!((o[1] - 0.28, o[2] + 0.55, 0.0), NodeStyle(r=0.0, label="a2", labelcolor="#000000"))
    return sc
end

levelcolor(h, hmax) =
    h == hmax ? "#2e7d32" : h == 1 ? "#7da3f2" : h == 2 ? "#4f81e8" : "#5cb85c"

"""
Balanced binary leaf tree over the lattice (recursive halving of the
row-major site order). Internal tensors float above the centroid of their
descendant sites at a height proportional to subtree height. With
`leaftensors=true` every site also gets an explicit leaf tensor above its
lattice dot (PRX style); otherwise the tree attaches directly to the dots.
"""
function binary_tree_scene(Lx=3, Ly=3; dz=0.52, leaftensors=true)
    sc = square_lattice_scene(Lx, Ly)
    heights = Dict{Symbol,Int}()
    counter = Ref(0)
    # NOTE: locals here must not share names with locals of the enclosing
    # function — assignments in a nested function to enclosing-scope names
    # are captured as one shared (boxed) variable and corrupt the recursion.
    function build(idxs)
        if length(idxs) == 1
            i = only(idxs)
            leaftensors || return site(i), 0
            lf = Symbol(:l, i)
            sp = sitepos(i, Lx, Ly)
            sc.pos[lf] = (sp[1], sp[2], dz)
            heights[lf] = 1
            push!(sc.edges, (lf, site(i), :leg))
            return lf, 1
        end
        counter[] += 1
        n = Symbol(:b, counter[])
        cut = length(idxs) ÷ 2
        child1 = build(idxs[1:cut])
        child2 = build(idxs[cut+1:end])
        h = max(child1[2], child2[2]) + 1
        ps = [sitepos(i, Lx, Ly) for i in idxs]
        sc.pos[n] = (sum(p -> p[1], ps) / length(ps), sum(p -> p[2], ps) / length(ps), h * dz)
        heights[n] = h
        push!(sc.edges, (n, child1[1], :tree), (n, child2[1], :tree))
        return n, h
    end
    hroot = build(1:Lx*Ly)[2]
    for (k, hk) in heights
        sc.style[k] = NodeStyle(r=hk == hroot ? 12.5 : 11.0, fill=levelcolor(hk, hroot))
    end
    return sc
end

"""
Comb / T3NS tree: one virtual tensor directly above each site, chained in
row-major site order; the chain jumps at the row boundaries show up as the
long arcing edges.
"""
function comb_scene(Lx=3, Ly=3; dz=0.95)
    sc = square_lattice_scene(Lx, Ly)
    for i in 1:Lx*Ly
        p = sitepos(i, Lx, Ly)
        t = Symbol(:t, i)
        sc.pos[t] = (p[1], p[2], dz)
        sc.style[t] = NodeStyle(r=11.0, fill="#4f81e8")
        push!(sc.edges, (t, site(i), :leg))
        i > 1 && push!(sc.edges, (Symbol(:t, i - 1), t, :tree))
    end
    return sc
end

# --- camera -------------------------------------------------------------------

"Orthographic projection: yaw θ around z, elevation α. Returns (sx, sy, depth)."
function project(p, θ, α)
    x1 = p[1] * cos(θ) - p[2] * sin(θ)
    y1 = p[1] * sin(θ) + p[2] * cos(θ)
    return (x1, -(p[3] * cos(α) + y1 * sin(α)), y1 * cos(α) - p[3] * sin(α))
end

# Quadratic bezier control point giving the PRX-style droop: branches leave
# the parent heading outward and sag into the child from above. Same-height
# edges (the comb chain) arc gently upward instead.
function ctrlpoint(pa, pb)
    if abs(pa[3] - pb[3]) < 1e-9
        d = hypot(pb[1] - pa[1], pb[2] - pa[2])
        return ((pa[1] + pb[1]) / 2, (pa[2] + pb[2]) / 2, pa[3] + 0.12 + 0.05d)
    end
    hi, lo = pa[3] >= pb[3] ? (pa, pb) : (pb, pa)
    return (0.25hi[1] + 0.75lo[1], 0.25hi[2] + 0.75lo[2], lo[3] + 0.35(hi[3] - lo[3]))
end

bez(a, q, b, t) = ntuple(k -> (1 - t)^2 * a[k] + 2(1 - t) * t * q[k] + t^2 * b[k], 3)

# --- SVG rendering (MSVG-safe primitives only) ----------------------------------

r1(x) = round(x; digits=1)

"Solid segment as a filled quad (survives rasterizers that drop line colors)."
function quadsvg(a, b, w, color; opacity="")
    len = max(hypot(b[1] - a[1], b[2] - a[2]), 1e-9)
    nx, ny = -(b[2] - a[2]) / len * w / 2, (b[1] - a[1]) / len * w / 2
    op = isempty(opacity) ? "" : " opacity=\"$opacity\""
    return "<polygon points=\"$(r1(a[1] + nx)),$(r1(a[2] + ny)) $(r1(b[1] + nx)),$(r1(b[2] + ny)) $(r1(b[1] - nx)),$(r1(b[2] - ny)) $(r1(a[1] - nx)),$(r1(a[2] - ny))\" fill=\"$color\"$op/>"
end

"Dashed segment as a chain of filled quads."
function dashquadsvg(a, b, w, color; dash=7.0, gap=5.0)
    len = max(hypot(b[1] - a[1], b[2] - a[2]), 1e-9)
    ux, uy = (b[1] - a[1]) / len, (b[2] - a[2]) / len
    out = String[]
    t = 0.0
    while t < len
        t2 = min(t + dash, len)
        push!(out, quadsvg((a[1] + ux * t, a[2] + uy * t), (a[1] + ux * t2, a[2] + uy * t2), w, color))
        t = t2 + gap
    end
    return join(out, "")
end

function halotextsvg(x, y, text, size, color; bold=false, italic=false, mono=false)
    fam = mono ? " font-family=\"Courier,monospace\"" : ""
    w = bold ? " font-weight=\"700\"" : ""
    it = italic ? " font-style=\"italic\"" : ""
    attrs = "text-anchor=\"middle\" font-size=\"$size\"$fam$w$it"
    out = String[]
    for (dx, dy) in ((-0.9, 0), (0.9, 0), (0, -0.9), (0, 0.9))
        push!(out, "<text x=\"$(r1(x + dx))\" y=\"$(r1(y + dy))\" $attrs fill=\"white\">$text</text>")
    end
    push!(out, "<text x=\"$(r1(x))\" y=\"$(r1(y))\" $attrs fill=\"$color\">$text</text>")
    return join(out, "")
end

"Largest horizontal radius of a scene — bounds the projected extent under rotation."
scene_radius(sc::Scene) = maximum(hypot(p[1], p[2]) for p in values(sc.pos); init=1.0)

function render_panel_svg(io, title, sc::Scene, θ, α; cx, cy, s=110, panel_w=500, panel_h=480)
    s = min(s, (panel_w / 2 - 12) / scene_radius(sc))   # fit at every yaw angle
    items = Tuple{Float64,String}[]   # (depth, svg fragment); far first
    P = Dict(k => project(v, θ, α) for (k, v) in sc.pos)
    scr(q) = (cx + s * q[1], cy + s * q[2])
    for (a, b, kind) in sc.edges
        pa, pb = scr(P[a]), scr(P[b])
        d = (P[a][3] + P[b][3]) / 2
        if kind == :tree
            q = ctrlpoint(sc.pos[a], sc.pos[b])
            n = 18
            prj = [project(bez(sc.pos[a], q, sc.pos[b], t / n), θ, α) for t in 0:n]
            for k in 1:n
                (x1, y1), (x2, y2) = scr(prj[k]), scr(prj[k+1])
                push!(items, ((prj[k][3] + prj[k+1][3]) / 2,
                    "<line x1=\"$(r1(x1))\" y1=\"$(r1(y1))\" x2=\"$(r1(x2))\" y2=\"$(r1(y2))\" stroke=\"#111111\" stroke-width=\"1.8\" stroke-linecap=\"round\"/>"))
            end
        elseif kind == :bond
            push!(items, (d + 0.05,
                "<line x1=\"$(r1(pa[1]))\" y1=\"$(r1(pa[2]))\" x2=\"$(r1(pb[1]))\" y2=\"$(r1(pb[2]))\" stroke=\"#111111\" stroke-width=\"1.5\" stroke-dasharray=\"7 5\"/>"))
        elseif kind == :leg
            push!(items, (d + 0.05,
                "<line x1=\"$(r1(pa[1]))\" y1=\"$(r1(pa[2]))\" x2=\"$(r1(pb[1]))\" y2=\"$(r1(pb[2]))\" stroke=\"#111111\" stroke-width=\"1.6\"/>"))
        elseif kind == :nn
            push!(items, (d + 0.05, quadsvg(pa, pb, 2.6, NN_RED)))
        elseif kind == :grid
            push!(items, (d + 0.5, dashquadsvg(pa, pb, 1.0, GRID_GRAY; dash=5.0, gap=4.0)))
        elseif kind == :arrow
            len = max(hypot(pb[1] - pa[1], pb[2] - pa[2]), 1e-9)
            ux, uy = (pb[1] - pa[1]) / len, (pb[2] - pa[2]) / len
            base = (pb[1] - 10ux, pb[2] - 10uy)
            head = "<polygon points=\"$(r1(pb[1])),$(r1(pb[2])) $(r1(base[1] - 4uy)),$(r1(base[2] + 4ux)) $(r1(base[1] + 4uy)),$(r1(base[2] - 4ux))\" fill=\"#111111\"/>"
            push!(items, (d, quadsvg(pa, base, 2.4, "#111111") * head))
        end
    end
    for (k, q) in P
        st = sc.style[k]
        x, y = scr(q)
        frag = String[]
        if st.r > 0
            st.halo && push!(frag, "<circle cx=\"$(r1(x))\" cy=\"$(r1(y))\" r=\"$(st.r + 3)\" fill=\"white\"/>")
            sw = st.strokewidth > 0 ? " stroke=\"$(st.stroke)\" stroke-width=\"$(st.strokewidth)\"" : ""
            push!(frag, "<circle cx=\"$(r1(x))\" cy=\"$(r1(y))\" r=\"$(st.r)\" fill=\"$(st.fill)\"$sw/>")
        end
        isempty(st.label) ||
            push!(frag, halotextsvg(x, y + st.r + 13, st.label, 11, st.labelcolor; bold=true, italic=true))
        isempty(st.sublabel) ||
            push!(frag, halotextsvg(x, y + st.r + 24, st.sublabel, 8.5, "#111111"; mono=true))
        isempty(frag) || push!(items, (q[3] - 0.02, join(frag, "")))
    end
    sort!(items; by=first, rev=true)
    foreach(x -> println(io, x[2]), items)
    println(io, "<text x=\"$cx\" y=\"$(panel_h - 22)\" text-anchor=\"middle\" font-size=\"17\" font-weight=\"700\">$title</text>")
end

function frame_svg(io, panels::Vector{Panel}, θ; α=deg2rad(28), title="", subtitle="",
    panel_w=500, panel_h=480, s=110)
    W = panel_w * length(panels)
    println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$W\" height=\"$panel_h\" viewBox=\"0 0 $W $panel_h\">")
    println(io, "<style>text{font-family:Helvetica,Arial,sans-serif;fill:#0f172a}</style>")
    println(io, "<rect width=\"$W\" height=\"$panel_h\" fill=\"white\"/>")
    isempty(title) ||
        println(io, "<text x=\"$(W ÷ 2)\" y=\"34\" text-anchor=\"middle\" font-size=\"22\" font-weight=\"700\">$title</text>")
    isempty(subtitle) ||
        println(io, "<text x=\"$(W ÷ 2)\" y=\"56\" text-anchor=\"middle\" font-size=\"12.5\" fill=\"#475569\">$subtitle</text>")
    for (i, panel) in enumerate(panels)
        render_panel_svg(io, panel.title, panel.variants[1][2], θ, α;
            cx=panel_w * i - panel_w ÷ 2, cy=panel_h - 130, s, panel_w, panel_h)
    end
    println(io, "</svg>")
end

"Static snapshot of the first variant of each panel."
function write_snapshot_svg(path, panels; θ=deg2rad(-32), kwargs...)
    open(io -> frame_svg(io, panels, θ; kwargs...), path, "w")
    return path
end

"Camera-orbit GIF via ImageMagick; returns the path or nothing if unavailable."
function write_orbit_gif(path, panels; θ0=deg2rad(-32), nframes=48, width=900, kwargs...)
    Sys.which("magick") === nothing && return nothing
    mktempdir() do dir
        frames = map(0:nframes-1) do f
            fp = joinpath(dir, "frame_$(lpad(f, 3, '0')).svg")
            open(io -> frame_svg(io, panels, θ0 + 2π * f / nframes; kwargs...), fp, "w")
            fp
        end
        run(`magick -delay 7 -loop 0 $frames -resize $(width)x $path`)
    end
    return path
end

# --- interactive HTML viewer ----------------------------------------------------

jsstr(x) = "\"" * replace(string(x), "\\" => "\\\\", "\"" => "\\\"") * "\""

function js_scene(sc::Scene)
    ids = sort(collect(keys(sc.pos)); by=string)
    idx = Dict(k => i - 1 for (i, k) in enumerate(ids))
    ns = join((begin
            st = sc.style[k]
            "{p:[$(join(sc.pos[k], ","))],r:$(st.r),f:$(jsstr(st.fill)),st:$(jsstr(st.stroke)),sw:$(st.strokewidth),halo:$(st.halo),lb:$(jsstr(st.label)),sub:$(jsstr(st.sublabel)),lc:$(jsstr(st.labelcolor))}"
        end for k in ids), ",")
    es = join(("[$(idx[a]),$(idx[b]),$(jsstr(kind))]" for (a, b, kind) in sc.edges), ",")
    return "{nodes:[$ns],edges:[$es]}"
end

js_panel(p::Panel) =
    "{title:$(jsstr(p.title)),variants:[" *
    join(("{label:$(jsstr(lb)),scene:$(js_scene(sc))}" for (lb, sc) in p.variants), ",") * "]}"

const HTML_TEMPLATE = raw"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8"/>
<title>__TITLE__</title>
<style>
  body{font-family:Helvetica,Arial,sans-serif;color:#0f172a;background:#fff;margin:24px;text-align:center}
  h1{font-size:22px;margin:0 0 4px}
  p.sub{font-size:13px;color:#475569;margin:2px 0 10px}
  .row{display:flex;justify-content:center;gap:24px;flex-wrap:wrap}
  figure{margin:0}
  canvas{border:1px solid #e2e8f0;border-radius:8px;cursor:grab;touch-action:none}
  figcaption{font-size:16px;font-weight:700;margin-top:8px}
  .opt{font-size:12.5px;color:#475569;display:block;margin:0 0 12px;cursor:pointer}
  select{font-size:12px;color:#334155;margin-top:5px}
</style>
</head>
<body>
<h1>__TITLE__</h1>
<p class="sub">__SUBTITLE__</p>
<label class="opt"><input type="checkbox" id="spin" checked/> auto-rotate</label>
<div class="row" id="rowdiv"></div>
<script>
var PANELS = [__PANELS__];
var SPIN = true;

function project(p, th, al) {
  var x1 = p[0] * Math.cos(th) - p[1] * Math.sin(th);
  var y1 = p[0] * Math.sin(th) + p[1] * Math.cos(th);
  return { x: x1,
           y: -(p[2] * Math.cos(al) + y1 * Math.sin(al)),
           d: y1 * Math.cos(al) - p[2] * Math.sin(al) };
}

function ctrlpoint(pa, pb) {
  if (Math.abs(pa[2] - pb[2]) < 1e-9) {
    var d = Math.hypot(pb[0] - pa[0], pb[1] - pa[1]);
    return [(pa[0] + pb[0]) / 2, (pa[1] + pb[1]) / 2, pa[2] + 0.12 + 0.05 * d];
  }
  var hi = pa[2] >= pb[2] ? pa : pb, lo = pa[2] >= pb[2] ? pb : pa;
  return [0.25 * hi[0] + 0.75 * lo[0], 0.25 * hi[1] + 0.75 * lo[1], lo[2] + 0.35 * (hi[2] - lo[2])];
}

function bez(a, q, b, t) {
  var u = 1 - t, out = [];
  for (var k = 0; k < 3; k++) out.push(u * u * a[k] + 2 * u * t * q[k] + t * t * b[k]);
  return out;
}

var EDGE_STYLE = {
  bond:  { color: "#111", width: 1.5, dash: [7, 5] },
  leg:   { color: "#111", width: 1.6, dash: [] },
  nn:    { color: "#d62728", width: 2.5, dash: [] },
  grid:  { color: "#c9ced6", width: 1.0, dash: [4, 4] },
  arrow: { color: "#111", width: 2.2, dash: [] }
};

function setup(canvas, view) {
  var ctx = canvas.getContext("2d");
  var dpr = window.devicePixelRatio || 1;
  var cw = canvas.width, ch = canvas.height;
  canvas.width = cw * dpr; canvas.height = ch * dpr;
  canvas.style.width = cw + "px"; canvas.style.height = ch + "px";
  ctx.scale(dpr, dpr);
  var th = -0.56, al = 0.49, s = 100, cx = cw / 2, cy = ch - 120;
  var dragging = false, px = 0, py = 0;

  canvas.addEventListener("pointerdown", function (e) {
    dragging = true; px = e.clientX; py = e.clientY;
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener("pointermove", function (e) {
    if (!dragging) return;
    th += (e.clientX - px) * 0.01;
    al = Math.min(1.4, Math.max(0.08, al + (e.clientY - py) * 0.006));
    px = e.clientX; py = e.clientY;
  });
  canvas.addEventListener("pointerup", function () { dragging = false; });

  function haloText(text, x, y, font, color) {
    ctx.font = font;
    ctx.textAlign = "center"; ctx.textBaseline = "middle";
    ctx.strokeStyle = "white"; ctx.lineWidth = 4; ctx.lineJoin = "round";
    ctx.strokeText(text, x, y);
    ctx.fillStyle = color;
    ctx.fillText(text, x, y);
  }

  function draw() {
    var scene = view.scene;
    var R = 1;
    scene.nodes.forEach(function (n) { R = Math.max(R, Math.hypot(n.p[0], n.p[1])); });
    s = Math.min(100, (cw / 2 - 12) / R);   // fit at every yaw angle
    if (SPIN && !dragging) th += 0.006;
    ctx.clearRect(0, 0, cw, ch);
    var items = [];
    var P = scene.nodes.map(function (n) { return project(n.p, th, al); });
    scene.edges.forEach(function (e) {
      var pa = scene.nodes[e[0]].p, pb = scene.nodes[e[1]].p, kind = e[2];
      if (kind === "tree") {
        var q = ctrlpoint(pa, pb), n = 18, prj = [];
        for (var t = 0; t <= n; t++) prj.push(project(bez(pa, q, pb, t / n), th, al));
        for (var k = 0; k < n; k++) (function (a, b) {
          items.push({ d: (a.d + b.d) / 2, fn: function () {
            ctx.beginPath();
            ctx.moveTo(cx + s * a.x, cy + s * a.y);
            ctx.lineTo(cx + s * b.x, cy + s * b.y);
            ctx.strokeStyle = "#111"; ctx.lineWidth = 1.8;
            ctx.lineCap = "round"; ctx.setLineDash([]);
            ctx.stroke();
          }});
        })(prj[k], prj[k + 1]);
      } else {
        var a = P[e[0]], b = P[e[1]];
        var stl = EDGE_STYLE[kind] || EDGE_STYLE.leg;
        var bias = kind === "grid" ? 0.5 : 0.05;
        items.push({ d: (a.d + b.d) / 2 + bias, fn: function () {
          ctx.beginPath();
          ctx.moveTo(cx + s * a.x, cy + s * a.y);
          ctx.lineTo(cx + s * b.x, cy + s * b.y);
          ctx.strokeStyle = stl.color; ctx.lineWidth = stl.width;
          ctx.setLineDash(stl.dash);
          ctx.stroke(); ctx.setLineDash([]);
          if (kind === "arrow") {
            var dx = (b.x - a.x), dy = (b.y - a.y);
            var L = Math.hypot(dx, dy) || 1e-9;
            var ux = dx / L, uy = dy / L;
            var tipx = cx + s * b.x, tipy = cy + s * b.y;
            ctx.beginPath();
            ctx.moveTo(tipx, tipy);
            ctx.lineTo(tipx - 10 * ux - 4 * uy, tipy - 10 * uy + 4 * ux);
            ctx.lineTo(tipx - 10 * ux + 4 * uy, tipy - 10 * uy - 4 * ux);
            ctx.closePath();
            ctx.fillStyle = stl.color; ctx.fill();
          }
        }});
      }
    });
    scene.nodes.forEach(function (n, i) {
      var q = P[i];
      items.push({ d: q.d - 0.02, fn: function () {
        var x = cx + s * q.x, y = cy + s * q.y;
        if (n.r > 0) {
          ctx.beginPath();
          ctx.arc(x, y, n.r, 0, 2 * Math.PI);
          if (n.halo) { ctx.strokeStyle = "white"; ctx.lineWidth = 6; ctx.stroke(); }
          ctx.fillStyle = n.f; ctx.fill();
          if (n.sw > 0) { ctx.strokeStyle = n.st; ctx.lineWidth = n.sw; ctx.stroke(); }
        }
        if (n.lb) haloText(n.lb, x, y + n.r + 13, "italic 700 11px Helvetica,Arial,sans-serif", n.lc);
        if (n.sub) haloText(n.sub, x, y + n.r + 25, "10px Courier,monospace", "#111");
      }});
    });
    items.sort(function (a, b) { return b.d - a.d; });
    items.forEach(function (it) { it.fn(); });
    requestAnimationFrame(draw);
  }
  requestAnimationFrame(draw);
}

var row = document.getElementById("rowdiv");
PANELS.forEach(function (p, i) {
  var fig = document.createElement("figure");
  var canvas = document.createElement("canvas");
  canvas.width = 440; canvas.height = 430;
  fig.appendChild(canvas);
  var cap = document.createElement("figcaption");
  cap.textContent = p.title;
  fig.appendChild(cap);
  var view = { scene: p.variants[0].scene };
  if (p.variants.length > 1) {
    var sel = document.createElement("select");
    p.variants.forEach(function (v, j) {
      var opt = document.createElement("option");
      opt.value = j; opt.textContent = v.label;
      sel.appendChild(opt);
    });
    sel.addEventListener("change", function (e) {
      view.scene = p.variants[+e.target.value].scene;
    });
    fig.appendChild(sel);
  }
  row.appendChild(fig);
  setup(canvas, view);
});

document.getElementById("spin").addEventListener("change", function (e) {
  SPIN = e.target.checked;
});
</script>
</body>
</html>
"""

"""
    write_viewer_html(path, panels; title, subtitle)

Write a self-contained interactive viewer: one auto-rotating, drag-to-orbit
canvas per panel, a global auto-rotate toggle, and a variant selector under
any panel with more than one scene variant.
"""
function write_viewer_html(path, panels; title="GraftVisual", subtitle="")
    open(path, "w") do io
        print(io, replace(HTML_TEMPLATE,
            "__TITLE__" => title,
            "__SUBTITLE__" => subtitle,
            "__PANELS__" => join(map(js_panel, panels), ",")))
    end
    return path
end

end # module

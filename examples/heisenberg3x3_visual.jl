root = normpath(joinpath(@__DIR__, ".."))
artifact_dir = joinpath(root, "artifacts")
mkpath(artifact_dir)

svg_path = joinpath(artifact_dir, "heisenberg3x3_visual.svg")
png_path = joinpath(artifact_dir, "heisenberg3x3_visual.png")
pdf_path = joinpath(artifact_dir, "heisenberg3x3_visual.pdf")

sites = Symbol.("s" .* string.(1:9))
bonds = Tuple{Int,Int}[]
for y in 0:2, x in 1:3
    i = x + 3y
    x < 3 && push!(bonds, (i, i + 1))
    y < 2 && push!(bonds, (i, i + 3))
end

esc(x) = replace(string(x), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
pt(p) = "$(round(p[1]; digits=1)),$(round(p[2]; digits=1))"

function line(io, a, b; color="#334155", width=2, dash=false)
    d = dash ? " stroke-dasharray=\"7 5\"" : ""
    println(io, "<line x1=\"$(a[1])\" y1=\"$(a[2])\" x2=\"$(b[1])\" y2=\"$(b[2])\" stroke=\"$color\" stroke-width=\"$width\" stroke-linecap=\"round\"$d/>")
end

function curve(io, a, b; color="#dc2626")
    midy = min(a[2], b[2]) - 28 - abs(a[1] - b[1]) / 14
    println(io, "<path d=\"M $(pt(a)) C $(a[1]),$midy $(b[1]),$midy $(pt(b))\" fill=\"none\" stroke=\"$color\" stroke-width=\"3\" stroke-dasharray=\"9 5\" stroke-linecap=\"round\" opacity=\"0.82\"/>")
end

function node(io, p, label; physical=false)
    if physical
        println(io, "<circle cx=\"$(p[1])\" cy=\"$(p[2])\" r=\"17\" fill=\"#e0f2fe\" stroke=\"#0369a1\" stroke-width=\"2\"/>")
    else
        println(io, "<rect x=\"$(p[1]-22)\" y=\"$(p[2]-12)\" width=\"44\" height=\"24\" rx=\"4\" fill=\"#f8fafc\" stroke=\"#475569\" stroke-width=\"1.5\"/>")
    end
    println(io, "<text x=\"$(p[1])\" y=\"$(p[2]+4)\" text-anchor=\"middle\" font-size=\"11\" font-weight=\"600\">$(esc(label))</text>")
end

function binary_layout(labels; x0=330, y0=130, dx=42, dy=78)
    pos = Dict{Symbol,Tuple{Float64,Float64}}()
    edges = Tuple{Symbol,Symbol}[]
    counter = Ref(0)
    function build(parent, labels, depth)
        if length(labels) == 1
            leaf = only(labels)
            pos[leaf] = (x0 + (parse(Int, string(leaf)[2:end]) - 1) * dx, y0 + 4dy)
            push!(edges, (parent, leaf))
            return pos[leaf][1]
        end
        counter[] += 1
        n = Symbol(:b, counter[])
        push!(edges, (parent, n))
        cut = length(labels) ÷ 2
        xs = (build(n, labels[1:cut], depth + 1), build(n, labels[cut+1:end], depth + 1))
        pos[n] = (sum(xs) / 2, y0 + depth * dy)
        return pos[n][1]
    end
    pos[:broot] = (x0 + 4dx, y0)
    build(:broot, labels[1:4], 1)
    build(:broot, labels[5:9], 1)
    return pos, edges
end

function comb_layout(labels; x0=735, y0=170, dx=45)
    pos = Dict{Symbol,Tuple{Float64,Float64}}()
    edges = Tuple{Symbol,Symbol}[]
    for i in eachindex(labels)
        sp = Symbol(:t, i)
        pos[sp] = (x0 + (i - 1) * dx, y0)
        pos[labels[i]] = (x0 + (i - 1) * dx, y0 + 210)
        push!(edges, (sp, labels[i]))
        i > 1 && push!(edges, (Symbol(:t, i - 1), sp))
    end
    return pos, edges
end

function draw_tree(io, title, pos, edges, xmid)
    println(io, "<text x=\"$xmid\" y=\"95\" text-anchor=\"middle\" font-size=\"18\" font-weight=\"700\">$(esc(title))</text>")
    for e in edges
        line(io, pos[e[1]], pos[e[2]])
    end
    for (i, j) in bonds
        curve(io, pos[sites[i]], pos[sites[j]])
    end
    for (k, p) in sort(collect(pos); by=x -> string(x[1]))
        node(io, p, k; physical=k in sites)
    end
end

open(svg_path, "w") do io
    println(io, raw"""<svg xmlns="http://www.w3.org/2000/svg" width="1180" height="560" viewBox="0 0 1180 560">""")
    println(io, raw"""<style>text{font-family:Helvetica,Arial,sans-serif;fill:#0f172a}</style>""")
    println(io, raw"""<rect width="1180" height="560" fill="white"/>""")
    println(io, raw"""<text x="590" y="42" text-anchor="middle" font-size="24" font-weight="700">3x3 Heisenberg open-boundary visual check</text>""")
    println(io, raw"""<text x="590" y="66" text-anchor="middle" font-size="13" fill="#475569">solid lines: tree/lattice edges; dashed red lines: Heisenberg nearest-neighbor bonds overlaid on tree leaves</text>""")

    println(io, raw"""<text x="125" y="112" text-anchor="middle" font-size="18" font-weight="700">model graph</text>""")
    lattice = Dict(sites[x + 1 + 3y] => (55 + 70x, 165 + 70y) for y in 0:2 for x in 0:2)
    for (i, j) in bonds
        line(io, lattice[sites[i]], lattice[sites[j]]; color="#2563eb", width=2.4)
    end
    for s in sites
        node(io, lattice[s], s; physical=true)
    end

    bpos, bedges = binary_layout(sites)
    draw_tree(io, "balanced binary leaf tree", bpos, bedges, 500)
    cpos, cedges = comb_layout(sites)
    draw_tree(io, "comb / T3NS tree", cpos, cedges, 915)
    println(io, "</svg>")
end

if Sys.which("magick") !== nothing
    run(`magick $svg_path $png_path`)
    run(`magick $svg_path $pdf_path`)
end

println(svg_path)
isfile(png_path) && println(png_path)
isfile(pdf_path) && println(pdf_path)

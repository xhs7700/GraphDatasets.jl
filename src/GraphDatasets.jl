module GraphDatasets

using GeneralGraphs
using Base.Filesystem
using Downloads
using Tar, TranscodingStreams, CodecBzip2, CodecZlib
using ProgressBars

export loadUndiKONECT, loadUndiSNAP
export loadPseudofractal, loadKoch, loadCayleyTree, loadHanoiExt, loadApollo, loadPseudoExt, load3CayleyTree, loadCorona

function callback(desc::AbstractString)
    pbar = Ref{ProgressBar}()
    init = Ref(false)
    now = Ref(0)
    return (total, new_now) -> begin
        if total != 0
            if init[]
                if new_now > now[]
                    update(pbar[], new_now - now[])
                    now[] = new_now
                end
            else
                pbar[] = ProgressBar(total=total)
                set_description(pbar[], desc)
                init[] = true
            end
        end
    end
end

function loadUndiKONECT(internal_name::AbstractString, name::AbstractString)
    url = "http://konect.cc/files/download.tsv.$internal_name.tar.bz2"
    bz_io = IOBuffer()
    Downloads.download(url, bz_io; progress=callback(internal_name))
    tar_bytes = transcode(Bzip2Decompressor, take!(bz_io))
    dir_path = mktempdir()
    Tar.extract(IOBuffer(tar_bytes), dir_path)
    file_dir = joinpath(dir_path, internal_name)
    file_path = nothing
    for file_name in readdir(file_dir)
        if startswith(file_name, "out.")
            file_path = joinpath(file_dir, file_name)
        end
    end
    if isnothing(file_path)
        error("cannot find graph file in $file_dir.")
    end
    g = GeneralGraph{Int}(() -> 1, name, file_path)
    rm(dir_path; recursive=true)
    return g
end

loadUndiKONECT(internal_name::AbstractString) = loadUndiKONECT(internal_name, internal_name)

function loadUndiSNAP(url::AbstractString, name::AbstractString)
    gzip_io = IOBuffer()
    Downloads.download(url, gzip_io; progress=callback(name))
    txt_bytes = transcode(GzipDecompressor, take!(gzip_io))
    return GeneralGraph{Int}(() -> 1, name, IOBuffer(txt_bytes))
end

loadUndiSNAP(url::AbstractString) = loadUndiSNAP(url, "UndiSNAP")

@doc raw"""
    loadPseudofractal(g)

Generate the GeneralGraph object of the pseudofractal network ``F_g``.

The number of vertices and edges of ``F_g`` is ``\left(3^{g+1} + 3\right) / 2`` and ``3^{g+1}``.

The Kemeny constant of ``F_g`` is

```math
\frac{5}{2}\times 3^g - \frac{5}{3}\times 2^g + \frac{1}{2}.
```
"""
loadPseudofractal(g::Int) = _loadPseudoExt(1, g, "Pseudofractal_$g")

@doc raw"""
    loadPseudoExt(m,g)

Generate the GeneralGraph object of the extended pseudofractal network ``F_{m,g}``.

The number of vertices and edges of ``F_{m,g}`` is ``3\times\left( (2m+1)^g + 1\right) / 2`` and ``3\times(2m+1)^g``.

The Kemeny constant of ``F_{m,g}`` is

```math
1 + \frac{1}{30m\left(2m+1\right)}\times\left( 8 + 25m + 18m^2 + (135m+90m^2)(1+2m)^g - (28m^2+120m+8)\left(\frac{2+4m}{2+m}\right)^g \right).
```
"""
loadPseudoExt(m::Int, g::Int) = _loadPseudoExt(m, g, "PseudoExt_$(m)_$g")

function _loadPseudoExt(m::Int, g::Int, name::AbstractString)
    edges = Tuple{Int,Int}[(1, 2), (1, 3), (2, 3)]
    N = 3
    for _ in 1:g
        new_edges = Tuple{Int,Int}[]
        sizehint!(new_edges, 2 * m * length(edges))
        for (u, v) in edges
            for _ in 1:m
                N += 1
                push!(new_edges, (u, N), (v, N))
            end
        end
        append!(edges, new_edges)
    end
    return GeneralGraph(name, Set(1:N), Dict(edge => one(Int) for edge in edges))
end

@doc raw"""
    loadCorona(q,g)

Generate the GeneralGraph object of the edge corona project of ``(q+2)``-Clique ``G_q(g)``.

The number of vertices and edges of ``G_q(g)`` is ``\frac{2}{q+3}\left(\frac{(q+1)(q+2)}{2}\right)^{g+1} + \frac{2q+4}{q+3}`` and ``\left(\frac{(q+1)(q+2)}{2}\right)^{g+1}``.

The Kemeny constant of ``G_q(g)`` is

```math
\left(\frac{(q+1)^2}{q+2} - \frac{3q+3}{2}\right)\times (q+1)^g + \frac{(q+1)(3q+7)}{2q+6}\times\left(\frac{(q+1)(q+2)}{2}\right)^g + \frac{q+1}{q+3}
```
"""
function loadCorona(q::Int, g::Int)
    q += 2
    edges = Tuple{Int,Int}[]
    N = q
    for u in 1:q, v in u+1:q
        push!(edges, (u, v))
    end
    for _ in 1:g
        new_edges = Tuple{Int,Int}[]
        sizehint!(new_edges, (2 * q + q * (q - 1) ÷ 2) * length(edges))
        for (u, v) in edges
            for x in 1:q
                push!(new_edges, (u, N + x), (v, N + x))
            end
            for x in 1:q, y in x+1:q
                push!(new_edges, (N + x, N + y))
            end
            N += q
        end
        append!(edges, new_edges)
    end
    return GeneralGraph(name, Set(1:N), Dict(edge => one(Int) for edge in edges))
end

@doc raw"""
    loadKoch(g)

Generate the GeneralGraph object of the Koch network ``M_g``.

The number of vertices and edges of ``M_g`` is ``2\times 4^g +1`` and ``3\times 4^g``.

The Kemeny constant of ``M_g`` is

```math
\left(1 + 2g\right)\times 4^g + \frac{1}{3}.
```
"""
function loadKoch(g::Int)
    triangles = Tuple{Int,Int,Int}[(1, 2, 3)]
    N = 3
    for _ in 1:g
        new_triangles = Tuple{Int,Int,Int}[]
        sizehint!(new_triangles, 3 * length(triangles))
        for triangle in triangles, u in triangle
            push!(new_triangles, (u, N + 1, N + 2))
            N += 2
        end
        append!(triangles, new_triangles)
    end
    edges = Dict{Tuple{Int,Int},Int}()
    for (x, y, z) in triangles
        edges[(x, y)] = 1
        edges[(x, z)] = 1
        edges[(y, z)] = 1
    end
    return GeneralGraph("Koch_$g", Set(1:N), edges)
end

@doc raw"""
    load3CayleyTree(g)

Generate the GeneralGraph object of the 3-Cayley Tree ``C_{3,g}``.

The number of vertices and edges of ``C_{3,g}`` is ``3\times 2^g - 2`` and ``3\times 2^g - 3``.

The Kemeny constant of ``C_{3,g}`` is

```math
\frac{3g\times 4^{g+1} - 13\times 2^{2g+1} + 35\times 2^g - 9}{2\left(2^g - 1\right)}.
```
"""
load3CayleyTree(g::Int) = loadCayleyTree(3, g)

function loadCayleyTree(b::Int, g::Int)
    edges = Dict{Tuple{Int,Int},Int}()
    leafs = Int[]
    for leaf in 2:b+1
        edges[(1, leaf)] = 1
        push!(leafs, leaf)
    end
    N = b + 1
    for _ in 2:g
        new_leafs = Int[]
        for leaf in leafs
            for new_leaf in N+1:N+b-1
                edges[(leaf, new_leaf)] = 1
                push!(new_leafs, new_leaf)
            end
            N += b - 1
        end
        leafs = new_leafs
    end
    return GeneralGraph("CayleyTree_$(b)_$(g)", Set(1:N), edges)
end

@doc raw"""
    loadHanoiExt(g)

Generate the GeneralGraph object of the extended Hanoi graph ``\tilde{H}_g``.

The number of vertices and edges of ``\tilde{H}_g`` is ``4\times 3^{g-1}`` and ``2\times 3^g``.

The Kemeny constant of ``\tilde{H}_g`` is

```math
\frac{32\times 5^g\times 3^{g-1}-64\times 3^{2g-2}-2\times 3^g}{10\left( 3^g+3^{g-1}-1 \right)}.
```
"""
function loadHanoiExt(g::Int)
    edges = Tuple{Int,Int}[(1, 2), (1, 3), (2, 3)]
    inc = 1
    for i in 2:g-1
        inc *= 3
        new_edges = Tuple{Int,Int}[]
        for (u, v) in edges
            push!(new_edges, (u + inc, v + inc), (u + 2 * inc, v + 2 * inc))
        end
        append!(edges, new_edges)
        for x in 0:2
            y = (x + 1) % 3
            u = parse(Int, string(x) * repeat(string(y), i - 1), base=3) + 1
            v = parse(Int, string(y) * repeat(string(x), i - 1), base=3) + 1
            push!(edges, (u, v))
        end
    end
    inc *= 3
    new_edges = Tuple{Int,Int}[]
    for (u, v) in edges
        push!(new_edges, (u + inc, v + inc), (u + 2 * inc, v + 2 * inc), (u + 3 * inc, v + 3 * inc))
    end
    append!(edges, new_edges)
    for x in 0:2
        y = (x + 1) % 3
        u = parse(Int, string(x) * repeat(string(y), g - 1), base=3) + 1
        v = parse(Int, string(y) * repeat(string(x), g - 1), base=3) + 1
        push!(edges, (u, v))
        u = parse(Int, repeat(string(x), g), base=3) + 1
        v = parse(Int, "10" * repeat(string(x), g - 1), base=3) + 1
        push!(edges, (u, v))
    end
    N = 4 * 3^(g - 1)
    return GeneralGraph("HanoiExt_$g", Set(1:N), Dict(edge => one(Int) for edge in edges))
end

@doc raw"""
    loadApollo(g)

Generate the GeneralGraph object of the Apollonian network ``A_g``.

The number of vertices and edges of ``A_g`` is ``2\times 3^g+2`` and ``6\times 3^g``.

The Kemeny constant of ``A_g`` is

```math
1+\frac{1}{12}\left( 32\times 3^g-16\left(\frac{9}{5}\right)^g +11 \right).
```
"""
function loadApollo(g::Int)
    triangles = NTuple{3,Int}[(1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4)]
    N = 4
    for _ in 1:g
        new_triangles = NTuple{3,Int}[]
        sizehint!(new_triangles, 3 * length(triangles))
        for (x, y, z) in triangles
            N += 1
            push!(new_triangles, (x, y, N), (x, z, N), (y, z, N))
        end
        append!(triangles, new_triangles)
    end
    edges = Dict{Tuple{Int,Int},Int}()
    for (x, y, z) in triangles
        edges[(x, y)] = 1
        edges[(x, z)] = 1
        edges[(y, z)] = 1
    end
    return GeneralGraph("Apollo_$g", Set(1:N), edges)
end

end # module GraphDatasets

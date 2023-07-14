module GraphDatasets

using GeneralGraphs
using Base.Filesystem
using Downloads
using Tar, TranscodingStreams, CodecBzip2, CodecZlib
using ProgressBars

export loadUndiKONECT, loadUndiSNAP
export loadPseudofractal, loadKoch, loadCayleyTree, loadHanoiExt, loadApollo, loadPseudoExt, load3CayleyTree

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

loadPseudofractal(g::Int) = _loadPseudoExt(1, g, "Pseudofractal_$g")

loadPseudoExt(m::Int, g::Int) = _loadPseudoExt(m, g, "PseudoExt_$(m)_$g")

function _loadPseudoExt(m::Int, g::Int, name::AbstractString)
    edges = Tuple{Int,Int}[(1, 2), (1, 3), (2, 3)]
    N = 3
    for _ in 1:g
        new_edges = Tuple{Int,Int}()
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

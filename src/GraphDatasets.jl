module GraphDatasets

using GeneralGraphs
using Base.Filesystem
using Downloads
using Tar, TranscodingStreams, CodecBzip2, CodecZlib
using ProgressBars

export loadUndiKONECT, loadUndiSNAP

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
    file_path = joinpath(dir_path, internal_name, "out.$internal_name")
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

end # module GraphDatasets

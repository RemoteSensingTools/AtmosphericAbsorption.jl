#=
One-time dev script: build data/tips2021.arrow — the latest HITRAN total internal
partition sums (TIPS-2021, Gamache et al. 2021) — from HITRAN's per-isotopologue
partition-sum files at https://hitran.org/data/Q/q{global_id}.txt (two columns: T[K],
Q(T)). The (mol, iso) → global_id map is data/hitran_global_ids.csv (HITRAN molparam
order). Output schema matches tips2017.arrow: a long (mol, iso, T, Q) table that
src/PartitionFunctions/tips2021.jl groups and splines on demand.

Run: julia --project=. data/_build_tips2021.jl   (needs network access to hitran.org)
=#
using Arrow, Downloads

const HERE = @__DIR__
const QURL = "https://hitran.org/data/Q"

# (mol, iso) → global_id
rows = Tuple{Int,Int,Int}[]
for ln in eachline(joinpath(HERE, "hitran_global_ids.csv"))
    startswith(ln, "mol") && continue
    a = split(strip(ln), ',')
    isempty(a[1]) && continue
    push!(rows, (parse(Int, a[1]), parse(Int, a[2]), parse(Int, a[3])))
end

# HITRAN ships Q on a 1 K grid to 5000 K; Q(T) is smooth, so we keep a physics-weighted
# subset (fine through the atmospheric range, coarser at high T — 296 K stays a node) to
# hold the table near tips2017's size without measurable spline error.
const TGRID = sort(unique(Float64.([1:5:300; 300:20:1000; 1000:100:5000])))

tmol = Int32[]; tiso = Int32[]; tT = Float64[]; tQ = Float64[]
missing_gid = Tuple{Int,Int,Int}[]
for (M, I, g) in rows
    buf = IOBuffer()
    ok = try
        Downloads.request("$QURL/q$(g).txt"; output = buf, throw = true); true
    catch
        false
    end
    ok || (push!(missing_gid, (M, I, g)); continue)
    Q = Dict{Float64,Float64}()
    for ln in eachline(seekstart(buf))
        s = split(strip(ln))
        length(s) < 2 && continue
        Tk = tryparse(Float64, s[1]); Qk = tryparse(Float64, s[2])
        (Tk === nothing || Qk === nothing) && continue
        Q[Tk] = Qk
    end
    isempty(Q) && (push!(missing_gid, (M, I, g)); continue)   # empty/HTML response
    for Tk in TGRID
        haskey(Q, Tk) || continue
        push!(tmol, M); push!(tiso, I); push!(tT, Tk); push!(tQ, Q[Tk])
    end
end

Arrow.write(joinpath(HERE, "tips2021.arrow"), (; mol = tmol, iso = tiso, T = tT, Q = tQ))
npairs = length(Set(zip(tmol, tiso)))
println("tips2021.arrow: ", length(tmol), " rows, ", npairs, " (mol,iso) pairs")
isempty(missing_gid) || println("no Q-file for: ", missing_gid)

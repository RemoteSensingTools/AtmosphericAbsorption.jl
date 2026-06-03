#=
One-time dev script: convert the static NetCDF lookup tables (iso_info.nc,
TIPS_2017.nc) shipped by vSmartMOM into lightweight Arrow files, so the package
needs only Arrow (not NetCDF/HDF5) at runtime. Kept for provenance.

Run once from the package environment (which still has NetCDF available):
    julia --project=. data/_build_arrow.jl
=#
using NetCDF, Arrow

here = @__DIR__

# --- iso_info: matrices [mol, iso] → long table -----------------------------------
iso = joinpath(here, "iso_info.nc")
weight   = ncread(iso, "mol_weight")     # [g/mol], fill -1
abund    = ncread(iso, "abundance")
gid      = ncread(iso, "global_id")
molname  = ncread(iso, "mol_name")       # [mol, k] strings

nmol, niso = size(weight)
mol = Int32[]; isn = Int32[]; mm = Float64[]; ab = Float64[]; g = Int32[]; nm = String[]
for m in 1:nmol, i in 1:niso
    gid[m, i] == -1 && continue
    push!(mol, m); push!(isn, i)
    push!(mm, weight[m, i]); push!(ab, abund[m, i]); push!(g, Int32(gid[m, i]))
    push!(nm, strip(String(molname[m, 1])))
end
Arrow.write(joinpath(here, "iso_info.arrow"),
            (; mol, iso = isn, molar_mass = mm, abundance = ab, global_id = g, mol_name = nm))
println("iso_info.arrow: ", length(mol), " (mol,iso) rows")

# --- TIPS-2017: (mol, iso, T) cube → long table (valid entries only) ---------------
tips = joinpath(here, "TIPS_2017.nc")
TT = ncread(tips, "TIPS_2017_T")
QQ = ncread(tips, "TIPS_2017_Q")
M, I, _ = size(TT)
tmol = Int32[]; tiso = Int32[]; tT = Float64[]; tQ = Float64[]
for m in 1:M, i in 1:I
    e = findfirst(==(-1), @view TT[m, i, :])
    n = e === nothing ? size(TT, 3) : e - 1
    n == 0 && continue
    for k in 1:n
        push!(tmol, m); push!(tiso, i); push!(tT, TT[m, i, k]); push!(tQ, QQ[m, i, k])
    end
end
Arrow.write(joinpath(here, "tips2017.arrow"),
            (; mol = tmol, iso = tiso, T = tT, Q = tQ))
println("tips2017.arrow: ", length(tmol), " (mol,iso,T) rows")

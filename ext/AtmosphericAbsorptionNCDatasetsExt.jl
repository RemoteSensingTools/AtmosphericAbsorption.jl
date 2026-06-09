module AtmosphericAbsorptionNCDatasetsExt

# Reader for AER ABSCO `.hdf` files (HDF5, read natively by NCDatasets). Maps the file's
# Gas_NN_Absorption(P, T, broadener, ν) cube — which NCDatasets returns transposed to Julia order
# (ν, broadener, T, P) — onto an `AbscoLUT`, converting Pressure Pa→hPa. Kept in an extension so the
# core package stays NetCDF-free; `using NCDatasets` enables `read_absco`.

using NCDatasets
using AtmosphericAbsorption: AbscoLUT
import AtmosphericAbsorption: read_absco
using AtmosphericAbsorption.Architectures: default_architecture, AbstractArchitecture

function read_absco(path::AbstractString; scale::Real = 1.0,
                    FT::Type{<:AbstractFloat} = Float32,
                    architecture::AbstractArchitecture = default_architecture())
    NCDataset(path) do ds
        # Locate the Gas_NN_Absorption cube by name (robust to how Gas_Index is stored — string,
        # zero-padded or numeric) and take the molecule id from that name.
        names = collect(keys(ds))
        i = findfirst(k -> occursin(r"^Gas_0*\d+_Absorption$", k), names)
        i === nothing && error("read_absco: no `Gas_NN_Absorption` variable found in $path")
        csname = names[i]
        molnum = parse(Int, match(r"Gas_(\d+)_Absorption", csname).captures[1])
        raw = Array(ds[csname])                                      # Julia order (ν, broadener, T, P)
        σ = eltype(raw) === FT ? raw : FT.(raw)
        scale == 1 || (σ = σ .* FT(scale))
        T   = Matrix{FT}(Array(ds["Temperature"]))                  # (T_idx, P)
        p   = Vector{FT}(Array(ds["Pressure"])) ./ 100              # Pa → hPa
        ν   = Vector{FT}(Array(ds["Wavenumber"]))                   # cm⁻¹
        vmr = haskey(ds, "Broadener_01_VMR") ? Vector{FT}(Array(ds["Broadener_01_VMR"])) : FT[0]
        # NCDatasets reverses the file's C dim order, so the cube must arrive as (ν, broadener, T, P);
        # fail loudly rather than silently mis-look-up if a file ever differs.
        size(σ) == (length(ν), length(vmr), size(T, 1), length(p)) ||
            error("read_absco: $csname has shape $(size(σ)); expected (ν, broadener, T, P) = " *
                  "$((length(ν), length(vmr), size(T, 1), length(p)))")
        return AbscoLUT(molnum, -1, ν, p, T, vmr, σ; architecture)
    end
end

end # module

"""
Build the six Float32 native-grid ABSCO LUTs used by the three OCO-2 bands.

Run with the package environment (NCDatasets must be available):

    ABSCO_ROOT=/path/to/v5.2_final \
    ABSCO_LUT_OUTPUT=/path/to/output \
    julia --project=. benchmark/build_oco2_absco_luts.jl

The products are CPU-portable `.absco` files. Load directly on a GPU with
`load_absco_lut(path; architecture=GPU())`.
"""

using AtmosphericAbsorption
using NCDatasets
using Printf
using AtmosphericAbsorption.Architectures: CPU

const ABSCO_ROOT = get(ENV, "ABSCO_ROOT",
    "/kiwi-data/Data/Spectroscopy/ABSCO_CS_Database/v5.2_final")
const OUTPUT_ROOT = get(ENV, "ABSCO_LUT_OUTPUT", joinpath(ABSCO_ROOT, "oco2_luts_f32"))

const PRODUCTS = (
    ("o2_v52.hdf",  :o2,   "o2_o2_v52_f32.absco"),
    ("co2_v52.hdf", :wco2, "co2_wco2_v52_f32.absco"),
    ("co2_v52.hdf", :sco2, "co2_sco2_v52_f32.absco"),
    ("h2o_v52.hdf", :o2,   "h2o_o2_v52_f32.absco"),
    ("h2o_v52.hdf", :wco2, "h2o_wco2_v52_f32.absco"),
    ("h2o_v52.hdf", :sco2, "h2o_sco2_v52_f32.absco"),
)

mkpath(OUTPUT_ROOT)
for (source_name, band, output_name) in PRODUCTS
    source = joinpath(ABSCO_ROOT, source_name)
    output = joinpath(OUTPUT_ROOT, output_name)
    elapsed = @elapsed begin
        lut = read_oco2_absco(source, band; FT=Float32, architecture=CPU(),
                              broadener_vmr=:all)
        save_absco_lut(output, lut)
    end
    @printf("%-27s shape=%-24s %7.1f MiB  %6.2f s\n",
            output_name, string(size(lut.σ)), filesize(output) / 2.0^20, elapsed)
end

println("Wrote Float32 OCO-2 ABSCO products to $OUTPUT_ROOT")

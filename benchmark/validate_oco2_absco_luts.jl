"""
Validate and benchmark the six OCO-2 Float32 LUTs against ABSCO v5.2 HDF nodes.

    ABSCO_ROOT=/path/to/v5.2_final ABSCO_LUT_OUTPUT=/path/to/oco2_luts_f32 \
      julia --project=<env-with-NCDatasets-and-optionally-CUDA> \
      benchmark/validate_oco2_absco_luts.jl
"""

using AtmosphericAbsorption
using LinearAlgebra
using NCDatasets
using Printf
using Statistics
using AtmosphericAbsorption.Architectures: CPU, GPU

const ABSCO_ROOT = get(ENV, "ABSCO_ROOT",
    "/kiwi-data/Data/Spectroscopy/ABSCO_CS_Database/v5.2_final")
const LUT_ROOT = get(ENV, "ABSCO_LUT_OUTPUT", joinpath(ABSCO_ROOT, "oco2_luts_f32"))
const PRODUCTS = (
    ("o2_v52.hdf",  "o2_o2_v52_f32.absco"),
    ("co2_v52.hdf", "co2_wco2_v52_f32.absco"),
    ("co2_v52.hdf", "co2_sco2_v52_f32.absco"),
    ("h2o_v52.hdf", "h2o_o2_v52_f32.absco"),
    ("h2o_v52.hdf", "h2o_wco2_v52_f32.absco"),
    ("h2o_v52.hdf", "h2o_sco2_v52_f32.absco"),
)

const HAS_CUDA = Base.find_package("CUDA") !== nothing
HAS_CUDA && @eval using CUDA
const CUDA_OK = HAS_CUDA && CUDA.functional()

@inline function bracket(nodes, x)
    length(nodes) == 1 && return 1, 1, 0.0
    i = clamp(searchsortedlast(nodes, x), 1, length(nodes) - 1)
    return i, i + 1, (x - nodes[i]) / (nodes[i + 1] - nodes[i])
end

# Independent scalar implementation of the documented four-axis linear interpolation. It uses
# Float64 arithmetic over stored Float32 values and does not call package interpolation internals.
function linear_reference(lut, x, pressure, temperature, water_vmr)
    ν, p, T, vmr, σ = lut.ν, lut.p, lut.T, lut.vmr, lut.σ
    x32, p32 = Float32(x), Float32(pressure)
    t32, v32 = Float32(temperature), Float32(water_vmr)
    (x32 < first(ν) || x32 > last(ν)) && return 0.0
    iν, iν1, fν = bracket(ν, x32)
    ip, ip1, fp = bracket(p, clamp(p32, first(p), last(p)))
    iv, iv1, fv = bracket(vmr, clamp(v32, first(vmr), last(vmr)))

    at_pressure(jp) = begin
        taxis = @view T[:, jp]
        it, it1, ft = bracket(taxis, clamp(t32, first(taxis), last(taxis)))
        at_wavenumber(jν) = begin
            low = (1 - ft) * Float64(σ[jν, iv, it, jp]) +
                  ft * Float64(σ[jν, iv, it1, jp])
            iv == iv1 && return low
            high = (1 - ft) * Float64(σ[jν, iv1, it, jp]) +
                   ft * Float64(σ[jν, iv1, it1, jp])
            return (1 - fv) * low + fv * high
        end
        return (1 - fν) * at_wavenumber(iν) + fν * at_wavenumber(iν1)
    end
    low = at_pressure(ip)
    ip == ip1 && return low
    return (1 - fp) * low + fp * at_pressure(ip1)
end

function raw_node(path, lut, iν, iv, it, ip)
    NCDataset(path) do ds
        νsource = Array(ds["Wavenumber"])
        source_index = argmin(abs.(νsource .- Float64(lut.ν[iν])))
        variable = only(filter(name -> startswith(name, "Gas_") &&
                                       endswith(name, "_Absorption"), collect(keys(ds))))
        return Float32(ds[variable][source_index, iv, it, ip])
    end
end

function median_query_ms(lut, grid; samples=25)
    ip = cld(length(lut.p), 2)
    p, T = lut.p[ip], lut.T[cld(size(lut.T, 1), 2), ip]
    compute_cross_section(lut, grid, p, T; vmr=0.03f0, interp=:linear)
    times = [@elapsed compute_cross_section(lut, grid, p, T;
                                             vmr=0.03f0, interp=:linear) for _ in 1:samples]
    return 1e3 * median(times)
end

function temperature_loo(lut)
    # This is stricter than production use: the operational LUT keeps every temperature node.
    σ = lut.σ
    relative_rms, transmission_rms = 0.0, 0.0
    pressures = unique(round.(Int, range(1, length(lut.p), length=7)))
    for ip in pressures, iv in eachindex(lut.vmr), it in 2:2:(size(lut.T, 1) - 1)
        truth = Float64.(@view σ[:, iv, it, ip])
        predicted = 0.5 .* (Float64.(@view σ[:, iv, it - 1, ip]) .+
                            Float64.(@view σ[:, iv, it + 1, ip]))
        norm(truth) == 0 && continue
        relative_rms = max(relative_rms, norm(predicted - truth) / norm(truth))
        column = 1 / maximum(truth)
        transmission_rms = max(transmission_rms,
            sqrt(mean((exp.(-column .* predicted) .- exp.(-column .* truth)).^2)))
    end
    return relative_rms, transmission_rms
end

failures = String[]
println("ABSCO v5.2 Float32 OCO-2 LUT validation; CUDA parity: ",
        CUDA_OK ? "enabled" : "not available")
for (source_name, lut_name) in PRODUCTS
    source_path, lut_path = joinpath(ABSCO_ROOT, source_name), joinpath(LUT_ROOT, lut_name)
    load_seconds = @elapsed lut = load_absco_lut(lut_path; architecture=CPU())
    iν, iv, it, ip = cld(length(lut.ν), 2), 2, 9, 32
    node_ok = raw_node(source_path, lut, iν, iv, it, ip) == lut.σ[iν, iv, it, ip]
    type_ok = eltype(lut) === Float32 && eltype(lut.σ) === Float32

    query_grid = collect(Float32, range(first(lut.ν) + 0.005f0,
                                        last(lut.ν) - 0.005f0, length=257))
    pressure = 0.37f0 * lut.p[21] + 0.63f0 * lut.p[22]
    temperature = 0.41f0 * lut.T[8, 21] + 0.59f0 * lut.T[9, 21]
    water_vmr = 0.041f0
    actual = Array(compute_cross_section(lut, query_grid, pressure, temperature;
                                         vmr=water_vmr, interp=:linear))
    reference = [linear_reference(lut, x, pressure, temperature, water_vmr)
                 for x in query_grid]
    reference_error = norm(Float64.(actual) - reference) / max(norm(reference), eps())
    native_ms = median_query_ms(lut, lut.ν)
    instrument_grid = collect(Float32, range(first(lut.ν), last(lut.ν), length=2048))
    instrument_ms = median_query_ms(lut, instrument_grid)
    loo_xsec, loo_transmission = temperature_loo(lut)

    gpu_error, gpu_ms = NaN, NaN
    if CUDA_OK
        gpu = load_absco_lut(lut_path; architecture=GPU())
        gpu_actual = Array(compute_cross_section(gpu, query_grid, pressure, temperature;
                                                 vmr=water_vmr, interp=:linear))
        gpu_error = norm(Float64.(gpu_actual) - Float64.(actual)) /
                    max(norm(Float64.(actual)), eps())
        gpu_ms = median_query_ms(gpu, instrument_grid)
        gpu = nothing
        CUDA.reclaim()
    end

    node_ok || push!(failures, "$lut_name: raw-node mismatch")
    type_ok || push!(failures, "$lut_name: not pure Float32")
    reference_error <= 2e-5 || push!(failures, "$lut_name: independent linear error")
    (!CUDA_OK || gpu_error <= 2e-5) || push!(failures, "$lut_name: CPU/GPU mismatch")
    loo_transmission <= 0.01 || push!(failures, "$lut_name: T-LOO transmission RMS > 1%")

    @printf("%-27s shape=%-24s load=%5.2fs CPU(native/2048)=%7.3f/%7.3fms ",
            lut_name, string(size(lut.σ)), load_seconds, native_ms, instrument_ms)
    @printf("node=%s linear=%.2e T-LOO(xsec/trans)=%.3f/%.4f",
            node_ok ? "exact" : "FAIL", reference_error, loo_xsec, loo_transmission)
    CUDA_OK && @printf(" GPU2048=%7.3fms parity=%.2e", gpu_ms, gpu_error)
    println()
    lut = nothing
    GC.gc()
end

isempty(failures) || error("ABSCO validation failed:\n  " * join(failures, "\n  "))
println("PASS: all operational gates satisfied")

using Test
using VofiJul

const CONV_RADIUS = 0.35
const CONV_DOMAIN_START = -0.5
const CONV_DOMAIN_SIZE = 1.0

function integrate_volume_over_grid(func, n)
    h = CONV_DOMAIN_SIZE / n
    spacing = (h, h, h)
    xex_tmp = zeros(Float64, 4)
    nex = [0, 0]
    npt = zeros(Int, 4)
    nvis = [0, 0]
    total = 0.0
    for i in 0:n-1, j in 0:n-1, k in 0:n-1
        xin = [CONV_DOMAIN_START + i * h,
               CONV_DOMAIN_START + j * h,
               CONV_DOMAIN_START + k * h]
        cc = vofi_get_cc(func, nothing, xin, spacing, xex_tmp, nex, npt, nvis, 3)
        total += cc
    end
    return total * h^3
end

function write_convergence_csv(path, rows)
    open(path, "w") do io
        println(io, "N,h,volume_numerical,volume_exact,abs_error,rel_error,observed_order")
        for r in rows
            println(io, string(r.N, ",", r.h, ",", r.volume_num, ",", r.volume_exact,
                               ",", r.abs_error, ",", r.rel_error, ",", r.observed_order))
        end
    end
    return path
end

function run_sphere_convergence(; resolutions = [6, 8, 12, 16, 24, 32],
                                 radius = CONV_RADIUS,
                                 outfile = joinpath(@__DIR__, "convergence_3d_sphere.csv"))
    function sphere_sdf(x, _)
        return sqrt(x[1]^2 + x[2]^2 + x[3]^2) - radius
    end
    exact = (4.0 / 3.0) * π * radius^3
    results = NamedTuple[]
    prev_err = nothing
    prev_h = nothing
    for n in resolutions
        h = CONV_DOMAIN_SIZE / n
        vol_num = integrate_volume_over_grid(sphere_sdf, n)
        abs_err = abs(vol_num - exact)
        rel_err = abs_err / exact
        order = if prev_err === nothing || prev_err <= 0 || abs_err <= 0
            missing
        else
            log(prev_err / abs_err) / log(prev_h / h)
        end
        push!(results, (; N = n, h = h, volume_num = vol_num, volume_exact = exact,
                         abs_error = abs_err, rel_error = rel_err, observed_order = order))
        prev_err = abs_err
        prev_h = h
    end
    write_convergence_csv(outfile, results)
    return results, outfile
end

@testset "3D sphere convergence sweep" begin
    results, outfile = run_sphere_convergence()
    errors = [r.abs_error for r in results]
#    @test issorted(errors; rev=true)
    @test results[end].rel_error < 0.05
    @test isfile(outfile)
end

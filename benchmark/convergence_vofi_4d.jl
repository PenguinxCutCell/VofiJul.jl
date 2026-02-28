using Test
using VofiJul

const CONV4D_DOMAIN_START = -0.5
const CONV4D_DOMAIN_SIZE = 1.0

function integrate_4d_grid(func; n::Int)
    h = CONV4D_DOMAIN_SIZE / n
    spacing = fill(h, 4)
    xex = zeros(Float64, 5)
    nex = [0, 0]
    npt = zeros(Int, 4)
    nvis = [0, 0]
    total = 0.0
    for i in 0:n-1, j in 0:n-1, k in 0:n-1, l in 0:n-1
        xin = [CONV4D_DOMAIN_START + i * h,
               CONV4D_DOMAIN_START + j * h,
               CONV4D_DOMAIN_START + k * h,
               CONV4D_DOMAIN_START + l * h]
        cc = vofi_get_cc(func, nothing, xin, spacing, xex, nex, npt, nvis, 4)
        total += cc
    end
    return total * h^4
end

const SPHERE_RADIUS_4D = 0.35
function hypersphere_sdf_4d(x, _)
    return sqrt(x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2) - SPHERE_RADIUS_4D
end

const HYPERCUBE_HALF_4D = 0.3
function hypercube_sdf_4d(x, _)
    return maximum(abs.(x)) - HYPERCUBE_HALF_4D
end

const ELLIPSOID_AXES_4D = (0.45, 0.3, 0.25, 0.5)
function ellipsoid_sdf_4d(x, _)
    return (x[1] / ELLIPSOID_AXES_4D[1])^2 +
           (x[2] / ELLIPSOID_AXES_4D[2])^2 +
           (x[3] / ELLIPSOID_AXES_4D[3])^2 +
           (x[4] / ELLIPSOID_AXES_4D[4])^2 - 1.0
end

const SIN_OFFSET_4D = 0.1
const SIN_AMPLITUDE_4D = 0.15
function sinuslab_sdf_4d(x, _)
    oscillation = sinpi(2 * x[1]) * sinpi(2 * x[2]) * sinpi(2 * x[3])
    target = SIN_OFFSET_4D + SIN_AMPLITUDE_4D * oscillation
    return x[4] - target
end

const CONV4D_CASES = [
    (name = "hypersphere", sdf = hypersphere_sdf_4d,
     exact = (π^2 / 2) * SPHERE_RADIUS_4D^4, rel_tol = 0.04),
    (name = "ellipsoid", sdf = ellipsoid_sdf_4d,
     exact = (π^2 / 2) * prod(ELLIPSOID_AXES_4D), rel_tol = 0.05),
    (name = "sinusoidal", sdf = sinuslab_sdf_4d,
     exact = 0.5 + SIN_OFFSET_4D, rel_tol = 0.02),
]

function write_convergence_csv(path, rows)
    open(path, "w") do io
        println(io, "shape,N,h,volume_numerical,volume_exact,abs_error,rel_error,observed_order")
        for r in rows
            ord_str = ismissing(r.observed_order) ? "" : string(r.observed_order)
            println(io, string(r.shape, ",", r.N, ",", r.h, ",", r.volume_num, ",",
                               r.volume_exact, ",", r.abs_error, ",", r.rel_error, ",", ord_str))
        end
    end
    return path
end

function run_convergence_4d(; resolutions = [4, 8, 16],
                              cases = CONV4D_CASES,
                              outfile = joinpath(@__DIR__, "convergence_4d.csv"))
    rows = NamedTuple[]
    prev = Dict{String, Tuple{Union{Missing, Float64}, Union{Missing, Float64}}}()
    for c in cases
        prev[c.name] = (missing, missing)
    end
    for n in resolutions
        println("Running 4D convergence tests at N = $(n)...")
        h = CONV4D_DOMAIN_SIZE / n
        for c in cases
            println("  Integrating shape: $(c.name)")
            vol = integrate_4d_grid(c.sdf; n = n)
            abs_err = abs(vol - c.exact)
            rel_err = abs_err / c.exact
            ph, pe = prev[c.name]
            order = (ismissing(pe) || pe <= 0.0 || abs_err <= 0.0) ? missing :
                    log(pe / abs_err) / log(ph / h)
            push!(rows, (; shape = c.name, N = n, h = h, volume_num = vol,
                          volume_exact = c.exact, abs_error = abs_err,
                          rel_error = rel_err, observed_order = order))
            prev[c.name] = (h, abs_err)
        end
    end
    write_convergence_csv(outfile, rows)
    return rows, outfile
end

@testset "4D convergence sweeps" begin
    rows, outfile = run_convergence_4d()
    @test isfile(outfile)
    for c in CONV4D_CASES
        errs = [r.abs_error for r in rows if r.shape == c.name]
        rels = [r.rel_error for r in rows if r.shape == c.name]
        @test !isempty(errs)
#        @test issorted(errs; rev = true)
        @test rels[end] < c.rel_tol
    end
end

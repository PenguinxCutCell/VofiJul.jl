using BenchmarkTools
using VofiJul
using Vofinit
HAVE_VOFINIT = true

# 2D circle SDF used for the benchmark
const R2 = 0.2
function circle_sdf_vofi(x, _)
    return sqrt(x[1]^2 + x[2]^2) - R2
end

# Vofinit wrappers (accept both 2-arg and 4-arg forms to match possible C wrappers)
if HAVE_VOFINIT
    function circle_vofinit_func(x::Float64, y::Float64)
        return sqrt(x^2 + y^2) - R2
    end
end

function integrate_with_vofijul_2d(N)
    h = 1.0 / N
    area = 0.0
    xex = zeros(Float64,4)
    nex = [0,0]
    npt = zeros(Int,4)
    nvis = [0,0]
    for i in 0:N-1, j in 0:N-1
        xin = [ -0.5 + i*h, -0.5 + j*h ]
        cc = vofi_get_cc(circle_sdf_vofi, nothing, xin, [h,h], xex, nex, npt, nvis, 2)
        area += cc
    end
    return area * h^2
end

function integrate_with_vofinit_2d(N)
    if !HAVE_VOFINIT
        error("Vofinit not available")
    end
    h = 1.0 / N
    area = 0.0
    xex = zeros(Float64,4)
    for i in 0:N-1, j in 0:N-1
        xin = [ -0.5 + i*h, -0.5 + j*h ]
        xin4 = [xin[1], xin[2]]
        h4 = [h, h]
        cc = Vofinit.getcc(circle_vofinit_func, xin4, h4, xex, Cint(2))
        area += cc
    end
    return area * h^2
end

function run_bench_2d(N=64)
    println("Benchmark (2D): VofiJul vofi_get_cc integrate over $(N)^2 cells")
    b1 = @benchmark integrate_with_vofijul_2d($N) samples=10
    println(b1)
    println("Allocation test: VofiJul integrate_with_vofijul_2d($(N))")
    @btime integrate_with_vofijul_2d($N)

    if HAVE_VOFINIT
        println("Benchmark (2D): Vofinit.getcc integrate over $(N)^2 cells")
        try
            b2 = @benchmark integrate_with_vofinit_2d($N)
            println(b2)
            println("Allocation test: Vofinit integrate_with_vofinit_2d($(N))")
            @btime integrate_with_vofinit_2d($N)
        catch err
            @warn "Vofinit 2D benchmark failed" err=err
        end
    end
end

run_bench_2d(16)
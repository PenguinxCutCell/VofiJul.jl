using BenchmarkTools
using VofiJul
using Vofinit

HAVE_VOFINIT = true

# Sphere SDF used for the benchmark
const R = 0.4
function sphere_sdf_vofi(x, _)
    return sqrt(x[1]^2 + x[2]^2 + x[3]^2) - R
end

# For Vofinit (if available) define a function matching what the wrapper expects.
# This may need to be adapted to the exact wrapper signature; we keep a best-effort wrapper.
if HAVE_VOFINIT
    # Provide a 4-arg function in case Vofinit expects (x,y,z,w)
    function sphere_vofinit_func(x::Float64, y::Float64, z::Float64, w::Float64)
        return sqrt(x^2 + y^2 + z^2) - R
    end
    # Provide a 3-arg wrapper in case the C wrapper calls with only (x,y,z)
    function sphere_vofinit_func(x::Float64, y::Float64, z::Float64)
        return sphere_vofinit_func(x, y, z, 0.0)
    end
end

function integrate_with_vofijul(N)
    h = 1.0 / N
    vol = 0.0
    xex = zeros(Float64,4)
    nex = [0,0]
    npt = zeros(Int,4)
    nvis = [0,0]
    for i in 0:N-1, j in 0:N-1, k in 0:N-1
        xin = [ -0.5 + i*h, -0.5 + j*h, -0.5 + k*h ]
        cc = vofi_get_cc(sphere_sdf_vofi, nothing, xin, [h,h,h], xex, nex, npt, nvis, 3)
        vol += cc
    end
    return vol * h^3
end

function integrate_with_vofinit(N)
    if !HAVE_VOFINIT
        error("Vofinit not available")
    end
    h = 1.0 / N
    vol = 0.0
    xex = zeros(Float64,4)
    # call pattern: Vofinit.getcc(func, x0, h0, xex, ndim0)
    # note: the wrapper may expect different argument types; this is a best-effort call
    cc = 0.0
    for i in 0:N-1, j in 0:N-1, k in 0:N-1
        xin = [ -0.5 + i*h, -0.5 + j*h, -0.5 + k*h ]
        try
            cc = Vofinit.getcc(sphere_vofinit_func, xin, [h,h,h], xex, Cint(3))
        catch err
            # if the wrapper fails, rethrow so the user can inspect
            rethrow(err)
        end
        vol += cc
    end
    return vol * h^3
end

function run_bench(N=16)
    println("Benchmark: VofiJul vofi_get_cc integrate over $(N)^3 cells")
    b1 = @benchmark integrate_with_vofijul($N) samples=10
    println(b1)
    if HAVE_VOFINIT
        println("Benchmark: Vofinit.getcc integrate over $(N)^3 cells")
        try
            b2 = @benchmark integrate_with_vofinit($N) samples=10
            println(b2)
        catch err
            @warn "Vofinit benchmark failed" err=err
        end
    end
end

run_bench(16)
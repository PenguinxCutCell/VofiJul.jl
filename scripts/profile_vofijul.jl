#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using VofiJul
using Profile

function build_case(case::Symbol)
    if case == :getcc3d
        function sphere_sdf(x, _)
            r = 0.4
            return sqrt(x[1]^2 + x[2]^2 + x[3]^2) - r
        end
        xin = [-0.5, -0.5, -0.5]
        h = [1.0, 1.0, 1.0]
        xex = zeros(Float64, 4)
        nex = [0,0]
        npt = [0,0,0,0]
        nvis = [0,0]
        return () -> vofi_get_cc(sphere_sdf, nothing, xin, h, xex, nex, npt, nvis, 3)

    elseif case == :getcc2d
        function circle_sdf(x, _)
            r = 0.4
            return sqrt(x[1]^2 + x[2]^2) - r
        end
        xin = [-0.5, -0.5]
        h = [1.0, 1.0]
        xex = zeros(Float64, 4)
        nex = [0,0]
        npt = [0,0,0,0]
        nvis = [0,0]
        return () -> vofi_get_cc(circle_sdf, nothing, xin, h, xex, nex, npt, nvis, 2)

    elseif case == :getcc4d
        function hypersphere_sdf(x, _)
            r = 0.35
            return sqrt(x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2) - r
        end
        xin = [-0.5, -0.5, -0.5, -0.5]
        h = [1.0, 1.0, 1.0, 1.0]
        xex = zeros(Float64, 5)
        nex = [0,0]
        npt = [0,0,0,0]
        nvis = [0,0]
        return () -> vofi_get_cc(hypersphere_sdf, nothing, xin, h, xex, nex, npt, nvis, 4)
        
    elseif case == :interface3d
        function sphere_sdf(x, _)
            r = 0.4
            return sqrt(x[1]^2 + x[2]^2 + x[3]^2) - r
        end
        xin = [-0.5, -0.5, -0.5]
        h = [1.0, 1.0, 1.0]
        return () -> vofi_interface_centroid(sphere_sdf, nothing, xin, h, 3)

    elseif case == :celltype3d
        function hyperplane(x, _)
            return x[1] + x[2] + x[3] - 0.5
        end
        xin = [0.0, 0.0, 0.0]
        h = [1.0, 1.0, 1.0]
        return () -> vofi_get_cell_type(hyperplane, nothing, xin, h, 3)

    elseif case == :gauss
        return () -> begin
            nodes = gauss_legendre_nodes(8)
            weights = gauss_legendre_weights(8)
            return (nodes, weights)
        end

    else
        error("Unknown case: $(case). Use one of :getcc3d, :getcc2d, :interface3d, :celltype3d, :gauss.")
    end
end

function benchmark_profile(; case::Symbol=:getcc3d, n_eval::Int=2000, profile_n::Int=800)
    f = build_case(case)

    # Warm-up compilation
    f()

    elapsed = @elapsed begin
        for _ in 1:n_eval
            f()
        end
    end

    total_alloc_bytes = @allocated begin
        for _ in 1:n_eval
            f()
        end
    end

    println("Case: $(case)")
    println("Runs: $(n_eval)")
    println("Total time (s): ", elapsed)
    println("Avg time (us/run): ", elapsed * 1e6 / n_eval)
    println("Total allocated (bytes): ", total_alloc_bytes)
    println("Avg allocated (bytes/run): ", total_alloc_bytes / n_eval)

    Profile.clear()
    @profile begin
        for _ in 1:profile_n
            f()
        end
    end

    println("\nCPU profile (flat, by sample count):")
    Profile.print(format=:flat, sortedby=:count, mincount=5)

    return nothing
end

case = length(ARGS) >= 1 ? Symbol(ARGS[1]) : :getcc3d
n_eval = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 2000
profile_n = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 800

benchmark_profile(; case=case, n_eval=n_eval, profile_n=profile_n)

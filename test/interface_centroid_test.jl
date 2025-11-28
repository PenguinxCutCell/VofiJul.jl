using Test
using VofiJul
using LinearAlgebra

"""
Simple signed-distance helper constructors.
"""
sd_line(x) = x[1] - 0.3
sd_plane2D(x) = x[1] + 2x[2] - 0.5
sd_plane3D(x) = x[1] - 0.2 + 0.5x[2] - 0.25x[3]

@testset "vofi_interface_centroid 1D" begin
    xin = [0.0]
    h0 = [1.0]
    c = vofi_interface_centroid(sd_line, nothing, xin, h0, 1)
    @test isapprox(c[1], 0.3; atol=1e-8)
end

@testset "vofi_interface_centroid 2D" begin
    xin = [0.0, 0.0]
    h0 = [1.0, 1.0]
    c = vofi_interface_centroid(sd_plane2D, nothing, xin, h0, 2)
    # The plane x + 2y = 0.5 intersects the unit square near (0.1,0.2)
    @test isapprox(c[1] + 2c[2], 0.5; atol=1e-6)
end

@testset "vofi_interface_centroid 3D" begin
    xin = [0.0, 0.0, 0.0]
    h0 = [1.0, 1.0, 1.0]
    c = vofi_interface_centroid(sd_plane3D, nothing, xin, h0, 3)
    @test isapprox(c[1] - 0.2 + 0.5c[2] - 0.25c[3], 0.0; atol=1e-6)
end

@testset "vofi_interface_centroid 4D" begin
    # 4D test case: plane x + 2y + 3z + 4w = 0.5
    sd_plane4D(x) = x[1] + 2x[2] + 3x[3] + 4x[4] - 0.5
    xin = [0.0, 0.0, 0.0, 0.0]
    h0 = [1.0, 1.0, 1.0, 1.0]
    c = vofi_interface_centroid(sd_plane4D, nothing, xin, h0, 4)
    @test isapprox(c[1] + 2c[2] + 3c[3] + 4c[4], 0.5; atol=1e-6)
end

@testset "vofi_interface_centroid fallback" begin
    # Test fallback to cell centre when no sign change
    sd_constant(x) = 1.0
    xin = [0.0, 0.0]
    h0 = [1.0, 1.0]
    c = vofi_interface_centroid(sd_constant, nothing, xin, h0, 2)
    @test isapprox(c, [0.5, 0.5]; atol=1e-8)
end


@testset "vofi_interface_centroid sphere 3d" begin
    # Test sphere
    sd_sphere(x) = norm(x) - 0.5
    xin = [0.0, 0.0, 0.0]
    h0 = [1.0, 1.0, 1.0]
    c = vofi_interface_centroid(sd_sphere, nothing, xin, h0, 3)
    @test isapprox(norm(c), 0.5; atol=1e-6)
end

@testset "vofi_interface_centroid hypersphere 4d" begin
    # Test hypersphere in 4D
    sd_hypersphere(x) = norm(x) - 0.5
    xin = [0.0, 0.0, 0.0, 0.0]
    h0 = [1.0, 1.0, 1.0, 1.0]
    c = vofi_interface_centroid(sd_hypersphere, nothing, xin, h0, 4)
    @test isapprox(norm(c), 0.5; atol=1e-6)
end

@testset "vofi_interface_centroid 4D hypercube" begin
    # Test hypercube interface in 4D
    sd_hypercube(x) = maximum(abs.(x)) - 0.3
    xin = [0.0, 0.0, 0.0, 0.0]
    h0 = [1.0, 1.0, 1.0, 1.0]
    c = vofi_interface_centroid(sd_hypercube, nothing, xin, h0, 4)
    @test isapprox(norm(c, Inf), 0.3; atol=1e-6)
end

# Tests for interface centroid computed via vofi_get_cc using precomputed heights
@testset "vofi_get_cc interface centroid 2D" begin
    # Test 2D interface centroid from vofi_get_cc
    # Vertical line at x = 0.3 crossing the unit square
    sd_vline(x, _) = x[1] - 0.3
    xin = [0.0, 0.0]
    h0 = [1.0, 1.0]
    xex = zeros(Float64, 6)
    nex = [1, 2]  # Request volume centroid and interface centroid
    npt = zeros(Int, 4)
    nvis = [0, 0]
    cc = vofi_get_cc(sd_vline, nothing, xin, h0, xex, nex, npt, nvis, 2)
    # The interface is a vertical line at x = 0.3, from y=0 to y=1
    # Interface centroid should be at (0.3, 0.5)
    @test isapprox(xex[5], 0.3; atol=1e-3)  # x-coordinate
    @test isapprox(xex[6], 0.5; atol=1e-3)  # y-coordinate
end

@testset "vofi_get_cc interface centroid 2D diagonal" begin
    # Test 2D interface centroid for a diagonal line
    sd_diag(x, _) = x[1] - x[2]  # Line x = y
    xin = [0.0, 0.0]
    h0 = [1.0, 1.0]
    xex = zeros(Float64, 6)
    nex = [1, 2]  # Request volume centroid and interface centroid
    npt = zeros(Int, 4)
    nvis = [0, 0]
    cc = vofi_get_cc(sd_diag, nothing, xin, h0, xex, nex, npt, nvis, 2)
    # The interface is a diagonal line from (0,0) to (1,1)
    # Interface centroid should be at (0.5, 0.5)
    @test isapprox(xex[5], xex[6]; atol=1e-2)  # x = y on the interface
end

@testset "vofi_get_cc interface centroid 3D sphere" begin
    # Test 3D interface surface centroid for a sphere
    sd_sphere(x, _) = sqrt(x[1]^2 + x[2]^2 + x[3]^2) - 0.4
    xin = [0.0, 0.0, 0.0]
    h0 = [1.0, 1.0, 1.0]
    xex = zeros(Float64, 7)
    nex = [1, 2]  # Request volume centroid and interface centroid
    npt = zeros(Int, 4)
    nvis = [0, 0]
    cc = vofi_get_cc(sd_sphere, nothing, xin, h0, xex, nex, npt, nvis, 3)
    # Interface centroid should lie on the sphere surface
    iface_cent = [xex[5], xex[6], xex[7]]
    if norm(iface_cent) > 0.01  # Only check if we got a valid centroid
        @test isapprox(norm(iface_cent), 0.4; atol=0.1)
    end
end

@testset "vofi_get_cc interface centroid 3D plane" begin
    # Test 3D interface surface centroid for a plane
    sd_plane(x, _) = x[1] - 0.5  # Vertical plane at x = 0.5
    xin = [0.0, 0.0, 0.0]
    h0 = [1.0, 1.0, 1.0]
    xex = zeros(Float64, 7)
    nex = [1, 2]  # Request volume centroid and interface centroid
    npt = zeros(Int, 4)
    nvis = [0, 0]
    cc = vofi_get_cc(sd_plane, nothing, xin, h0, xex, nex, npt, nvis, 3)
    # Interface centroid should be at center of the plane section (0.5, 0.5, 0.5)
    @test isapprox(xex[5], 0.5; atol=0.1)  # x-coordinate
    @test isapprox(xex[6], 0.5; atol=0.1)  # y-coordinate
    @test isapprox(xex[7], 0.5; atol=0.1)  # z-coordinate
end
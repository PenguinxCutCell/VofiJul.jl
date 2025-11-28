using Test
using VofiJul

# Test sphere SDF for 3D
sphere_sdf(x, _) = sqrt(x[1]^2 + x[2]^2 + x[3]^2) - 0.4

# Test circle SDF for 2D
circle_sdf(x, _) = sqrt(x[1]^2 + x[2]^2) - 0.4

# Test line SDF for 1D
line_sdf(x, _) = abs(x[1]) - 0.4

@testset "vofi_get_cc_threaded" begin
    # 3D test - sphere
    n = 10
    h = 1.0 / n
    cells_3d = [[x, y, z] for x in -0.5:h:0.5-h for y in -0.5:h:0.5-h for z in -0.5:h:0.5-h]
    
    cc_values = vofi_get_cc_threaded(sphere_sdf, nothing, cells_3d, [h, h, h], 3)
    
    @test length(cc_values) == n^3
    total_vol = sum(cc_values) * h^3
    exact_vol = 4/3 * π * 0.4^3
    @test total_vol ≈ exact_vol atol=0.02
    
    # 2D test - circle
    cells_2d = [[x, y] for x in -0.5:h:0.5-h for y in -0.5:h:0.5-h]
    cc_values_2d = vofi_get_cc_threaded(circle_sdf, nothing, cells_2d, [h, h], 2)
    
    @test length(cc_values_2d) == n^2
    total_area = sum(cc_values_2d) * h^2
    exact_area = π * 0.4^2
    @test total_area ≈ exact_area atol=0.02
    
    # 1D test - line
    cells_1d = [[x] for x in -0.5:h:0.5-h]
    cc_values_1d = vofi_get_cc_threaded(line_sdf, nothing, cells_1d, [h], 1)
    
    @test length(cc_values_1d) == n
    total_length = sum(cc_values_1d) * h
    exact_length = 2 * 0.4
    @test total_length ≈ exact_length atol=0.02
end

@testset "vofi_integrate_threaded" begin
    # 3D sphere integration
    result_3d = vofi_integrate_threaded(sphere_sdf, nothing, 
                                        [-0.5, -0.5, -0.5], [0.5, 0.5, 0.5],
                                        10, 3)
    exact_vol = 4/3 * π * 0.4^3
    @test result_3d.total ≈ exact_vol atol=0.02
    @test size(result_3d.cc_values) == (10, 10, 10)
    @test result_3d.cell_volume ≈ (1.0/10)^3
    
    # 2D circle integration
    result_2d = vofi_integrate_threaded(circle_sdf, nothing,
                                        [-0.5, -0.5], [0.5, 0.5],
                                        10, 2)
    exact_area = π * 0.4^2
    @test result_2d.total ≈ exact_area atol=0.02
    @test size(result_2d.cc_values) == (10, 10)
    
    # 1D line integration
    result_1d = vofi_integrate_threaded(line_sdf, nothing,
                                        [-0.5], [0.5],
                                        10, 1)
    exact_length = 2 * 0.4
    @test result_1d.total ≈ exact_length atol=0.02
    @test size(result_1d.cc_values) == (10,)
end

@testset "vofi_get_cell_type_threaded" begin
    # Test cells that are fully inside, outside, or cut
    cells = [
        [0.0, 0.0, 0.0],        # center - likely cut
        [-0.45, -0.45, -0.45],  # corner - likely outside
        [0.1, 0.1, 0.1],        # near center - likely inside
    ]
    h = [0.1, 0.1, 0.1]
    
    types = vofi_get_cell_type_threaded(sphere_sdf, nothing, cells, h, 3)
    
    @test length(types) == 3
    @test all(t -> t in (-1, 0, 1), types)  # Valid cell types
end

@testset "Results consistency" begin
    # Verify threaded and non-threaded give same results
    cells = [[x, y] for x in -0.3:0.1:0.2 for y in -0.3:0.1:0.2]
    h = [0.1, 0.1]
    
    # Non-threaded
    cc_sequential = [begin
        xex = zeros(Float64, 4)
        vofi_get_cc(circle_sdf, nothing, c, h, xex, [0, 0], zeros(Int, 4), [0, 0], 2)
    end for c in cells]
    
    # Threaded
    cc_threaded = vofi_get_cc_threaded(circle_sdf, nothing, cells, h, 2)
    
    @test length(cc_sequential) == length(cc_threaded)
    @test all(cc_sequential .≈ cc_threaded)
end

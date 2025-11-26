using Test
using VofiJul

@testset "Interface centroid tests" begin
    # Test 1D: Interface position
    @testset "1D interface centroid" begin
        # Simple line at x = 0.25
        line_sdf(x, _) = x[1] - 0.25
        measure, centroid = vofi_get_interface_centroid(line_sdf, nothing, [0.0], [1.0], 1)
        @test measure == 0.0  # Point has zero measure
        @test centroid[1] ≈ 0.25 atol=1e-10
        
        # Line at x = 0.75
        line_sdf2(x, _) = x[1] - 0.75
        measure2, centroid2 = vofi_get_interface_centroid(line_sdf2, nothing, [0.0], [1.0], 1)
        @test centroid2[1] ≈ 0.75 atol=1e-10
        
        # No interface (fully inside)
        neg_func(x, _) = -1.0
        measure_neg, centroid_neg = vofi_get_interface_centroid(neg_func, nothing, [0.0], [1.0], 1)
        @test measure_neg == 0.0
        @test centroid_neg[1] == 0.0
        
        # No interface (fully outside)
        pos_func(x, _) = 1.0
        measure_pos, centroid_pos = vofi_get_interface_centroid(pos_func, nothing, [0.0], [1.0], 1)
        @test measure_pos == 0.0
        @test centroid_pos[1] == 0.0
    end
    
    # Test 2D: Interface length and centroid
    @testset "2D interface centroid" begin
        # Vertical line at x = 0.5
        vline_sdf(x, _) = x[1] - 0.5
        measure, centroid = vofi_get_interface_centroid(vline_sdf, nothing, [0.0, 0.0], [1.0, 1.0], 2)
        @test measure ≈ 1.0 atol=1e-6  # Line length = 1
        @test centroid[1] ≈ 0.5 atol=1e-6  # x-coordinate of centroid
        @test centroid[2] ≈ 0.5 atol=1e-6  # y-coordinate at midpoint of line
        
        # Horizontal line at y = 0.5
        hline_sdf(x, _) = x[2] - 0.5
        measure_h, centroid_h = vofi_get_interface_centroid(hline_sdf, nothing, [0.0, 0.0], [1.0, 1.0], 2)
        @test measure_h ≈ 1.0 atol=1e-6  # Line length = 1
        @test centroid_h[1] ≈ 0.5 atol=1e-6  # x-coordinate at midpoint
        @test centroid_h[2] ≈ 0.5 atol=1e-6  # y-coordinate of line
        
        # Circle centered at origin - centroid should be at (0, 0)
        r = 0.4
        circle_sdf(x, _) = sqrt(x[1]^2 + x[2]^2) - r
        # Test in a single cut cell containing part of the circle
        measure_c, centroid_c = vofi_get_interface_centroid(circle_sdf, nothing, [0.0, 0.0], [0.5, 0.5], 2)
        @test measure_c > 0  # Should have interface length
        # For circle arc in first quadrant, centroid should be in first quadrant
        @test centroid_c[1] >= 0
        @test centroid_c[2] >= 0
        
        # No interface (fully inside)
        neg_func2d(x, _) = -1.0
        measure_neg, centroid_neg = vofi_get_interface_centroid(neg_func2d, nothing, [0.0, 0.0], [1.0, 1.0], 2)
        @test measure_neg == 0.0
        @test all(centroid_neg .== 0.0)
    end
    
    # Test 3D: Interface area and centroid
    @testset "3D interface centroid" begin
        # Plane x = 0.5 
        plane_sdf(x, _) = x[1] - 0.5
        measure, centroid = vofi_get_interface_centroid(plane_sdf, nothing, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 3)
        @test measure ≈ 1.0 atol=1e-6  # Plane area = 1
        @test centroid[1] ≈ 0.5 atol=1e-6  # x-coordinate of plane
        @test centroid[2] ≈ 0.5 atol=1e-6  # y-coordinate at center
        @test centroid[3] ≈ 0.5 atol=1e-6  # z-coordinate at center
        
        # Plane y = 0.5
        plane_y_sdf(x, _) = x[2] - 0.5
        measure_y, centroid_y = vofi_get_interface_centroid(plane_y_sdf, nothing, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 3)
        @test measure_y ≈ 1.0 atol=1e-6  # Plane area = 1
        @test centroid_y[1] ≈ 0.5 atol=1e-6
        @test centroid_y[2] ≈ 0.5 atol=1e-6
        @test centroid_y[3] ≈ 0.5 atol=1e-6
        
        # Sphere centered at origin
        r = 0.4
        sphere_sdf(x, _) = sqrt(x[1]^2 + x[2]^2 + x[3]^2) - r
        # Test in a single cut cell
        measure_s, centroid_s = vofi_get_interface_centroid(sphere_sdf, nothing, [0.0, 0.0, 0.0], [0.5, 0.5, 0.5], 3)
        @test measure_s > 0  # Should have interface area
        # For sphere surface in first octant, centroid should be in first octant
        @test centroid_s[1] >= 0
        @test centroid_s[2] >= 0
        @test centroid_s[3] >= 0
        
        # No interface (fully inside)
        neg_func3d(x, _) = -1.0
        measure_neg, centroid_neg = vofi_get_interface_centroid(neg_func3d, nothing, [0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 3)
        @test measure_neg == 0.0
        @test all(centroid_neg .== 0.0)
    end
    
    # Test integration over multiple cells - verify interface centroid is computed consistently
    @testset "Multi-cell integration consistency" begin
        # For a circle centered at (0,0), summing interface lengths * centroid / total length
        # should give approximately (0, 0) for the global interface centroid
        r = 0.35
        circle_sdf(x, _) = sqrt(x[1]^2 + x[2]^2) - r
        
        n = 10
        h = 1.0 / n
        xmin = -0.5
        
        total_length = 0.0
        weighted_centroid = [0.0, 0.0]
        
        for i in 0:n-1, j in 0:n-1
            xin = [xmin + i * h, xmin + j * h]
            measure, centroid = vofi_get_interface_centroid(circle_sdf, nothing, xin, [h, h], 2)
            if measure > 0
                total_length += measure
                weighted_centroid[1] += measure * centroid[1]
                weighted_centroid[2] += measure * centroid[2]
            end
        end
        
        # Expected perimeter: 2 * pi * r
        expected_perimeter = 2 * π * r
        @test total_length ≈ expected_perimeter atol=0.1
        
        # Global centroid should be approximately (0, 0)
        if total_length > 0
            global_centroid = weighted_centroid / total_length
            @test global_centroid[1] ≈ 0.0 atol=0.02
            @test global_centroid[2] ≈ 0.0 atol=0.02
        end
    end
end

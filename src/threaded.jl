# Threaded computation support for VofiJul
# 
# This module provides thread-safe parallel computation functions for
# integrating over multiple cells using Julia's multi-threading capabilities.

"""
    vofi_get_cc_threaded(impl_func, par, cells, h0, ndim0; 
                         nex=[0,0], npt=zeros(Int,4), nvis=[0,0]) -> Vector{Float64}

Compute volume fractions for multiple cells in parallel using threads.

# Arguments
- `impl_func`: Implicit function; returns `f(x, par)` (negative inside the reference phase).
- `par`: User data passed to `impl_func` (or `nothing`).
- `cells`: Vector of cell minimum corners (each element is a vector of coordinates).
- `h0`: Cell edge lengths (shared for all cells, or a vector of vectors for per-cell sizes).
- `ndim0`: 1, 2, 3, or 4 for the problem dimension.

# Keyword Arguments
- `nex`: Flags to compute centroid/interface; e.g. `[1,0]`. Default `[0,0]`.
- `npt`: User hints for quadrature points. Default `zeros(Int,4)`.
- `nvis`: Tecplot export flags. Default `[0,0]`.

# Returns
- Vector of volume fractions (cc values), one per cell.

# Example
```julia
using VofiJul

sphere_sdf(x, _) = sqrt(sum(abs2, x)) - 0.4

# Create grid of cells
n = 20
h = 1.0 / n
cells = [[x, y, z] for x in -0.5:h:0.5-h for y in -0.5:h:0.5-h for z in -0.5:h:0.5-h]

# Compute in parallel
cc_values = vofi_get_cc_threaded(sphere_sdf, nothing, cells, [h, h, h], 3)
total_volume = sum(cc_values) * h^3
```
"""
function vofi_get_cc_threaded(impl_func, par, cells::AbstractVector, h0, ndim0::Integer;
                               nex=[0,0], npt=zeros(Int,4), nvis=[0,0])
    ncells = length(cells)
    results = Vector{vofi_real}(undef, ncells)
    
    # Pre-allocate thread-local xex buffers - use max threads available
    max_threads = max(Threads.nthreads(), Threads.maxthreadid())
    xex_buffers = [zeros(vofi_real, 5) for _ in 1:max_threads]
    
    # Determine if h0 is shared or per-cell
    h0_is_shared = !isa(h0, AbstractVector{<:AbstractVector})
    
    Threads.@threads for i in 1:ncells
        tid = Threads.threadid()
        xex = xex_buffers[tid]
        fill!(xex, 0.0)
        
        cell_h0 = h0_is_shared ? h0 : h0[i]
        @inbounds results[i] = vofi_get_cc(impl_func, par, cells[i], cell_h0, xex, nex, npt, nvis, ndim0)
    end
    
    return results
end

"""
    vofi_integrate_threaded(impl_func, par, grid_min, grid_max, ncells, ndim0;
                            nex=[0,0], npt=zeros(Int,4), nvis=[0,0]) -> NamedTuple

Integrate over a regular Cartesian grid in parallel.

# Arguments
- `impl_func`: Implicit function; returns `f(x, par)` (negative inside the reference phase).
- `par`: User data passed to `impl_func` (or `nothing`).
- `grid_min`: Minimum corner of the grid domain.
- `grid_max`: Maximum corner of the grid domain.
- `ncells`: Number of cells in each dimension (Int or Vector{Int}).
- `ndim0`: 1, 2, 3, or 4 for the problem dimension.

# Keyword Arguments
- `nex`: Flags to compute centroid/interface; e.g. `[1,0]`. Default `[0,0]`.
- `npt`: User hints for quadrature points. Default `zeros(Int,4)`.
- `nvis`: Tecplot export flags. Default `[0,0]`.

# Returns
A NamedTuple with:
- `total`: Total integrated volume/area/length.
- `cc_values`: Matrix of volume fractions indexed by cell position.
- `cell_volume`: Volume of each cell.

# Example
```julia
using VofiJul

sphere_sdf(x, _) = sqrt(sum(abs2, x)) - 0.4

result = vofi_integrate_threaded(sphere_sdf, nothing, 
                                  [-0.5, -0.5, -0.5], [0.5, 0.5, 0.5],
                                  20, 3)
println("Computed volume: ", result.total)
println("Exact volume: ", 4/3 * π * 0.4^3)
```
"""
function vofi_integrate_threaded(impl_func, par, grid_min, grid_max, ncells, ndim0::Integer;
                                  nex=[0,0], npt=zeros(Int,4), nvis=[0,0])
    # Handle scalar ncells
    if isa(ncells, Integer)
        ncells_vec = fill(ncells, ndim0)
    else
        ncells_vec = collect(ncells)[1:ndim0]
    end
    
    # Compute cell sizes
    h = [(grid_max[d] - grid_min[d]) / ncells_vec[d] for d in 1:ndim0]
    
    # Compute cell volume
    cell_vol = prod(h)
    
    # Create array for results
    cc_array = Array{vofi_real}(undef, ncells_vec...)
    
    # Generate all cell indices
    total_cells = prod(ncells_vec)
    
    # Pre-allocate thread-local xex buffers
    max_threads = max(Threads.nthreads(), Threads.maxthreadid())
    xex_buffers = [zeros(vofi_real, 5) for _ in 1:max_threads]
    
    # Compute in parallel using linear indexing
    Threads.@threads for linear_idx in 1:total_cells
        tid = Threads.threadid()
        xex = xex_buffers[tid]
        fill!(xex, 0.0)
        
        # Convert linear index to Cartesian indices (0-based)
        idx = linear_idx - 1
        cell_indices = Vector{Int}(undef, ndim0)
        for d in 1:ndim0
            cell_indices[d] = idx % ncells_vec[d]
            idx = idx ÷ ncells_vec[d]
        end
        
        # Compute cell corner
        xin = [grid_min[d] + cell_indices[d] * h[d] for d in 1:ndim0]
        
        # Compute volume fraction
        cc = vofi_get_cc(impl_func, par, xin, h, xex, nex, npt, nvis, ndim0)
        
        # Store result using Cartesian indices (1-based)
        @inbounds cc_array[(cell_indices .+ 1)...] = cc
    end
    
    # Compute total
    total = sum(cc_array) * cell_vol
    
    return (total=total, cc_values=cc_array, cell_volume=cell_vol)
end

"""
    vofi_get_cell_type_threaded(impl_func, par, cells, h0, ndim0) -> Vector{Int}

Determine cell types for multiple cells in parallel.

Returns 1 if fully inside (f < 0), 0 if fully outside, -1 if cut by the interface.

# Arguments
- `impl_func`: Implicit function; returns `f(x, par)` (negative inside the reference phase).
- `par`: User data passed to `impl_func` (or `nothing`).
- `cells`: Vector of cell minimum corners.
- `h0`: Cell edge lengths (shared for all cells).
- `ndim0`: 1, 2, 3, or 4 for the problem dimension.

# Example
```julia
using VofiJul

sphere_sdf(x, _) = sqrt(sum(abs2, x)) - 0.4
cells = [[-0.5, -0.5, -0.5], [0.0, 0.0, 0.0], [0.3, 0.3, 0.3]]
types = vofi_get_cell_type_threaded(sphere_sdf, nothing, cells, [0.1, 0.1, 0.1], 3)
```
"""
function vofi_get_cell_type_threaded(impl_func, par, cells::AbstractVector, h0, ndim0::Integer)
    ncells = length(cells)
    results = Vector{Int}(undef, ncells)
    
    Threads.@threads for i in 1:ncells
        @inbounds results[i] = vofi_get_cell_type(impl_func, par, cells[i], h0, ndim0)
    end
    
    return results
end

"""
    vofi_interface_centroid_threaded(impl_func, par, cells, h0, ndim0; 
                                      tol=1e-10, max_iter=50) -> Vector{Union{Nothing, Vector{Float64}}}

Compute interface centroids for multiple cells in parallel.

Returns a vector where each element is either the centroid coordinates or `nothing` if no interface.

# Arguments
- `impl_func`: Implicit function; returns `f(x, par)` (negative inside the reference phase).
- `par`: User data passed to `impl_func` (or `nothing`).
- `cells`: Vector of cell minimum corners.
- `h0`: Cell edge lengths.
- `ndim0`: 1, 2, or 3 for the problem dimension.

# Keyword Arguments
- `tol`: Tolerance for the 1D root-finder (default `1e-10`).
- `max_iter`: Maximum iterations for the bisection root-finder (default `50`).
"""
function vofi_interface_centroid_threaded(impl_func, par, cells::AbstractVector, h0, ndim0::Integer;
                                           tol=1e-10, max_iter=50)
    ncells = length(cells)
    results = Vector{Union{Nothing, Vector{Float64}}}(undef, ncells)
    
    Threads.@threads for i in 1:ncells
        @inbounds results[i] = vofi_interface_centroid(impl_func, par, cells[i], h0, ndim0; 
                                                        tol=tol, max_iter=max_iter)
    end
    
    return results
end

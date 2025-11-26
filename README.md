# VofiJul.jl

![CI](https://github.com/PenguinxCutCell/VofiJul.jl/actions/workflows/ci.yml/badge.svg)
![Coverage](https://codecov.io/gh/PenguinxCutCell/VofiJul.jl/branch/main/graph/badge.svg)

Julia port of the [VOFI](https://github.com/vofi-dev/vofi) library for initializing volume fractions from an analytic implicit surface `f(x,y,z)=0`. Cells are treated as line segments (1D), rectangles (2D), or cuboids (3D), the reference phase is in the region `f < 0`, and integration follows the original VOFI algorithms.

## Usage (preview)

```julia
using VofiJul

# implicit function f(x) -> Float64 (negative inside)
sphere_sdf(x, _) = sqrt(sum(abs2, x)) - 0.4

# 3D example
xex = zeros(Float64, 4)                      # centroid / interface info
cc  = vofi_get_cc(sphere_sdf, nothing,
                  [-0.5, -0.5, -0.5],        # cell min corner
                  [1.0, 1.0, 1.0],           # cell sizes
                  xex, [1, 0], [0, 0, 0, 0], # centroids on
                  [0, 0], 3)                 # ndim=3

# 1D example
line_sdf(x, _) = x[1] - 0.25
xex_1d = zeros(Float64, 4)
cc_1d = vofi_get_cc(line_sdf, nothing,
                    [0.0],                   # cell min corner
                    [1.0],                   # cell size
                    xex_1d, [1, 0], [0, 0],  # centroid on
                    [0, 0], 1)               # ndim=1

# Interface centroid example
circle_sdf(x, _) = sqrt(x[1]^2 + x[2]^2) - 0.25
measure, centroid = vofi_get_interface_centroid(circle_sdf, nothing,
                                                 [0.0, 0.0],  # cell min corner
                                                 [0.5, 0.5],  # cell sizes
                                                 2)           # ndim=2
# measure: interface length (2D) or area (3D)
# centroid: position of interface centroid
```

See `test/runtests.jl` for more examples, including Cartesian integrations in 1D, 2D, and 3D.

## Main routines

`vofi_get_cc(impl_func, par, xin, h0, xex, nex, npt, nvis, ndim0)`
- `impl_func`: implicit function; returns `f(x, par)` (negative inside the reference phase).
- `par`: user data passed to `impl_func` (or `nothing`).
- `xin`: minimum corner of the cell.
- `h0`: cell edge lengths.
- `xex`: output buffer for centroid/interface data (length >= 4).
- `nex`: flags to compute centroid/interface (1D point, 2D length, or 3D area); e.g. `[1,0]`.
- `npt`: user hints for quadrature points (can be zeros to auto-select).
- `nvis`: Tecplot export flags (set to zeros to disable).
- `ndim0`: 1, 2, or 3 for the problem dimension.

`vofi_get_cell_type(impl_func, par, xin, h0, ndim0)`
- Same `impl_func`, `par`, `xin`, `h0`, `ndim0` as above.
- Returns `1` if the cell is fully inside (f < 0), `0` if fully outside, `-1` if cut.

`vofi_get_interface_centroid(impl_func, par, xin, h0, ndim0)`
- Same `impl_func`, `par`, `xin`, `h0`, `ndim0` as above.
- Returns a tuple `(interface_measure, interface_centroid)` where:
  - `interface_measure`: 0 for 1D (point), arc length for 2D, surface area for 3D
  - `interface_centroid`: position of the interface centroid (1D: interface point, 2D/3D: centroid coordinates)
- For cells that are fully inside (f < 0) or fully outside (f > 0), returns `(0.0, zeros(ndim0))`.


## Installation

The package is a pure Julia port; add it as a local dev package:

```shell
julia --project -e 'using Pkg; Pkg.develop(path="."); Pkg.test()'
```

## Credits

All credits for the original algorithms, design, and C implementation go to the VOFI authors:
Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, Ruben Scardovelli, Philip Yecko, and Stephane Zaleski. This Julia port is a reimplementation for the Julia ecosystem.

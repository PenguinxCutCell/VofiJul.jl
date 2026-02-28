# VofiJul.jl

![CI](https://github.com/PenguinxCutCell/VofiJul.jl/actions/workflows/ci.yml/badge.svg)
![Coverage](https://codecov.io/gh/PenguinxCutCell/VofiJul.jl/branch/main/graph/badge.svg)

Julia port of the [VOFI](https://github.com/vofi-dev/vofi) library for initializing volume fractions from an analytic implicit surface `f(x,y,z)=0`. Cells are treated as line segments (1D), rectangles (2D), cuboids (3D) or hyperectangles (4D), and the reference phase is in the region `f < 0`, and integration follows the original VOFI algorithms.

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
```

See `test/runtests.jl` for more examples, including Cartesian integrations in 1D, 2D, 3D, and 4D.

## Main routines

`vofi_get_cc(impl_func, par, xin, h0, xex, nex, npt, nvis, ndim0)`
- `impl_func`: implicit function; returns `f(x, par)` (negative inside the reference phase).
- `par`: user data passed to `impl_func` (or `nothing`).
- `xin`: minimum corner of the cell.
- `h0`: cell edge lengths.
- `xex`: output buffer for centroid/interface data (length ≥ 4).
- `nex`: flags to compute centroid/interface (1D point, 2D length, 3D area or 4D volume); e.g. `[1,0]`.
- `npt`: user hints for quadrature points (can be zeros to auto-select).
- `nvis`: Tecplot export flags (set to zeros to disable).
- `ndim0`: 1, 2, 3 or 4 for the problem dimension.

`vofi_get_cell_type(impl_func, par, xin, h0, ndim0)`
- Same `impl_func`, `par`, `xin`, `h0`, `ndim0` as above.
- Returns `1` if the cell is fully inside (f < 0), `0` if fully outside, `-1` if cut.

`vofi_interface_centroid(impl_func, par, xin, h0, ndim0; tol=1e-10, max_iter=50)`
- `impl_func`: implicit function; returns `f(x, par)` (negative inside the reference phase).
- `par`: user data passed to `impl_func` (or `nothing`).
- `xin`: minimum corner of the cell.
- `h0`: cell edge lengths.
- `ndim0`: 1, 2, 3 or 4 for the problem dimension.
- `tol`: tolerance for the 1D root-finder along the estimated normal (default `1e-10`).
- `max_iter`: maximum iterations for the bisection root-finder (default `50`).


## Installation

The package is a pure Julia port; add it as a local dev package:

```shell
julia --project -e 'using Pkg; Pkg.develop(path=\".\"); Pkg.test()'
```

## Current VOFI 

VOFI library work in 1D, 2D, and 3D, but the original C implementation is not designed for 4D. This Julia port extends the algorithms to support 4D hypercubes, which are relevant for certain applications (e.g. space-time cut cells).

## Credits

All credits for the original algorithms, design, and C implementation go to the VOFI authors:
Andrea Chierici, Leonardo Chirco, Vincent Le Chenadec, Ruben Scardovelli, Philip Yecko, and Stéphane Zaleski. This Julia port is a reimplementation for the Julia ecosystem. 

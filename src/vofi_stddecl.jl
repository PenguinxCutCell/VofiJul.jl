"""
Core constants, helper routines, and lightweight data structures that mirror
`vofi/include/vofi_stddecl.h` from the original C implementation.
"""
const vofi_real = Float64
const vofi_creal = Float64
const vofi_int = Int
const vofi_cint = Int
const vofi_void_cptr = Any
const vofi_int_cpt = Vector{vofi_int}
const Integrand = Function

const EPS_M = 1.5e-7
const EPS_LOC = 1.5e-7
const EPS_E = 5.0e-7
const EPS_SEGM = 1.0e-12
const EPS_ROOT = 1.0e-14
const EPS_NOT0 = 1.0e-90
const NEAR_EDGE_RATIO = 2.0e-2
const MAX_ITER_ROOT = 15
const MAX_ITER_MINI = 50
const NDIM = 3
const NVER = 4
const NSE = 2
const NSEG = 10
const NGLM = 20

# Type aliases for StaticArrays with common sizes
const SVec2 = SVector{NSE, vofi_real}
const SVec3 = SVector{NDIM, vofi_real}
const SVec4 = SVector{4, vofi_real}
const SVecNSEG = SVector{NSEG, vofi_real}
const SVecNSEGP1 = SVector{NSEG + 1, vofi_real}
const SVecNGLMP2 = SVector{NGLM + 2, vofi_real}
const SIVec2 = SVector{NSE, vofi_int}
const SIVec3 = SVector{NDIM, vofi_int}
const SIVecNSEG = SVector{NSEG, vofi_int}
const SMat2x2 = SMatrix{NSE, NSE, vofi_real, NSE * NSE}
const SMat3x3 = SMatrix{NDIM, NDIM, vofi_real, NDIM * NDIM}
const SMat2x3 = SMatrix{NSE, NDIM, vofi_int, NSE * NDIM}

@inline MIN(a, b) = min(a, b)
@inline MAX(a, b) = max(a, b)
@inline SGN0P(a) = a < 0 ? -1 : 1
@inline Sq(a) = a * a
@inline Sq2(a) = a[1] * a[1] + a[2] * a[2]
@inline Sq3(a) = a[1] * a[1] + a[2] * a[2] + a[3] * a[3]
@inline Sqd3(a, b) = Sq(a[1] - b[1]) + Sq(a[2] - b[2]) + Sq(a[3] - b[3])

macro SHFT4(a, b, c, d)
    quote
        $(esc(a)) = $(esc(b))
        $(esc(b)) = $(esc(c))
        $(esc(c)) = $(esc(d))
    end
end

macro CPSF(s, t, f, g)
    quote
        $(esc(s)) = $(esc(t))
        $(esc(f)) = $(esc(g))
    end
end

mutable struct MinData
    xval::MVector{NDIM, vofi_real}
    fval::vofi_real
    sval::vofi_real
    isc::MVector{NDIM, vofi_int}
    ipt::vofi_int
    function MinData(xval, fval, sval, isc, ipt)
        return new(xval, fval, sval, isc, ipt)
    end
end

function MinData(; xval = nothing,
                 fval = zero(vofi_real),
                 sval = zero(vofi_real),
                 isc = nothing,
                 ipt = zero(vofi_int))
    xv = xval === nothing ? @MVector(zeros(vofi_real, NDIM)) :
         MVector{NDIM, vofi_real}(xval)
    iscv = isc === nothing ? @MVector(zeros(vofi_int, NDIM)) :
            MVector{NDIM, vofi_int}(isc)
    return MinData(xv, fval, sval, iscv, ipt)
end

mutable struct MinData4D
    xval::Vector{vofi_real}
    fval::vofi_real
    sval::vofi_real
    span::vofi_real
    isc::Vector{vofi_int}
    function MinData4D(xval, fval, sval, span, isc)
        return new(xval, fval, sval, span, isc)
    end
end

function MinData4D(; xval = nothing,
                   fval = zero(vofi_real),
                   sval = zero(vofi_real),
                   span = zero(vofi_real),
                   isc = nothing)
    xv = xval === nothing ? zeros(vofi_real, 4) : copy(xval)
    iscv = isc === nothing ? zeros(vofi_int, 4) : copy(isc)
    return MinData4D(xv, fval, sval, span, iscv)
end

mutable struct XFSP4D
    edges::Vector{MinData4D}
    sectors::Vector{MinData4D}
    ipt::vofi_int
    function XFSP4D(edges, sectors, ipt)
        return new(edges, sectors, ipt)
    end
end

function XFSP4D()
    edges = [MinData4D() for _ in 1:8]
    sectors = [MinData4D() for _ in 1:2]
    return XFSP4D(edges, sectors, 0)
end

mutable struct DirData
    ind1::vofi_int
    ind2::vofi_int
    swt1::vofi_int
    swt2::vofi_int
    consi::vofi_int
    function DirData(ind1, ind2, swt1, swt2, consi)
        return new(ind1, ind2, swt1, swt2, consi)
    end
end

DirData(; ind1 = 0, ind2 = 0, swt1 = 0, swt2 = 0, consi = 0) =
    DirData(ind1, ind2, swt1, swt2, consi)

mutable struct LenData
    np0::vofi_int
    f_sign::vofi_int
    xt0::MVector{NGLM + 2, vofi_real}
    ht0::MVector{NGLM + 2, vofi_real}
    htp::MVector{NGLM + 2, vofi_real}
    function LenData(np0, f_sign, xt0, ht0, htp)
        length(xt0) == NGLM + 2 || throw(ArgumentError("xt0 must have NGLM + 2 entries"))
        length(ht0) == NGLM + 2 || throw(ArgumentError("ht0 must have NGLM + 2 entries"))
        length(htp) == NGLM + 2 || throw(ArgumentError("htp must have NGLM + 2 entries"))
        return new(np0, f_sign, xt0, ht0, htp)
    end
end

function LenData(; np0 = 0, f_sign = 1,
                 xt0 = nothing, ht0 = nothing, htp = nothing)
    xt = xt0 === nothing ? @MVector(zeros(vofi_real, NGLM + 2)) : MVector{NGLM + 2, vofi_real}(xt0)
    ht = ht0 === nothing ? @MVector(zeros(vofi_real, NGLM + 2)) : MVector{NGLM + 2, vofi_real}(ht0)
    hp = htp === nothing ? @MVector(zeros(vofi_real, NGLM + 2)) : MVector{NGLM + 2, vofi_real}(htp)
    return LenData(np0, f_sign, xt, ht, hp)
end

Base.copy!(dest::MinData, src::MinData) = begin
    dest.xval .= src.xval
    dest.fval = src.fval
    dest.sval = src.sval
    dest.isc .= src.isc
    dest.ipt = src.ipt
    dest
end

Base.copy!(dest::LenData, src::LenData) = begin
    dest.np0 = src.np0
    dest.f_sign = src.f_sign
    dest.xt0 .= src.xt0
    dest.ht0 .= src.ht0
    dest.htp .= src.htp
    dest
end

Base.copy!(dest::MinData4D, src::MinData4D) = begin
    dest.xval .= src.xval
    dest.fval = src.fval
    dest.sval = src.sval
    dest.span = src.span
    dest.isc .= src.isc
    dest
end

@inline function call_integrand(func::F, par, coords::AbstractVector) where {F}
    if par === nothing && applicable(func, coords)
        return func(coords)
    end
    return func(coords, par)
end

struct VofiCache
    x0::MVector{NDIM, vofi_real}
    hvec::MVector{NDIM, vofi_real}
    
    # 1D
    f0_1D::MVector{NSE, vofi_real}
    
    # 2D
    pdir2::MVector{NDIM, vofi_real}
    sdir2::MVector{NDIM, vofi_real}
    f02D::MMatrix{NSE, NSE, vofi_real, 4}
    xfsp_single::MinData
    base2::MVector{NSEG + 1, vofi_real}
    nsect2::MVector{NSEG, Int}
    ndire2::MVector{NSEG, Int}
    centroid2::MVector{NDIM + 1, vofi_real}
    xhp1::LenData
    xhp2::LenData
    basei_inner2D::MVector{NSEG + 1, Int}
    x1_inner2D::MVector{NDIM, vofi_real}
    x2_inner2D::MVector{NDIM, vofi_real}
    fse_inner2D::MVector{NSE, vofi_real}
    n02D::MMatrix{NSE, NSE, Int, 4}
    fc2D::MMatrix{NDIM, NDIM, vofi_real, 9}
    
    # 3D
    pdir3::MVector{NDIM, vofi_real}
    sdir3::MVector{NDIM, vofi_real}
    tdir3::MVector{NDIM, vofi_real}
    f03D::MArray{Tuple{NSE, NSE, NSE}, vofi_real, 3, 8}
    fc3D::MArray{Tuple{NDIM, NDIM, NDIM}, vofi_real, 3, 27}
    fd3D::MMatrix{NDIM, NDIM, vofi_real, 9}
    n03D::MArray{Tuple{NSE, NSE, NSE}, Int, 3, 8}
    x3D::MVector{NDIM, vofi_real}
    xfsp1::MinData
    xfsp2::MinData
    xfsp3::MinData
    xfsp4::MinData
    xfsp5::MinData
    xfsp_tuple::Tuple{MinData, MinData, MinData, MinData, MinData}
    base3::MVector{NSEG + 1, vofi_real}
    centroid3::MVector{NDIM + 1, vofi_real}
    xhpn1::LenData
    xhpn2::LenData
    xhpo1::LenData
    xhpo2::LenData
    xhpn_edge1::LenData
    xhpn_edge2::LenData
    xfs_vol::MinData
    nsect_vol::MVector{NSEG, Int}
    ndire_vol::MVector{NSEG, Int}
    xfsl_limits3D::MinData
    xfsl_check_plane::MinData
    xfsl_inner2D::MinData
    xfsl_edge2D::MinData
    xfsl_cb1::MinData
    xfsl_cb2::MinData
    xfsl_cb3::MinData
    xfsl_cb4::MinData
    xfsl_cb5::MinData
    
    # area
    x1_area::MVector{NDIM, vofi_real}
    x20_area::MVector{NDIM, vofi_real}
    x21_area::MVector{NDIM, vofi_real}
    s0_area::MVector{4, vofi_real}
    fse_area::MVector{NSE, vofi_real}
    
    # volume
    x1_vol::MVector{NDIM, vofi_real}
    base_int_vol::MVector{NSEG + 1, vofi_real}
    xmidt_vol::MVector{NGLM + 2, vofi_real}
    xedge_vol::MVector{NDIM, vofi_real}
    
    # 4D
    x0_4::Vector{vofi_real}
    h4::Vector{vofi_real}
    pdir4::Vector{vofi_real}
    sdir4::Vector{vofi_real}
    tdir4::Vector{vofi_real}
    qdir4::Vector{vofi_real}
    f04D::Array{vofi_real, 4}
    xfsp4D::XFSP4D
    base4::Vector{vofi_real}
    centroid4::Vector{vofi_real}
    hvec4::Vector{vofi_real}
    hh4::Vector{vofi_real}
    x4::Vector{vofi_real}
    fgrad4::Vector{vofi_real}
    n04D::Array{Int, 4}
    mags4::Vector{vofi_real}
    order4::Vector{Int}
end

function VofiCache()
    VofiCache(
        @MVector(zeros(vofi_real, NDIM)),
        @MVector(zeros(vofi_real, NDIM)),
        
        @MVector(zeros(vofi_real, NSE)),
        
        @MVector(zeros(vofi_real, NDIM)),
        @MVector(zeros(vofi_real, NDIM)),
        @MMatrix(zeros(vofi_real, NSE, NSE)),
        MinData(),
        @MVector(zeros(vofi_real, NSEG + 1)),
        @MVector(zeros(Int, NSEG)),
        @MVector(zeros(Int, NSEG)),
        @MVector(zeros(vofi_real, NDIM + 1)),
        LenData(),
        LenData(),
        @MVector(zeros(Int, NSEG + 1)),
        @MVector(zeros(vofi_real, NDIM)),
        @MVector(zeros(vofi_real, NDIM)),
        @MVector(zeros(vofi_real, NSE)),
        @MMatrix(zeros(Int, NSE, NSE)),
        @MMatrix(zeros(vofi_real, NDIM, NDIM)),
        
        @MVector(zeros(vofi_real, NDIM)),
        @MVector(zeros(vofi_real, NDIM)),
        @MVector(zeros(vofi_real, NDIM)),
        @MArray(zeros(vofi_real, NSE, NSE, NSE)),
        @MArray(zeros(vofi_real, NDIM, NDIM, NDIM)),
        @MMatrix(zeros(vofi_real, NDIM, NDIM)),
        @MArray(zeros(Int, NSE, NSE, NSE)),
        @MVector(zeros(vofi_real, NDIM)),
        MinData(), MinData(), MinData(), MinData(), MinData(),
        (MinData(), MinData(), MinData(), MinData(), MinData()),
        @MVector(zeros(vofi_real, NSEG + 1)),
        @MVector(zeros(vofi_real, NDIM + 1)),
        LenData(), LenData(), LenData(), LenData(), LenData(), LenData(),
        MinData(),
        @MVector(zeros(Int, NSEG)),
        @MVector(zeros(Int, NSEG)),
        MinData(), MinData(), MinData(), MinData(),
        MinData(), MinData(), MinData(), MinData(), MinData(),
        
        @MVector(zeros(vofi_real, NDIM)),
        @MVector(zeros(vofi_real, NDIM)),
        @MVector(zeros(vofi_real, NDIM)),
        @MVector(zeros(vofi_real, 4)),
        @MVector(zeros(vofi_real, NSE)),
        
        @MVector(zeros(vofi_real, NDIM)),
        @MVector(zeros(vofi_real, NSEG + 1)),
        @MVector(zeros(vofi_real, NGLM + 2)),
        @MVector(zeros(vofi_real, NDIM)),
        
        zeros(vofi_real, 4),
        zeros(vofi_real, 4),
        zeros(vofi_real, 4),
        zeros(vofi_real, 4),
        zeros(vofi_real, 4),
        zeros(vofi_real, 4),
        zeros(vofi_real, NSE, NSE, NSE, NSE),
        XFSP4D(),
        zeros(vofi_real, NSEG + 1),
        zeros(vofi_real, 5),
        zeros(vofi_real, 4),
        zeros(vofi_real, 4),
        zeros(vofi_real, 4),
        zeros(vofi_real, 4),
        zeros(Int, NSE, NSE, NSE, NSE),
        zeros(vofi_real, 4),
        zeros(Int, 4)
    )
end

function get_vofi_cache()
    get!(task_local_storage(), :vofi_cache) do
        VofiCache()
    end::VofiCache
end

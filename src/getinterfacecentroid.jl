"""
    vofi_get_interface_centroid(impl_func, par, xin, h0, ndim0)

Compute interface centroid and interface measure (length in 2D, area in 3D, volume in 4D) for a cell.

# Arguments
- `impl_func`: implicit function; returns `f(x, par)` (negative inside the reference phase)
- `par`: user data passed to `impl_func` (or `nothing`)
- `xin`: minimum corner of the cell
- `h0`: cell edge lengths
- `ndim0`: 1, 2, 3, or 4 for the problem dimension

# Returns
A tuple `(interface_measure, interface_centroid)` where:
- `interface_measure`: 0 for 1D (point), arc length for 2D, surface area for 3D, hypervolume for 4D
- `interface_centroid`: position of the interface centroid (1D: interface point, 2D/3D/4D: centroid coordinates)

For cells that are fully inside (f < 0) or fully outside (f > 0), returns `(0.0, zeros(ndim0))`.

# Example
```julia
# 2D circle
circle_sdf(x, _) = sqrt(x[1]^2 + x[2]^2) - 0.25
measure, centroid = vofi_get_interface_centroid(circle_sdf, nothing, [-0.5, -0.5], [1.0, 1.0], 2)
```
"""
function vofi_get_interface_centroid(impl_func, par, xin, h0, ndim0)
    x0 = @MVector zeros(vofi_real, NDIM)
    hvec = pad_to_ndim(h0)
    
    if ndim0 == 1
        x0[1] = xin[1]
        x0[2] = 0.0
        x0[3] = 0.0
        f0_1D = @MVector zeros(vofi_real, NSE)
        icc = vofi_order_dirs_1D(impl_func, par, x0, hvec, f0_1D)
        
        if icc >= 0
            # Cell is fully inside or fully outside - no interface
            return (0.0, zeros(vofi_real, 1))
        end
        
        # Find interface position (zero crossing)
        interface_pos = vofi_get_interface_position_1D(impl_func, par, x0, hvec, f0_1D)
        return (0.0, [interface_pos])  # 0 measure for point interface
        
    elseif ndim0 == 2
        x0[1] = xin[1]
        x0[2] = xin[2]
        xfsp_single = MinData()
        pdir = @MVector zeros(vofi_real, NDIM)
        sdir = @MVector zeros(vofi_real, NDIM)
        f02D = @MMatrix zeros(vofi_real, NSE, NSE)
        
        icc = vofi_order_dirs_2D(impl_func, par, x0, hvec, pdir, sdir, f02D, xfsp_single)
        
        if icc >= 0
            # Cell is fully inside or fully outside - no interface
            return (0.0, zeros(vofi_real, 2))
        end
        
        base = @MVector zeros(vofi_real, NSEG + 1)
        nsect = @MVector zeros(Int, NSEG)
        ndire = @MVector zeros(Int, NSEG)
        nsub = vofi_get_limits_2D(impl_func, par, x0, hvec, f02D, xfsp_single, base,
                                  pdir, sdir, nsect, ndire)
        centroid = @MVector zeros(vofi_real, NDIM + 1)
        xhp1 = LenData()
        xhp2 = LenData()
        npt = zeros(Int, 4)
        
        # Compute area (this populates xhp1, xhp2 with interface data)
        vofi_get_area(impl_func, par, x0, hvec, base, pdir, sdir, (xhp1, xhp2),
                     centroid, 0, npt, nsub, xfsp_single.ipt, nsect, ndire)
        
        # Compute interface length and centroid
        arc_length, interface_centroid = vofi_interface_length_with_centroid(
            impl_func, par, x0, hvec, pdir, sdir, (xhp1, xhp2), 0, true)
        
        return (arc_length, interface_centroid)
        
    elseif ndim0 == 3
        x0[1] = xin[1]
        x0[2] = xin[2]
        x0[3] = xin[3]
        pdir = @MVector zeros(vofi_real, NDIM)
        sdir = @MVector zeros(vofi_real, NDIM)
        tdir = @MVector zeros(vofi_real, NDIM)
        f03D = @MArray zeros(vofi_real, NSE, NSE, NSE)
        
        xfsp1 = MinData()
        xfsp2 = MinData()
        xfsp3 = MinData()
        xfsp4 = MinData()
        xfsp5 = MinData()
        xfsp = (xfsp1, xfsp2, xfsp3, xfsp4, xfsp5)
        
        icc = vofi_order_dirs_3D(impl_func, par, x0, hvec, pdir, sdir, tdir, f03D, xfsp)
        
        if icc >= 0
            # Cell is fully inside or fully outside - no interface
            return (0.0, zeros(vofi_real, 3))
        end
        
        base_ext = @MVector zeros(vofi_real, NSEG + 1)
        nsub_ext = vofi_get_limits_3D(impl_func, par, x0, hvec, f03D, xfsp, base_ext, pdir, sdir, tdir)
        
        # Compute interface area and centroid using integration
        surface_area, interface_centroid = vofi_get_interface_area_centroid_3D(
            impl_func, par, x0, hvec, base_ext, pdir, sdir, tdir, nsub_ext, xfsp[5].ipt)
        
        return (surface_area, interface_centroid)
        
    elseif ndim0 == 4
        # Use dynamic arrays for 4D (NDIM is 3; do not change it)
        x0_4 = Vector{vofi_real}(undef, 4)
        for i in 1:4
            x0_4[i] = xin[i]
        end
        length(h0) >= 4 || throw(ArgumentError("h0 must provide 4 entries when ndim0 == 4"))
        h4 = Vector{vofi_real}(undef, 4)
        for i in 1:4
            h4[i] = vofi_real(h0[i])
        end
        pdir = zeros(vofi_real, 4)
        sdir = zeros(vofi_real, 4)
        tdir = zeros(vofi_real, 4)
        qdir = zeros(vofi_real, 4)
        f04D = zeros(vofi_real, NSE, NSE, NSE, NSE)
        xfsp = XFSP4D()

        icc = vofi_order_dirs_4D(impl_func, par, x0_4, h4, pdir, sdir, tdir, qdir, f04D, xfsp)
        if icc >= 0
            # Cell is fully inside or fully outside - no interface
            return (0.0, zeros(vofi_real, 4))
        end

        base = zeros(vofi_real, NSEG + 1)
        nsub = vofi_get_limits_4D(impl_func, par, x0_4, h4, f04D, xfsp, base, pdir, sdir, tdir, qdir)
        
        # Compute interface volume and centroid using integration
        interface_volume, interface_centroid = vofi_get_interface_volume_centroid_4D(
            impl_func, par, x0_4, h4, base, pdir, sdir, tdir, qdir, nsub, xfsp.ipt)
        
        return (interface_volume, interface_centroid)
        
    else
        throw(ArgumentError("ndim0 must be 1, 2, 3, or 4 for interface centroid computation"))
    end
end

"""
    vofi_get_interface_position_1D(impl_func, par, x0, h0, f0)

Find the interface position (zero crossing) in 1D.
"""
function vofi_get_interface_position_1D(impl_func, par, x0, h0, f0)
    # f0 has two values: f0[1] at x0[1] and f0[2] at x0[1] + h0[1]
    if f0[1] * f0[2] >= 0
        # No zero crossing
        return x0[1] + 0.5 * h0[1]
    end
    
    # Find the zero crossing using linear interpolation as initial guess
    denom = abs(f0[1]) + abs(f0[2])
    if denom < EPS_NOT0
        frac = 0.5
    else
        frac = abs(f0[1]) / denom
    end
    x_zero = x0[1] + frac * h0[1]
    
    # Refine with Newton-Raphson
    x1 = @MVector zeros(vofi_real, NDIM)
    x1[2] = x0[2]
    x1[3] = x0[3]
    
    for _ in 1:MAX_ITER_ROOT
        x1[1] = x_zero
        f_val = call_integrand(impl_func, par, x1)
        
        if abs(f_val) < EPS_ROOT
            break
        end
        
        h_eps = EPS_LOC * h0[1]
        x1[1] = x_zero + h_eps
        f_plus = call_integrand(impl_func, par, x1)
        f_deriv = (f_plus - f_val) / h_eps
        
        if abs(f_deriv) < EPS_NOT0
            break
        end
        
        delta = -f_val / f_deriv
        x_zero += delta
        x_zero = max(x0[1], min(x0[1] + h0[1], x_zero))
        
        if abs(delta) < EPS_ROOT * h0[1]
            break
        end
    end
    
    return x_zero
end

"""
    vofi_get_interface_area_centroid_3D(impl_func, par, x0, h0, base_ext, pdir, sdir, tdir, nsub_ext, nptmp)

Compute interface area and centroid for 3D case.
"""
function vofi_get_interface_area_centroid_3D(impl_func, par, x0, h0, base_ext, pdir, sdir, tdir, nsub_ext, nptmp)
    x1 = @MVector zeros(vofi_real, NDIM)
    base_int = @MVector zeros(vofi_real, NSEG + 1)
    xmidt = @MVector zeros(vofi_real, NGLM + 2)
    
    hp = hs = ht = zero(vofi_real)
    for i in 1:NDIM
        hp += pdir[i] * h0[i]
        hs += sdir[i] * h0[i]
        ht += tdir[i] * h0[i]
    end
    hm = maximum(h0)
    
    xhpn1 = LenData()
    xhpn2 = LenData()
    xhpo1 = LenData()
    xhpo2 = LenData()
    xfs = MinData()
    nsect = @MVector zeros(Int, NSEG)
    ndire = @MVector zeros(Int, NSEG)
    
    total_area = 0.0
    centroid_acc = zeros(vofi_real, NDIM)
    npt = zeros(Int, 4)
    centroid_tmp = @MVector zeros(vofi_real, NDIM + 1)
    
    for nt in 1:nsub_ext
        dt = base_ext[nt + 1] - base_ext[nt]
        mdpt = 0.5 * (base_ext[nt + 1] + base_ext[nt])
        for i in 1:NDIM
            x1[i] = x0[i] + tdir[i] * mdpt
            xfs.isc[i] = 0
        end
        
        sect_hexa = vofi_check_plane(impl_func, par, x1, h0, xfs, base_int, pdir, sdir, nsect, ndire)
        if sect_hexa == 0
            continue
        end
        
        dt_scaled = 18 * dt / hm
        if !isfinite(dt_scaled) || dt_scaled < 0
            dt_scaled = 0.0
        end
        nexpt = min(20, Int(floor(dt_scaled)) + 3)
        ptx_ext = gauss_legendre_nodes(nexpt)
        ptw_ext = gauss_legendre_weights(nexpt)
        
        xhpo1 = LenData()
        xhpo2 = LenData()
        xhpn_edge1 = LenData()
        xhpn_edge2 = LenData()
        xmidt[1] = base_ext[nt]
        xmidt[nexpt + 2] = base_ext[nt + 1]
        
        for k in 1:nexpt
            xit = mdpt + 0.5 * dt * ptx_ext[k]
            xmidt[k + 1] = xit
            for i in 1:NDIM
                x1[i] = x0[i] + tdir[i] * xit
            end
            
            nsub_int = vofi_get_limits_inner_2D(impl_func, par, x1, h0, xfs, base_int, pdir, sdir, nsect, ndire, sect_hexa)
            xhpn1 = LenData()
            xhpn2 = LenData()
            vofi_get_area(impl_func, par, x1, h0, base_int, pdir, sdir, (xhpn1, xhpn2),
                         centroid_tmp, 0, npt, nsub_int, nptmp, nsect, ndire)
            
            vofi_end_points(impl_func, par, x1, h0, pdir, sdir, (xhpn1, xhpn2))
            
            if k == 1
                xedge = @MVector zeros(vofi_real, NDIM)
                for i in 1:NDIM
                    xedge[i] = x0[i] + tdir[i] * xmidt[1]
                end
                nintmp = vofi_get_limits_edge_2D(impl_func, par, xedge, h0, xfs, base_int, pdir, sdir, nsub_int)
                vofi_edge_points(impl_func, par, xedge, h0, base_int, pdir, sdir, (xhpo1, xhpo2), (xhpn1.np0, xhpn2.np0), nintmp, nsect, ndire)
                vofi_end_points(impl_func, par, xedge, h0, pdir, sdir, (xhpo1, xhpo2))
            elseif k > 1 && k < nexpt
                area_contrib, centroid_contrib = vofi_interface_surface_with_centroid(
                    impl_func, par, x0, h0, xmidt, pdir, sdir, tdir, (xhpn1, xhpn2), (xhpo1, xhpo2), k, nexpt, 0, true)
                total_area += area_contrib
                for i in 1:NDIM
                    centroid_acc[i] += area_contrib * centroid_contrib[i]
                end
                copy!(xhpo1, xhpn1)
                copy!(xhpo2, xhpn2)
            else
                xedge = @MVector zeros(vofi_real, NDIM)
                for i in 1:NDIM
                    xedge[i] = x0[i] + tdir[i] * xmidt[nexpt + 2]
                end
                nintmp = vofi_get_limits_edge_2D(impl_func, par, xedge, h0, xfs, base_int, pdir, sdir, nsub_int)
                vofi_edge_points(impl_func, par, xedge, h0, base_int, pdir, sdir, (xhpn_edge1, xhpn_edge2), (xhpo1.np0, xhpo2.np0), nintmp, nsect, ndire)
                vofi_end_points(impl_func, par, xedge, h0, pdir, sdir, (xhpn_edge1, xhpn_edge2))
                area_contrib, centroid_contrib = vofi_interface_surface_with_centroid(
                    impl_func, par, x0, h0, xmidt, pdir, sdir, tdir, (xhpn_edge1, xhpn_edge2), (xhpo1, xhpo2), k + 1, nexpt, 0, true)
                total_area += area_contrib
                for i in 1:NDIM
                    centroid_acc[i] += area_contrib * centroid_contrib[i]
                end
            end
        end
    end
    
    # Normalize centroid
    interface_centroid = zeros(vofi_real, NDIM)
    if total_area > 0
        for i in 1:NDIM
            interface_centroid[i] = centroid_acc[i] / total_area
        end
    end
    
    return total_area, interface_centroid
end

"""
    vofi_get_interface_volume_centroid_4D(impl_func, par, x0, h0, base, pdir, sdir, tdir, qdir, nsub, nptmp)

Compute interface volume and centroid for 4D case.
The interface in 4D is a 3D hypersurface, so the "interface measure" is a 3D volume.
"""
function vofi_get_interface_volume_centroid_4D(impl_func, par, x0, h0, base, pdir, sdir, tdir, qdir, nsub, nptmp)
    ax_p = axis_index(pdir)
    ax_s = axis_index(sdir)
    ax_t = axis_index(tdir)
    ax_q = axis_index(qdir)
    hp = axis_length(pdir, h0)
    hs = axis_length(sdir, h0)
    ht = axis_length(tdir, h0)
    hq = axis_length(qdir, h0)
    hm = maximum(h0)

    # Build a 3D slice coordinate system
    xin3 = Vector{vofi_real}(undef, NDIM)
    xin3[1] = x0[ax_p]
    xin3[2] = x0[ax_s]
    xin3[3] = x0[ax_t]
    h3 = Vector{vofi_real}(undef, NDIM)
    h3[1] = hp
    h3[2] = hs
    h3[3] = ht

    xbuf = similar(x0)
    q_current = Ref(x0[ax_q])
    
    # Create a 3D slice function that evaluates the 4D function at a fixed q
    slice_func = let impl_func = impl_func, par = par, x0 = x0,
                     ax_p = ax_p, ax_s = ax_s, ax_t = ax_t, ax_q = ax_q,
                     xbuf = xbuf, q_current = q_current
        function (coords)
            for i in 1:length(x0)
                xbuf[i] = x0[i]
            end
            xbuf[ax_p] = coords[1]
            xbuf[ax_s] = coords[2]
            xbuf[ax_t] = coords[3]
            xbuf[ax_q] = q_current[]
            return call_integrand(impl_func, par, xbuf)
        end
    end

    total_volume = 0.0
    centroid_acc = zeros(vofi_real, 4)
    q_origin = x0[ax_q]
    max_nodes = NGLM

    for ns in 1:nsub
        dq = base[ns + 1] - base[ns]
        if dq <= EPS_LOC
            continue
        end
        mdpt = 0.5 * (base[ns + 1] + base[ns])
        nquad = clamp(Int(floor(18 * dq / hm)) + 3, 3, 20)
        if nptmp > 0
            nquad = min(nquad, nptmp)
        end
        nquad = min(nquad, max_nodes)
        nodes = gauss_legendre_nodes(nquad)
        weights = gauss_legendre_weights(nquad)
        
        seg_vol = 0.0
        seg_xp = seg_xs = seg_xt = seg_xq = 0.0
        
        for k in 1:nquad
            xi = mdpt + 0.5 * dq * nodes[k]
            xi = clamp(xi, 0.0, hq)
            q_abs = q_origin + xi
            q_current[] = q_abs
            
            # Get interface area and centroid for the 3D slice
            surface_area, centroid3 = vofi_get_interface_centroid(slice_func, nothing, xin3, h3, 3)
            
            w = weights[k]
            seg_vol += w * surface_area
            
            if surface_area > 0
                # Accumulate weighted centroid contributions
                seg_xp += w * surface_area * centroid3[1]
                seg_xs += w * surface_area * centroid3[2]
                seg_xt += w * surface_area * centroid3[3]
                seg_xq += w * surface_area * q_abs
            end
        end
        
        factor = 0.5 * dq
        total_volume += factor * seg_vol
        centroid_acc[1] += factor * seg_xp
        centroid_acc[2] += factor * seg_xs
        centroid_acc[3] += factor * seg_xt
        centroid_acc[4] += factor * seg_xq
    end

    # Normalize centroid and map back to original coordinate indices
    interface_centroid = zeros(vofi_real, 4)
    if total_volume > 0
        centroid_acc ./= total_volume
        interface_centroid[ax_p] = centroid_acc[1]
        interface_centroid[ax_s] = centroid_acc[2]
        interface_centroid[ax_t] = centroid_acc[3]
        interface_centroid[ax_q] = centroid_acc[4]
    end

    return total_volume, interface_centroid
end

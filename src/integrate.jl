function vofi_get_length_1D(impl_func, par, x0, h0, f0, xex, ncen)
    # For 1D, we just need to find the zero crossing
    # f0 has two values: f0[1] at x0[1] and f0[2] at x0[1] + h0[1]
    
    # If both have same sign, something went wrong - should have been caught earlier
    if f0[1] * f0[2] >= 0
        # No zero crossing - return full length if negative, 0 otherwise
        if f0[1] < 0 && f0[2] < 0
            if ncen > 0
                xex[1] = x0[1] + 0.5 * h0[1]
            end
            return h0[1]
        else
            return 0.0
        end
    end
    
    # Find the zero crossing using linear interpolation as initial guess
    # then refine with Newton-Raphson
    denom = abs(f0[1]) + abs(f0[2])
    if denom < EPS_NOT0
        # Both values very close to zero - use midpoint
        frac = 0.5
    else
        frac = abs(f0[1]) / denom
    end
    x_zero = x0[1] + frac * h0[1]
    
    # Refine the zero crossing location
    x1 = @MVector zeros(vofi_real, NDIM)
    x1[2] = x0[2]
    x1[3] = x0[3]
    
    # Newton-Raphson iteration
    for iter in 1:MAX_ITER_ROOT
        x1[1] = x_zero
        f_val = call_integrand(impl_func, par, x1)
        
        if abs(f_val) < EPS_ROOT
            break
        end
        
        # Compute numerical derivative
        h_eps = EPS_LOC * h0[1]
        x1[1] = x_zero + h_eps
        f_plus = call_integrand(impl_func, par, x1)
        f_deriv = (f_plus - f_val) / h_eps
        
        if abs(f_deriv) < EPS_NOT0
            break
        end
        
        # Newton step
        delta = -f_val / f_deriv
        x_zero += delta
        
        # Keep within bounds
        x_zero = max(x0[1], min(x0[1] + h0[1], x_zero))
        
        if abs(delta) < EPS_ROOT * h0[1]
            break
        end
    end
    
    # The length inside (negative region) is from x0[1] to x_zero if f0[1] < 0
    # or from x_zero to x0[1] + h0[1] if f0[2] < 0
    if f0[1] < 0
        length_inside = x_zero - x0[1]
        if ncen > 0
            # Centroid is at the midpoint of the negative region
            xex[1] = x0[1] + 0.5 * length_inside
        end
    else
        length_inside = (x0[1] + h0[1]) - x_zero
        if ncen > 0
            # Centroid is at the midpoint of the negative region
            xex[1] = x_zero + 0.5 * length_inside
        end
    end
    
    return length_inside
end

function vofi_get_area(impl_func, par, x0, h0, base, pdir, sdir, xhp, centroid, ncen, npt, nsub, nptmp, nsect, ndire)
    x1 = @MVector zeros(vofi_real, NDIM)
    x20 = @MVector zeros(vofi_real, NDIM)
    x21 = @MVector zeros(vofi_real, NDIM)
    s0 = @MVector zeros(vofi_real, 4)
    fse = @MVector zeros(vofi_real, NSE)
    area = 0.0
    hp = 0.0
    hs = 0.0
    for i in 1:NDIM
        x1[i] = x0[i] + pdir[i] * h0[i]
        hp += pdir[i] * h0[i]
        hs += sdir[i] * h0[i]
    end
    hm = maximum(h0)
    xp = 0.0
    xs = 0.0
    it0 = 1
    for ns in 1:nsub
        ds = base[ns + 1] - base[ns]
        mdpt = 0.5 * (base[ns + 1] + base[ns])
        if nsect[ns] > 0
            al = ds * hp
            area += al
            if ncen > 0
                xp += 0.5 * hp * al
                xs += mdpt * al
            end
        elseif nsect[ns] < 0
            npts = Int(clamp(floor(18 * ds / hm) + 3, 3, 20))
            npts = min(nptmp, npts)
            if 3 <= npt[2] <= 20
                npts = min(npt[2], npts)
            end
            if 3 <= npt[1] <= 20
                npts = max(npt[1], npts)
            end
            xhp[it0].np0 = npts
            f_sign = ndire[ns]
            xhp[it0].f_sign = f_sign
            j = npts - 3
            pts = gauss_legendre_nodes(j + 3)
            wts = gauss_legendre_weights(j + 3)

            quada = 0.0
            quadp = 0.0
            quads = 0.0
            a1 = a2 = b1 = b2 = 0.0
            xhp[it0].ht0[1] = xhp[it0].htp[1] = 0.0
            xhp[it0].xt0[1] = base[ns]
            xhp[it0].xt0[npts + 2] = base[ns + 1]
            xhp[it0].xt0[2] = mdpt + 0.5 * ds * pts[1]
            for i in 1:NDIM
                tmp = sdir[i] * xhp[it0].xt0[2]
                x20[i] = x0[i] + tmp
                x21[i] = x1[i] + tmp
            end
            fse[1] = call_integrand(impl_func, par, x20)
            fse[2] = call_integrand(impl_func, par, x21)
            s0[1] = hp
            if abs(fse[1]) < abs(fse[2])
                s0[2] = 0.0
                s0[3] = fse[1]
            else
                s0[2] = hp
                s0[3] = fse[2]
            end
            s0[4] = (fse[2] - fse[1]) / hp
            for k in 1:npts
                xhp[it0].ht0[k + 1] = vofi_get_segment_zero(impl_func, par, x20, pdir, s0, f_sign)
                xhp[it0].htp[k + 1] = s0[4]
                quada += wts[k] * xhp[it0].ht0[k + 1]
                if ncen > 0
                    quadp += wts[k] * 0.5 * xhp[it0].ht0[k + 1]^2
                    quads += wts[k] * xhp[it0].ht0[k + 1] * xhp[it0].xt0[k + 1]
                end
                if k < npts
                    xhp[it0].xt0[k + 2] = mdpt + 0.5 * ds * pts[k + 1]
                    s0[2] = xhp[it0].ht0[k + 1]
                    s0[4] = xhp[it0].htp[k + 1]
                    if k > 1
                        dxm1 = xhp[it0].xt0[k + 1] - xhp[it0].xt0[k]
                        dxp1 = xhp[it0].xt0[k + 2] - xhp[it0].xt0[k + 1]
                        a1 = (xhp[it0].ht0[k + 1] - xhp[it0].ht0[k]) / dxm1
                        s0[2] += a1 * dxp1
                        b1 = (xhp[it0].htp[k + 1] - xhp[it0].htp[k]) / dxm1
                        s0[4] += b1 * dxp1
                        if k > 2
                            dxm2 = xhp[it0].xt0[k + 1] - xhp[it0].xt0[k - 1]
                            dxp2 = xhp[it0].xt0[k + 2] - xhp[it0].xt0[k]
                            s0[2] += (a1 - a2) * dxp1 * dxp2 / dxm2
                            s0[4] += (b1 - b2) * dxp1 * dxp2 / dxm2
                        end
                    end
                    if f_sign < 0
                        s0[2] = hp - s0[2]
                    end
                    ratio = s0[2] / hp
                    if ratio < NEAR_EDGE_RATIO
                        s0[2] = 0.0
                    elseif ratio > 1 - NEAR_EDGE_RATIO
                        s0[2] = hp
                    end
                    for i in 1:NDIM
                        x20[i] = x0[i] + sdir[i] * xhp[it0].xt0[k + 2]
                        x21[i] = x20[i] + pdir[i] * s0[2]
                    end
                    s0[3] = call_integrand(impl_func, par, x21)
                end
                a2 = a1
                b2 = b1
            end
            quada *= 0.5 * ds
            area += quada
            if ncen > 0 && quada > 0
                quadp = 0.5 * ds * quadp / quada
                quads = 0.5 * ds * quads / quada
                if f_sign < 0
                    quadp = hp - quadp
                end
                xp += quadp * quada
                xs += quads * quada
            end
            it0 += 1
        end
    end
    centroid[1] = xp
    centroid[2] = xs
    return area
end

function vofi_get_volume(impl_func, par, x0, h0, base_ext, pdir, sdir, tdir,
                         centroid, nex, npt, nsub_ext, nptmp, nvis)
    x1 = @MVector zeros(vofi_real, NDIM)
    base_int = @MVector zeros(vofi_real, NSEG + 1)
    xmidt = @MVector zeros(vofi_real, NGLM + 2)
    volume = 0.0
    surfer = 0.0
    # Accumulators for interface centroid (weighted by surface area, in local coords)
    iface_cent_t = 0.0
    iface_cent_s = 0.0
    iface_cent_p = 0.0
    hp = hs = ht = 0.0
    for i in 1:NDIM
        hp += pdir[i] * h0[i]
        hs += sdir[i] * h0[i]
        ht += tdir[i] * h0[i]
    end
    hm = maximum(h0)
    xp = xs = xt = 0.0
    xhpn1 = LenData()
    xhpn2 = LenData()
    xhpo1 = LenData()
    xhpo2 = LenData()
    xfs = MinData()
    nsect = @MVector zeros(Int, NSEG)
    ndire = @MVector zeros(Int, NSEG)

    for nt in 1:nsub_ext
        dt = base_ext[nt + 1] - base_ext[nt]
        mdpt = 0.5 * (base_ext[nt + 1] + base_ext[nt])
        for i in 1:NDIM
            x1[i] = x0[i] + tdir[i] * mdpt
            xfs.isc[i] = 0
        end
        sect_hexa = vofi_check_plane(impl_func, par, x1, h0, xfs, base_int, pdir, sdir,
                                     nsect, ndire)
        if sect_hexa == 0
            if nsect[1] == 1
                vol = dt * hs * hp
                volume += vol
                if nex[1] > 0
                    xp += 0.5 * hp * vol
                    xs += 0.5 * hs * vol
                    xt += mdpt * vol
                end
            end
            continue
        end

        dt_scaled = 18 * dt / hm
        if !isfinite(dt_scaled) || dt_scaled < 0
            dt_scaled = 0.0
        end
        nexpt = min(20, Int(floor(dt_scaled)) + 3)
        if 3 <= npt[4] <= 20
            nexpt = min(npt[4], nexpt)
        end
        if 3 <= npt[3] <= 20
            nexpt = max(npt[3], nexpt)
        end
        ptx_ext = gauss_legendre_nodes(nexpt)
        ptw_ext = gauss_legendre_weights(nexpt)

        quadv = quadp = quads = quadt = 0.0
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
            nsub_int = vofi_get_limits_inner_2D(impl_func, par, x1, h0, xfs, base_int,
                                                pdir, sdir, nsect, ndire, sect_hexa)
            # Reset xhpn structures for this iteration
            xhpn1 = LenData()
            xhpn2 = LenData()
            area = vofi_get_area(impl_func, par, x1, h0, base_int, pdir, sdir, (xhpn1, xhpn2),
                                 centroid, nex[1], npt, nsub_int, nptmp, nsect, ndire)
            if nvis[1] > 0
                tecplot_heights(x1, h0, pdir, sdir, (xhpn1, xhpn2))
            end
            if nex[2] > 0
                vofi_end_points(impl_func, par, x1, h0, pdir, sdir, (xhpn1, xhpn2))
                if k == 1
                    xedge = @MVector zeros(vofi_real, NDIM)
                    for i in 1:NDIM
                        xedge[i] = x0[i] + tdir[i] * xmidt[1]
                    end
                    nintmp = vofi_get_limits_edge_2D(impl_func, par, xedge, h0, xfs,
                                                     base_int, pdir, sdir, nsub_int)
                    vofi_edge_points(impl_func, par, xedge, h0, base_int, pdir, sdir,
                                     (xhpo1, xhpo2), (xhpn1.np0, xhpn2.np0), nintmp, nsect, ndire)
                    vofi_end_points(impl_func, par, xedge, h0, pdir, sdir, (xhpo1, xhpo2))
                elseif k > 1 && k < nexpt
                    surf_contrib, (ct, cs, cp) = vofi_interface_surface_and_centroid(impl_func, par, x0, h0, xmidt, pdir,
                                                     sdir, tdir, (xhpn1, xhpn2), (xhpo1, xhpo2), k, nexpt, nvis[2])
                    surfer += surf_contrib
                    iface_cent_t += ct
                    iface_cent_s += cs
                    iface_cent_p += cp
                    copy!(xhpo1, xhpn1)
                    copy!(xhpo2, xhpn2)
                else
                    xedge = @MVector zeros(vofi_real, NDIM)
                    for i in 1:NDIM
                        xedge[i] = x0[i] + tdir[i] * xmidt[nexpt + 2]
                    end
                    nintmp = vofi_get_limits_edge_2D(impl_func, par, xedge, h0, xfs,
                                                     base_int, pdir, sdir, nsub_int)
                    vofi_edge_points(impl_func, par, xedge, h0, base_int, pdir, sdir,
                                     (xhpn_edge1, xhpn_edge2), (xhpo1.np0, xhpo2.np0), nintmp, nsect, ndire)
                    vofi_end_points(impl_func, par, xedge, h0, pdir, sdir, (xhpn_edge1, xhpn_edge2))
                    surf_contrib, (ct, cs, cp) = vofi_interface_surface_and_centroid(impl_func, par, x0, h0, xmidt, pdir,
                                                     sdir, tdir, (xhpn_edge1, xhpn_edge2), (xhpo1, xhpo2), k + 1, nexpt, nvis[2])
                    surfer += surf_contrib
                    iface_cent_t += ct
                    iface_cent_s += cs
                    iface_cent_p += cp
                end
            end
            quadv += ptw_ext[k] * area
            quadp += ptw_ext[k] * centroid[1]
            quads += ptw_ext[k] * centroid[2]
            quadt += ptw_ext[k] * area * xit
        end
        quadv *= 0.5 * dt
        volume += quadv
        if nex[1] > 0
            xp += 0.5 * dt * quadp
            xs += 0.5 * dt * quads
            xt += 0.5 * dt * quadt
        end
    end

    centroid[1] = xp
    centroid[2] = xs
    centroid[3] = xt
    centroid[4] = surfer
    # Store normalized interface centroid in elements 5-7 (local coords: t, s, p)
    if length(centroid) >= 7 && surfer > EPS_NOT0
        centroid[5] = iface_cent_t / surfer  # t-coordinate
        centroid[6] = iface_cent_s / surfer  # s-coordinate
        centroid[7] = iface_cent_p / surfer  # p-coordinate (height)
    elseif length(centroid) >= 7
        centroid[5] = 0.0
        centroid[6] = 0.0
        centroid[7] = 0.0
    end
    return volume
end

function vofi_get_hypervolume(impl_func, par, x0, h0, base, pdir, sdir, tdir, qdir,
                              centroid, nex, npt, nsub, nptmp, nvis)
    ax_p = axis_index(pdir)
    ax_s = axis_index(sdir)
    ax_t = axis_index(tdir)
    ax_q = axis_index(qdir)
    hp = axis_length(pdir, h0)
    hs = axis_length(sdir, h0)
    ht = axis_length(tdir, h0)
    hq = axis_length(qdir, h0)
    prod3 = hp * hs * ht
    hm = maximum(h0)

    xin3 = Vector{vofi_real}(undef, NDIM)
    xin3[1] = x0[ax_p]
    xin3[2] = x0[ax_s]
    xin3[3] = x0[ax_t]
    h3 = Vector{vofi_real}(undef, NDIM)
    h3[1] = hp
    h3[2] = hs
    h3[3] = ht
    # Extend xex3 to receive interface centroid: vol_centroid(3) + surface(1) + iface_centroid(3)
    xex3 = zeros(vofi_real, 7)
    nex_slice = zeros(Int, 2)
    want_centroid = nex[1] > 0
    want_surface = (length(nex) >= 2) && nex[2] > 0
    want_iface_centroid = want_surface && nex[2] > 1
    if want_centroid
        nex_slice[1] = 1
    end
    if want_iface_centroid
        nex_slice[2] = 2  # Request interface centroid from 3D slices
    elseif want_surface
        nex_slice[2] = 1
    end
    nvis_slice = zeros(Int, 2)

    xbuf = similar(x0)
    q_current = Ref(x0[ax_q])
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

    hypervolume = 0.0
    xp_acc = xs_acc = xt_acc = xq_acc = 0.0
    surface_acc = 0.0
    # Accumulators for interface centroid (weighted by surface area)
    iface_xp_acc = iface_xs_acc = iface_xt_acc = iface_xq_acc = 0.0
    q_origin = x0[ax_q]
    max_nodes = NGLM

    for ns in 1:nsub
        dq = base[ns + 1] - base[ns]
        if dq <= EPS_LOC
            continue
        end
        mdpt = 0.5 * (base[ns + 1] + base[ns])
        nquad = clamp(Int(floor(18 * dq / hm)) + 3, 3, 20)
        if length(npt) >= 4 && 3 <= npt[4] <= 20
            nquad = min(nquad, npt[4])
        end
        if nptmp > 0
            nquad = min(nquad, nptmp)
        end
        nquad = min(nquad, max_nodes)
        nodes = gauss_legendre_nodes(nquad)
        weights = gauss_legendre_weights(nquad)
        seg_vol = seg_xp = seg_xs = seg_xt = seg_xq = seg_surface = 0.0
        seg_iface_xp = seg_iface_xs = seg_iface_xt = seg_iface_xq = 0.0
        for k in 1:nquad
            xi = mdpt + 0.5 * dq * nodes[k]
            xi = clamp(xi, 0.0, hq)
            q_abs = q_origin + xi
            q_current[] = q_abs
            fill!(xex3, 0.0)
            cc = vofi_get_cc(slice_func, nothing, xin3, h3, xex3, nex_slice, npt, nvis_slice, 3)
            vol3 = cc * prod3
            w = weights[k]
            seg_vol += w * vol3
            if want_centroid && vol3 > 0
                seg_xp += w * vol3 * xex3[1]
                seg_xs += w * vol3 * xex3[2]
                seg_xt += w * vol3 * xex3[3]
                seg_xq += w * vol3 * q_abs
            end
            if want_surface && nex_slice[2] > 0
                slice_surf = xex3[4]
                seg_surface += w * slice_surf
                # Accumulate interface centroid weighted by surface area
                if want_iface_centroid && slice_surf > 0
                    seg_iface_xp += w * slice_surf * xex3[5]
                    seg_iface_xs += w * slice_surf * xex3[6]
                    seg_iface_xt += w * slice_surf * xex3[7]
                    seg_iface_xq += w * slice_surf * q_abs
                end
            end
        end
        factor = 0.5 * dq
        hypervolume += factor * seg_vol
        if want_centroid
            xp_acc += factor * seg_xp
            xs_acc += factor * seg_xs
            xt_acc += factor * seg_xt
            xq_acc += factor * seg_xq
        end
        if want_surface
            surface_acc += factor * seg_surface
        end
        if want_iface_centroid
            iface_xp_acc += factor * seg_iface_xp
            iface_xs_acc += factor * seg_iface_xs
            iface_xt_acc += factor * seg_iface_xt
            iface_xq_acc += factor * seg_iface_xq
        end
    end

    centroid[1] = xp_acc
    centroid[2] = xs_acc
    centroid[3] = xt_acc
    centroid[4] = xq_acc
    centroid[5] = surface_acc
    # Store normalized interface centroid in elements 6-9 if requested
    if length(centroid) >= 9 && surface_acc > EPS_NOT0
        centroid[6] = iface_xp_acc / surface_acc
        centroid[7] = iface_xs_acc / surface_acc
        centroid[8] = iface_xt_acc / surface_acc
        centroid[9] = iface_xq_acc / surface_acc
    elseif length(centroid) >= 9
        centroid[6] = 0.0
        centroid[7] = 0.0
        centroid[8] = 0.0
        centroid[9] = 0.0
    end
    return hypervolume
end

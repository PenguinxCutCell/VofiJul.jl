function vofi_get_cc(impl_func::F, par, xin, h0, xex, nex, npt, nvis, ndim0)::vofi_real where {F}
    # Ensure we have enough slots for barycenter coords plus interface measure
    nex_interface = length(nex) >= 2 ? nex[2] : 0
    required_len = (ndim0 == 4 && nex_interface > 0) ? 5 : max(4, ndim0)
    if length(xex) < required_len
        resize!(xex, required_len)
    end
    fill!(xex, 0)
    
    cache = get_vofi_cache()
    x0 = cache.x0
    fill!(x0, 0)
    hvec = pad_to_ndim!(cache.hvec, h0)
    
    if ndim0 == 1
        x0[1] = xin[1]
        x0[2] = 0.0
        x0[3] = 0.0
        f0_1D = cache.f0_1D
        fill!(f0_1D, 0)
        icc = vofi_order_dirs_1D(impl_func, par, x0, hvec, f0_1D)
        if icc >= 0
            cc = vofi_real(icc)
            if icc > 0 && nex[1] > 0
                xex[1] = x0[1] + 0.5 * hvec[1]
            end
            return cc
        end
        # Interface crosses the cell - compute the fraction
        length_frac = vofi_get_length_1D(impl_func, par, x0, hvec, f0_1D, xex, nex[1])
        cc = length_frac / hvec[1]
        return cc
    elseif ndim0 == 2
        x0[1] = xin[1]
        x0[2] = xin[2]
        xfsp_single = cache.xfsp_single
        pdir = cache.pdir2
        sdir = cache.sdir2
        fill!(pdir, 0)
        fill!(sdir, 0)
        f02D = cache.f02D
        fill!(f02D, 0)
        icc = vofi_order_dirs_2D(impl_func, par, x0, hvec, pdir, sdir, f02D, xfsp_single)
        if icc >= 0
            cc = vofi_real(icc)
            if icc > 0 && nex[1] > 0
                for i in 1:NSE
                    xex[i] = x0[i] + 0.5 * hvec[i]
                end
            end
            return cc
        end
        base = cache.base2
        fill!(base, 0)
        nsect = cache.nsect2
        fill!(nsect, 0)
        ndire = cache.ndire2
        fill!(ndire, 0)
        nsub = vofi_get_limits_2D(impl_func, par, x0, hvec, f02D, xfsp_single, base,
                                  pdir, sdir, nsect, ndire)
        centroid = cache.centroid2
        fill!(centroid, 0)
        xhp1 = cache.xhp1
        xhp2 = cache.xhp2
        area = vofi_get_area(impl_func, par, x0, hvec, base, pdir, sdir, (xhp1, xhp2),
                             centroid, nex[1], npt, nsub, xfsp_single.ipt, nsect, ndire)
        cc = area / (hvec[1] * hvec[2])
        if nvis[1] > 0
            tecplot_heights(x0, hvec, pdir, sdir, (xhp1, xhp2))
        end
        if nex[1] > 0 && area > 0
            centroid[1] /= area
            centroid[2] /= area
            centroid[3] = 0.0
            for i in 1:2
                xex[i] = x0[i] + centroid[1] * pdir[i] + centroid[2] * sdir[i]
            end
        end
        if nex[2] > 0
            xex[4] = vofi_interface_length(impl_func, par, x0, hvec, pdir, sdir, (xhp1, xhp2), nvis[2])
        end
        return cc
    elseif ndim0 == 3
        x0[1] = xin[1]
        x0[2] = xin[2]
        x0[3] = xin[3]
        pdir = cache.pdir3
        sdir = cache.sdir3
        tdir = cache.tdir3
        fill!(pdir, 0)
        fill!(sdir, 0)
        fill!(tdir, 0)
        f03D = cache.f03D
        fill!(f03D, 0)
        xfsp1 = cache.xfsp1
        xfsp2 = cache.xfsp2
        xfsp3 = cache.xfsp3
        xfsp4 = cache.xfsp4
        xfsp5 = cache.xfsp5
        xfsp = cache.xfsp_tuple

        icc = vofi_order_dirs_3D(impl_func, par, x0, hvec, pdir, sdir, tdir, f03D, xfsp)
        if icc >= 0
            cc = vofi_real(icc)
            if icc > 0 && nex[1] > 0
                for i in 1:NDIM
                    xex[i] = x0[i] + 0.5 * hvec[i]
                end
            end
            return cc
        end

        base = cache.base3
        fill!(base, 0)
        nsub = vofi_get_limits_3D(impl_func, par, x0, hvec, f03D, xfsp, base, pdir, sdir, tdir)
        centroid = cache.centroid3
        fill!(centroid, 0)
        volume = vofi_get_volume(impl_func, par, x0, hvec, base, pdir, sdir, tdir, centroid,
                                 nex, npt, nsub, xfsp[5].ipt, nvis)
        cc = volume / (hvec[1] * hvec[2] * hvec[3])
        if nex[1] > 0 && volume > 0
            centroid[1] /= volume
            centroid[2] /= volume
            centroid[3] /= volume
            for i in 1:NDIM
                xex[i] = x0[i] + centroid[1] * pdir[i] + centroid[2] * sdir[i] + centroid[3] * tdir[i]
            end
        end
        if nex[2] > 0
            xex[4] = centroid[4]
        end
        return cc
    elseif ndim0 == 4
        # Use dynamic arrays for 4D (NDIM is 3; do not change it)
        x0_4 = cache.x0_4
        for i in 1:4
            x0_4[i] = xin[i]
        end
        length(h0) >= 4 || throw(ArgumentError("h0 must provide 4 entries when ndim0 == 4"))
        h4 = cache.h4
        for i in 1:4
            h4[i] = vofi_real(h0[i])
        end
        pdir = cache.pdir4
        sdir = cache.sdir4
        tdir = cache.tdir4
        qdir = cache.qdir4
        fill!(pdir, 0)
        fill!(sdir, 0)
        fill!(tdir, 0)
        fill!(qdir, 0)
        f04D = cache.f04D
        fill!(f04D, 0)
        xfsp = cache.xfsp4D

        icc = vofi_order_dirs_4D(impl_func, par, x0_4, h4, pdir, sdir, tdir, qdir, f04D, xfsp)
        if icc >= 0
            cc = vofi_real(icc)
            if icc > 0 && nex[1] > 0
                for i in 1:4
                    xex[i] = x0_4[i] + 0.5 * h4[i]
                end
            end
            return cc
        end

        base = cache.base4
        fill!(base, 0)
        nsub = vofi_get_limits_4D(impl_func, par, x0_4, h4, f04D, xfsp, base, pdir, sdir, tdir, qdir)
        centroid = cache.centroid4
        fill!(centroid, 0)
        hypervolume = vofi_get_hypervolume(impl_func, par, x0_4, h4, base,
                          pdir, sdir, tdir, qdir,
                          centroid, nex, npt, nsub,
                          xfsp.ipt, nvis)
        cc = hypervolume / (h4[1] * h4[2] * h4[3] * h4[4])
        if nex[1] > 0 && hypervolume > 0
            for k in 1:4
                centroid[k] /= hypervolume
            end
            # centroid[1..4] are absolute coordinates in permuted (p,s,t,q) space
            # Map them back to original coordinate indices
            ax_p = axis_index(pdir)
            ax_s = axis_index(sdir)
            ax_t = axis_index(tdir)
            ax_q = axis_index(qdir)
            xex[ax_p] = centroid[1]
            xex[ax_s] = centroid[2]
            xex[ax_t] = centroid[3]
            xex[ax_q] = centroid[4]
        end
        if nex[2] > 0
            xex[5] = centroid[5]
        end
        return cc
    else
        throw(ArgumentError("ndim0 must be 1, 2, 3, or 4"))
    end
end

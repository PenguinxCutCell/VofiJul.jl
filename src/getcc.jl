function vofi_get_cc(impl_func, par, xin, h0, xex, nex, npt, nvis, ndim0)
    xex .= 0
    x0 = @MVector zeros(vofi_real, NDIM)
    hvec = pad_to_ndim(h0)
    if ndim0 == 1
        x0[1] = xin[1]
        x0[2] = 0.0
        x0[3] = 0.0
        f0_1D = @MVector zeros(vofi_real, NSE)
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
        xfsp_single = MinData()  # Only need one for 2D
        pdir = @MVector zeros(vofi_real, NDIM)
        sdir = @MVector zeros(vofi_real, NDIM)
        f02D = @MMatrix zeros(vofi_real, NSE, NSE)
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
        base = @MVector zeros(vofi_real, NSEG + 1)
        nsect = @MVector zeros(Int, NSEG)
        ndire = @MVector zeros(Int, NSEG)
        nsub = vofi_get_limits_2D(impl_func, par, x0, hvec, f02D, xfsp_single, base,
                                  pdir, sdir, nsect, ndire)
        centroid = @MVector zeros(vofi_real, NDIM + 1)
        xhp1 = LenData()
        xhp2 = LenData()
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
        pdir = @MVector zeros(vofi_real, NDIM)
        sdir = @MVector zeros(vofi_real, NDIM)
        tdir = @MVector zeros(vofi_real, NDIM)
        f03D = @MArray zeros(vofi_real, NSE, NSE, NSE)
        # Pre-allocate 5 MinData structs without array allocation
        xfsp1 = MinData()
        xfsp2 = MinData()
        xfsp3 = MinData()
        xfsp4 = MinData()
        xfsp5 = MinData()
        xfsp = (xfsp1, xfsp2, xfsp3, xfsp4, xfsp5)

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

        base = @MVector zeros(vofi_real, NSEG + 1)
        nsub = vofi_get_limits_3D(impl_func, par, x0, hvec, f03D, xfsp, base, pdir, sdir, tdir)
        centroid = @MVector zeros(vofi_real, NDIM + 1)
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
            cc = vofi_real(icc)
            if icc > 0 && nex[1] > 0
                for i in 1:4
                    xex[i] = x0_4[i] + 0.5 * h4[i]
                end
            end
            return cc
        end

        base = zeros(vofi_real, NSEG + 1)
        nsub = vofi_get_limits_4D(impl_func, par, x0_4, h4, f04D, xfsp, base, pdir, sdir, tdir, qdir)
        centroid = zeros(vofi_real, 5)
        hypervolume = vofi_get_hypervolume(impl_func, par, x0_4, h4, base,
                          pdir, sdir, tdir, qdir,
                          centroid, nex, npt, nsub,
                          xfsp.ipt, nvis)
        cc = hypervolume / (h4[1] * h4[2] * h4[3] * h4[4])
        if nex[1] > 0 && hypervolume > 0
            for k in 1:4
                centroid[k] /= hypervolume
            end
            for i in 1:4
                xex[i] = x0_4[i] + centroid[1] * pdir[i] + centroid[2] * sdir[i] +
                                 centroid[3] * tdir[i] + centroid[4] * qdir[i]
            end
        end
        if nex[2] > 0
            xex[5] = centroid[5]
        end
        return cc
    else
        throw(ArgumentError("ndim0 must be 1, 2, 3, or 4"))
    end
end

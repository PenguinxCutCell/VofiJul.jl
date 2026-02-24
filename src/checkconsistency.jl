function vofi_check_side_consistency(impl_func, par, x0, dir, fse, h0)
    fs = fse[1] + fse[2]
    consi = fs > 0 ? 1 : fs < 0 ? -1 : 0
    if consi != 0
        dh = max(EPS_M * h0, EPS_ROOT)
        f0 = abs(fse[1])
        f1 = abs(fse[2])
        ft = f0
        if f0 > f1
            ft = f1
            dh = h0 - dh
        end
        xs = similar(x0)
        for i in 1:NDIM
            xs[i] = x0[i] + dh * dir[i]
        end
        fs_val = consi * call_integrand(impl_func, par, xs)
        if fs_val >= ft
            consi = 0
        end
    end
    return consi
end

function vofi_check_face_consistency(impl_func, par, x0, h0, dir1, dir2, fv)
    ipsc = DirData(0, 0, 0, 0, 0)
    h1 = 0.0
    h2 = 0.0
    for i in 1:NDIM
        h1 += dir1[i] * h0[i]
        h2 += dir2[i] * h0[i]
    end
    sumf = fv[1] + fv[2] + fv[3] + fv[4]
    ipsc.consi = sumf > 0 ? 1 : sumf < 0 ? -1 : 0
    if ipsc.consi != 0
        dh1 = max(EPS_M * h1, EPS_ROOT)
        dh2 = max(EPS_M * h2, EPS_ROOT)
        f0 = abs(sumf)
        fl = abs.(fv)
        is1 = 1
        is2 = 1
        if fl[1] < f0
            f0 = fl[1]
            is1 = 1
            is2 = 1
        end
        if fl[2] < f0
            f0 = fl[2]
            ipsc.ind1 = 1
            is1 = -1
            is2 = 1
        end
        if fl[3] < f0
            f0 = fl[3]
            ipsc.ind1 = 0
            ipsc.ind2 = 1
            is1 = 1
            is2 = -1
        end
        if fl[4] < f0
            f0 = fl[4]
            ipsc.ind1 = 1
            ipsc.ind2 = 1
            is1 = -1
            is2 = -1
        end

        xx = @MVector zeros(vofi_real, NDIM)
        x1 = similar(xx)
        x2 = similar(xx)
        for i in 1:NDIM
            xx[i] = x0[i] + h1 * ipsc.ind1 * dir1[i] + h2 * ipsc.ind2 * dir2[i]
            x1[i] = xx[i] + dh1 * is1 * dir1[i]
            x2[i] = xx[i] + dh2 * is2 * dir2[i]
        end
        consi = 0
        f1 = ipsc.consi * call_integrand(impl_func, par, x1)
        if f1 < f0
            consi = ipsc.consi
            ipsc.swt1 = 1
        end
        f2 = ipsc.consi * call_integrand(impl_func, par, x2)
        if f2 < f0
            consi = ipsc.consi
            ipsc.swt2 = 1
        end
        ipsc.consi = consi
    end
    return ipsc
end

function vofi_check_line_consistency(impl_func, par, x0, dir, h0, n, xfs::MinData)
    consi = 0
    dh = max(EPS_M * h0, EPS_ROOT)
    xs = similar(x0)
    for i in 1:NDIM
        xs[i] = x0[i] + (1 - 2 * n) * h0 * dir[i]
    end
    f1 = call_integrand(impl_func, par, xs)
    for i in 1:NDIM
        xs[i] = x0[i] + (1 - 2 * n) * dh * dir[i]
    end
    f0 = call_integrand(impl_func, par, xs)
    if f0 * f1 <= 0
        consi = 1
        xfs.xval .= xs
        xfs.fval = f0
        xfs.sval = dh
        xfs.isc[1] = 1
        xfs.isc[2] = 1
    end
    return consi
end

function vofi_check_edge_consistency(impl_func, par, fse, x0, base, dir, h0, nsub)
    xs = similar(x0)
    s0 = @MVector zeros(vofi_real, 4)
    dh = max(EPS_M * h0, EPS_ROOT)
    if abs(fse[1]) < abs(fse[2])
        for i in 1:NDIM
            xs[i] = x0[i] + dh * dir[i]
        end
        fse[1] = call_integrand(impl_func, par, xs)
        if fse[1] * fse[2] > 0
            base[nsub + 1] = 0.0
        else
            f2neg = fse[1] < 0 ? 1 : -1
            s0[1] = h0 - dh
            s0[2] = 0.0
            s0[3] = fse[1]
            s0[4] = (fse[2] - fse[1]) / s0[1]
            dhl = vofi_get_segment_zero(impl_func, par, xs, dir, s0, f2neg)
            if f2neg < 0
                dhl = s0[1] - dhl
            end
            base[nsub + 1] = dhl + dh
        end
    else
        for i in 1:NDIM
            xs[i] = x0[i] + (h0 - dh) * dir[i]
        end
        fse[2] = call_integrand(impl_func, par, xs)
        if fse[1] * fse[2] > 0
            base[nsub + 1] = h0
        else
            f2neg = fse[1] < 0 ? 1 : -1
            s0[1] = h0 - dh
            s0[2] = s0[1]
            s0[3] = fse[2]
            s0[4] = (fse[2] - fse[1]) / s0[1]
            dhl = vofi_get_segment_zero(impl_func, par, xs, dir, s0, f2neg)
            if f2neg < 0
                dhl = s0[1] - dhl
            end
            base[nsub + 1] = dhl
        end
    end
    return nothing
end

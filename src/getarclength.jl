function vofi_interface_length(impl_func, par, x0, h0, pdir, sdir, xhhp, ipf)
    arc, _ = vofi_interface_length_and_centroid(impl_func, par, x0, h0, pdir, sdir, xhhp, ipf)
    return arc
end

"""
    vofi_interface_length_and_centroid(impl_func, par, x0, h0, pdir, sdir, xhhp, ipf)

Compute both the interface arc length and the interface centroid using precomputed heights.
Returns `(arc_length, (cent_p, cent_s))` where the centroid tuple is in local (pdir, sdir) coordinates.

This is more robust than bisection-based methods as it uses the already computed height
function values from the Gauss-Legendre integration.
"""
function vofi_interface_length_and_centroid(impl_func, par, x0, h0, pdir, sdir, xhhp, ipf)
    hp = zero(vofi_real)
    for i in 1:NDIM
        hp += pdir[i] * h0[i]
    end
    arc = 0.0
    # Accumulators for centroid computation (in local pdir, sdir coordinates)
    cent_p = 0.0  # weighted sum of p-coordinate (height direction)
    cent_s = 0.0  # weighted sum of s-coordinate (secondary direction)
    s0 = @MVector zeros(vofi_real, 4)
    s0[1] = hp
    nseg = xhhp[min(end, 2)].np0 > 0 ? 2 : (xhhp[1].np0 > 0 ? 1 : 0)
    x20 = @MVector zeros(vofi_real, NDIM)
    x21 = @MVector zeros(vofi_real, NDIM)
    sqrt3 = sqrt(3.0)

    for it0 in 1:nseg
        seg = xhhp[it0]
        npt = seg.np0
        npt > 0 || continue
        f_sign = seg.f_sign
        xt0 = seg.xt0
        ht0 = seg.ht0
        htp = seg.htp
        idx(i) = i + 1

        j0, j1, j2, j3 = 0, 1, 2, 3
        for _ in 1:NSE
            dx1 = xt0[idx(j1)] - xt0[idx(j2)]
            dx2 = xt0[idx(j2)] - xt0[idx(j3)]
            dx12 = xt0[idx(j1)] - xt0[idx(j3)]
            dc1 = xt0[idx(j0)] - xt0[idx(j1)]
            dc2 = xt0[idx(j0)] - xt0[idx(j2)]
            a1 = (ht0[idx(j1)] - ht0[idx(j2)]) / dx1
            a2 = (ht0[idx(j2)] - ht0[idx(j3)]) / dx2
            b1 = (htp[idx(j1)] - htp[idx(j2)]) / dx1
            b2 = (htp[idx(j2)] - htp[idx(j3)]) / dx2
            s0[2] = ht0[idx(j1)] + a1 * dc1 + (a1 - a2) * dc1 * dc2 / dx12
            s0[4] = htp[idx(j1)] + b1 * dc1 + (b1 - b2) * dc1 * dc2 / dx12
            if f_sign < 0
                s0[2] = hp - s0[2]
            end
            ratio = s0[2] / hp
            if ratio < NEAR_EDGE_RATIO
                s0[2] = 0.0
            elseif ratio > 1.0 - NEAR_EDGE_RATIO
                s0[2] = hp
            end
            for i in 1:NDIM
                x20[i] = x0[i] + sdir[i] * xt0[idx(j0)]
                x21[i] = x20[i] + pdir[i] * s0[2]
            end
            s0[3] = call_integrand(impl_func, par, x21)
            ht0[idx(j0)] = vofi_get_segment_zero(impl_func, par, x20, pdir, s0, f_sign)
            htp[idx(j0)] = s0[4]
            j0 = npt + 1
            j1 = npt
            j2 = npt - 1
            j3 = npt - 2
        end

        x0b = xt0[idx(0)]
        h0b = ht0[idx(0)]
        hpb = htp[idx(0)]
        xm = 0.5 * (x0b + xt0[idx(1)])
        xl = x0b
        hl = h0b
        xr = xt0[idx(1)]
        hr = ht0[idx(1)]
        k = 0

        for j in 0:npt
            dx1 = x0b - xt0[idx(k + 1)]
            dx2 = xt0[idx(k + 1)] - xt0[idx(k + 2)]
            dx12 = x0b - xt0[idx(k + 2)]
            dc1 = xm - x0b
            dc2 = xm - xt0[idx(k + 1)]
            a1 = (h0b - ht0[idx(k + 1)]) / dx1
            b1 = (hpb - htp[idx(k + 1)]) / dx1
            a2 = (ht0[idx(k + 1)] - ht0[idx(k + 2)]) / dx2
            b2 = (htp[idx(k + 1)] - htp[idx(k + 2)]) / dx2
            s0[2] = h0b + a1 * dc1 + (a1 - a2) * dc1 * dc2 / dx12
            s0[4] = hpb + b1 * dc1 + (b1 - b2) * dc1 * dc2 / dx12
            if f_sign < 0
                s0[2] = hp - s0[2]
            end
            ratio = s0[2] / hp
            if ratio < NEAR_EDGE_RATIO
                s0[2] = 0.0
            elseif ratio > 1.0 - NEAR_EDGE_RATIO
                s0[2] = hp
            end
            for i in 1:NDIM
                x20[i] = x0[i] + sdir[i] * xm
                x21[i] = x20[i] + pdir[i] * s0[2]
            end
            s0[3] = call_integrand(impl_func, par, x21)
            hm = vofi_get_segment_zero(impl_func, par, x20, pdir, s0, f_sign)
            hpm = s0[4]
            hsum = hl + hr
            hc = 0.5 * hsum + sqrt3 / 3.0 * (2 * hm - hsum)
            xc = xm
            d1 = sqrt((xl - xc)^2 + (hl - hc)^2)
            d2 = sqrt((xr - xc)^2 + (hr - hc)^2)
            arc += d1 + d2
            
            # Compute centroid contribution: midpoint of each segment weighted by segment length
            # Segment 1: from (hl, xl) to (hc, xc)
            mid_p1 = 0.5 * (hl + hc)
            mid_s1 = 0.5 * (xl + xc)
            cent_p += d1 * mid_p1
            cent_s += d1 * mid_s1
            
            # Segment 2: from (hc, xc) to (hr, xr)
            mid_p2 = 0.5 * (hc + hr)
            mid_s2 = 0.5 * (xc + xr)
            cent_p += d2 * mid_p2
            cent_s += d2 * mid_s2

            x0b = xm
            h0b = hm
            hpb = hpm
            k = min(j, npt - 1)
            xl = xt0[idx(k + 1)]
            hl = ht0[idx(k + 1)]
            xr = xt0[idx(k + 2)]
            hr = ht0[idx(k + 2)]
            xm = 0.5 * (xl + xr)
        end
    end

    # Normalize the centroid by total arc length
    if arc > EPS_NOT0
        cent_p /= arc
        cent_s /= arc
    else
        # Fallback: no interface or very small interface
        cent_p = 0.0
        cent_s = 0.0
    end
    
    return arc, (cent_p, cent_s)
end

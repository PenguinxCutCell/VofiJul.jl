function vofi_get_segment_zero(impl_func::F, par, x0, dir, s0, f_sign) where {F}
    xs = @MVector zeros(vofi_real, NDIM)
    sl = 0.0
    sr = s0[1]
    not_conv = true
    dsold = s0[1]
    dss = dsold
    ss = s0[2]
    fs = f_sign * s0[3]
    fps = f_sign * s0[4]
    sv = [ss, ss, ss]
    fv = [fs, fps, 0.0]
    gensec = 0
    fl = -EPS_SEGM
    fr = EPS_SEGM

    if fs < 0.0
        sl = ss
        fl = fs
    elseif fs > 0.0
        sr = ss
        fr = fs
    else
        not_conv = false
    end

    iter = 0
    while not_conv && iter < MAX_ITER_ROOT
        if ((ss - sr) * fps - fs) * ((ss - sl) * fps - fs) > 0.0 ||
           abs(2.0 * fs) > abs(dsold * fps)
            dsold = dss
            dss = 0.5 * (sr - sl)
            ss = sl + dss
            gensec = 0
        else
            dsold = dss
            dss = fs / fps
            ss -= dss
        end
        iter += 1
        if abs(dss) < EPS_ROOT
            not_conv = false
            s0[4] = f_sign * fps
        end

        if not_conv
            for i in 1:NDIM
                xs[i] = x0[i] + ss * dir[i]
            end
            fs = f_sign * call_integrand(impl_func, par, xs)
            fps = (fs - fv[1]) / (ss - sv[3])
            if fs < 0.0
                sl = ss
                fl = fs
            elseif fs > 0.0
                sr = ss
                fr = fs
            else
                not_conv = false
                s0[4] = f_sign * fps
            end
            sv = [sv[2], sv[3], ss]
            ds2 = sv[3] - sv[1]
            if gensec > 0 && abs(ds2) > EPS_ROOT
                fv[3] = (fps - fv[2]) / ds2
            else
                fv[3] = 0.0
            end
            fv = [fs, fps, fv[3]]
            gensec = 1
            fps = fv[2] + fv[3] * (sv[3] - sv[2])
        end
    end

    if !not_conv
        sz = f_sign * ss + 0.5 * (1 - f_sign) * s0[1]
    else
        s1 = 0.0
        f1 = f_sign * call_integrand(impl_func, par, x0)
        s2 = s0[1]
        for i in 1:NDIM
            xs[i] = x0[i] + s2 * dir[i]
        end
        f2 = f_sign * call_integrand(impl_func, par, xs)
        if f1 * f2 <= 0.0
            if sl > 0.0
                s1 = sl
                f1 = fl
            end
            if sr < s2
                s2 = sr
                f2 = fr
            end
            ss = s1 - f1 * (s2 - s1) / (f2 - f1)
            sz = f_sign * ss + 0.5 * (1 - f_sign) * s0[1]
            s0[4] = f_sign * fps
        else
            sz = f_sign * s0[1]  # no zero found
            s0[4] = 0.0
        end
    end

    return sz
end

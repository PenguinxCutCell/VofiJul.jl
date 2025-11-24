function vofi_check_boundary_line(impl_func, par, x0, h0, f0, xfs_pt, n0)
    nx = @MVector ones(vofi_int, NSE)
    ny = @MVector ones(vofi_int, NSE)
    sidedirx = (1.0, 0.0, 0.0)
    sidediry = (0.0, 1.0, 0.0)
    fse = @MVector zeros(vofi_real, NSE)
    x1 = @MVector zeros(vofi_real, NDIM)
    xfsl = MinData()
    check_dir = -1

    for i in 0:1
        for j in 0:1
            if n0[i + 1, j + 1] > 0
                if ny[i + 1] > 0
                    ny[i + 1] = 0
                    fse[1] = f0[i + 1, 1]
                    fse[2] = f0[i + 1, 2]
                    for k in 1:NDIM
                        x1[k] = x0[k] + i * sidedirx[k] * h0[1]
                    end
                    consi = vofi_check_side_consistency(impl_func, par, x1, sidediry, fse, h0[2])
                    if consi != 0
                        f2pos = consi
                        sign_change = vofi_get_segment_min(impl_func, par, x1, sidediry, fse, xfsl, h0[2], f2pos)
                        if sign_change != 0
                            copy!(xfs_pt, xfsl)
                            xfs_pt.isc[1] = 1
                            xfs_pt.isc[i + 2] = 1
                            check_dir = 0
                        end
                    end
                end
                if nx[j + 1] > 0
                    nx[j + 1] = 0
                    fse[1] = f0[1, j + 1]
                    fse[2] = f0[2, j + 1]
                    for k in 1:NDIM
                        x1[k] = x0[k] + j * sidediry[k] * h0[2]
                    end
                    consi = vofi_check_side_consistency(impl_func, par, x1, sidedirx, fse, h0[1])
                    if consi != 0
                        f2pos = consi
                        sign_change = vofi_get_segment_min(impl_func, par, x1, sidedirx, fse, xfsl, h0[1], f2pos)
                        if sign_change != 0
                            copy!(xfs_pt, xfsl)
                            xfs_pt.isc[1] = 1
                            xfs_pt.isc[j + 2] = 1
                            check_dir = 1
                        end
                    end
                end
                n0[i + 1, j + 1] = 0
            end
        end
    end

    return check_dir
end

function vofi_check_boundary_surface(impl_func, par, x0, h0, f0, xfs, n0)
    nx = @MVector ones(vofi_int, NSE)
    ny = @MVector ones(vofi_int, NSE)
    nz = @MVector ones(vofi_int, NSE)
    sidedirx = (1.0, 0.0, 0.0)
    sidediry = (0.0, 1.0, 0.0)
    sidedirz = (0.0, 0.0, 1.0)
    fve = @MVector zeros(vofi_real, NVER)
    x1 = @MVector zeros(vofi_real, NDIM)
    xfsl = MinData()
    check_dir = -1

    for i in 0:1
        for j in 0:1
            for k in 0:1
                if n0[i + 1, j + 1, k + 1] > 0
                    if nx[i + 1] > 0
                        nx[i + 1] = 0
                        fve[1] = f0[i + 1, 1, 1]
                        fve[2] = f0[i + 1, 2, 1]
                        fve[3] = f0[i + 1, 1, 2]
                        fve[4] = f0[i + 1, 2, 2]
                        for m in 1:NDIM
                            x1[m] = x0[m] + i * sidedirx[m] * h0[1]
                        end
                        ipsc = vofi_check_face_consistency(impl_func, par, x1, h0, sidediry, sidedirz, fve)
                        if ipsc.consi != 0
                            sign_change = vofi_get_face_min(impl_func, par, x1, h0, sidediry, sidedirz, fve, xfsl, ipsc)
                            if sign_change != 0
                                copy!(xfs[1], xfsl)
                                xfs[1].isc[1] = 1
                                xfs[1].isc[i + 2] = 1
                                check_dir = 0
                            end
                        end
                    end
                    if ny[j + 1] > 0
                        ny[j + 1] = 0
                        fve[1] = f0[1, j + 1, 1]
                        fve[2] = f0[2, j + 1, 1]
                        fve[3] = f0[1, j + 1, 2]
                        fve[4] = f0[2, j + 1, 2]
                        for m in 1:NDIM
                            x1[m] = x0[m] + j * sidediry[m] * h0[2]
                        end
                        ipsc = vofi_check_face_consistency(impl_func, par, x1, h0, sidedirx, sidedirz, fve)
                        if ipsc.consi != 0
                            sign_change = vofi_get_face_min(impl_func, par, x1, h0, sidedirx, sidedirz, fve, xfsl, ipsc)
                            if sign_change != 0
                                copy!(xfs[2], xfsl)
                                xfs[2].isc[1] = 1
                                xfs[2].isc[j + 2] = 1
                                check_dir = 0
                            end
                        end
                    end
                    if nz[k + 1] > 0
                        nz[k + 1] = 0
                        fve[1] = f0[1, 1, k + 1]
                        fve[2] = f0[2, 1, k + 1]
                        fve[3] = f0[1, 2, k + 1]
                        fve[4] = f0[2, 2, k + 1]
                        for m in 1:NDIM
                            x1[m] = x0[m] + k * sidedirz[m] * h0[3]
                        end
                        ipsc = vofi_check_face_consistency(impl_func, par, x1, h0, sidedirx, sidediry, fve)
                        if ipsc.consi != 0
                            sign_change = vofi_get_face_min(impl_func, par, x1, h0, sidedirx, sidediry, fve, xfsl, ipsc)
                            if sign_change != 0
                                copy!(xfs[3], xfsl)
                                xfs[3].isc[1] = 1
                                xfs[3].isc[k + 2] = 1
                                check_dir = 0
                            end
                        end
                    end
                    n0[i + 1, j + 1, k + 1] = 0
                end
            end
        end
    end

    return check_dir
end

function vofi_check_secondary_side(impl_func, par, x0, h0, pdir, sdir, f0, xfs_pt, fth)
    hs = zero(vofi_real)
    for i in 1:NDIM
        hs += sdir[i] * h0[i]
    end
    x1 = @MVector zeros(vofi_real, NDIM)
    fse = @MVector zeros(vofi_real, NSE)
    xfsl = MinData()

    for k in 0:1
        fse[1] = f0[k + 1, 1]
        fse[2] = f0[k + 1, 2]
        if fse[1] * fse[2] < 0
            xfs_pt.isc[1] = 1
            xfs_pt.isc[k + 2] = -1
        else
            if abs(fse[1]) > fth && abs(fse[2]) > fth
                continue
            end
            for i in 1:NDIM
                x1[i] = x0[i] + k * pdir[i] * h0[i]
            end
            consi = vofi_check_side_consistency(impl_func, par, x1, sdir, fse, hs)
            if consi != 0
                f2pos = consi
                sign_change = vofi_get_segment_min(impl_func, par, x1, sdir, fse, xfsl, hs, f2pos)
                if sign_change != 0
                    copy!(xfs_pt, xfsl)
                    xfs_pt.isc[1] = 1
                    xfs_pt.isc[k + 2] = 1
                end
            end
        end
    end
end

function vofi_check_secter_face(impl_func, par, x0, h0, pdir, sdir, tdir, f0, xfs_pt, fth)
    xfs_pt.isc .= 0
    x1 = @MVector zeros(vofi_real, NDIM)
    fve = @MVector zeros(vofi_real, NVER)
    xfsl = MinData()

    for m in 0:1
        np0 = nm0 = 0
        fve[1] = f0[m + 1, 1, 1]; np0 += fve[1] > 0; nm0 += fve[1] < 0
        fve[2] = f0[m + 1, 2, 1]; np0 += fve[2] > 0; nm0 += fve[2] < 0
        fve[3] = f0[m + 1, 1, 2]; np0 += fve[3] > 0; nm0 += fve[3] < 0
        fve[4] = f0[m + 1, 2, 2]; np0 += fve[4] > 0; nm0 += fve[4] < 0

        if nm0 * np0 == 0
            if !(abs(fve[1]) > fth && abs(fve[2]) > fth && abs(fve[3]) > fth && abs(fve[4]) > fth)
                for i in 1:NDIM
                    x1[i] = x0[i] + m * pdir[i] * h0[i]
                end
                ipsc = vofi_check_face_consistency(impl_func, par, x1, h0, sdir, tdir, fve)
                if ipsc.consi != 0
                    sign_change = vofi_get_face_min(impl_func, par, x1, h0, sdir, tdir, fve, xfsl, ipsc)
                    if sign_change != 0
                        copy!(xfs_pt, xfsl)
                        xfs_pt.isc[1] = 1
                        xfs_pt.isc[m + 1] = 1
                    end
                end
            end
        end
    end
    return nothing
end

function vofi_check_tertiary_side(impl_func, par, x0, h0, pdir, sdir, tdir, f0, xfs, fth)
    ht = zero(vofi_real)
    for i in 1:NDIM
        ht += tdir[i] * h0[i]
    end
    for idx in 1:4
        xfs[idx].isc .= 0
    end
    x1 = @MVector zeros(vofi_real, NDIM)
    fse = @MVector zeros(vofi_real, NSE)

    for m in 0:1
        for n in 0:1
            fse[1] = f0[m + 1, n + 1, 1]
            fse[2] = f0[m + 1, n + 1, 2]
            l0 = 2 * m + n + 1
            if fse[1] * fse[2] < 0
                xfs[l0].isc[1] = 1
                xfs[l0].isc[2] = -1
            else
                if !(abs(fse[1]) > fth && abs(fse[2]) > fth)
                    for i in 1:NDIM
                        x1[i] = x0[i] + m * pdir[i] * h0[i] + n * sdir[i] * h0[i]
                    end
                    consi = vofi_check_side_consistency(impl_func, par, x1, tdir, fse, ht)
                    if consi != 0
                        xfsl = MinData()
                        f2pos = consi
                        sign_change = vofi_get_segment_min(impl_func, par, x1, tdir, fse, xfsl, ht, f2pos)
                        if sign_change != 0
                            xfs[l0].xval .= xfsl.xval
                            xfs[l0].sval = xfsl.sval
                            xfs[l0].fval = xfsl.fval
                            xfs[l0].isc[1] = 1
                            xfs[l0].isc[2] = 1
                        end
                    end
                end
            end
        end
    end
    return nothing
end

function vofi_check_boundary_hypersurface(impl_func, par, x0::Vector{Float64}, h0::Vector{Float64}, f0::Array{Float64, 4}, xfs, n0::Array{Int64, 4})
    # For 4D hypercube, we have 8 cubic (3D) faces
    # Each face is defined by fixing one coordinate at 0 or 1
    nx = ones(Int, 2)
    ny = ones(Int, 2)
    nz = ones(Int, 2)
    nw = ones(Int, 2)
    
    sidedirx = [1.0, 0.0, 0.0, 0.0]
    sidediry = [0.0, 1.0, 0.0, 0.0]
    sidedirz = [0.0, 0.0, 1.0, 0.0]
    sidedirw = [0.0, 0.0, 0.0, 1.0]
    
    fcube = Array{Float64}(undef, 2, 2, 2)  # 8 vertices of a cube face
    x1 = Vector{Float64}(undef, 4)
    check_dir = -1
    
    # Check faces perpendicular to x-axis (i fixed)
    for i in 0:1
        if nx[i + 1] > 0
            nx[i + 1] = 0
            for j in 0:1, k in 0:1, l in 0:1
                if n0[i + 1, j + 1, k + 1, l + 1] > 0
                    # Extract the cubic face values
                    for jj in 0:1, kk in 0:1, ll in 0:1
                        fcube[jj + 1, kk + 1, ll + 1] = f0[i + 1, jj + 1, kk + 1, ll + 1]
                    end
                    
                    # Set position on the face
                    x1[1] = x0[1] + i * h0[1]
                    x1[2] = x0[2]
                    x1[3] = x0[3]
                    x1[4] = x0[4]
                    
                    # Check if this face has sign changes
                    if has_sign_change_3d(fcube)
                        check_dir = 0
                        break
                    end
                end
            end
        end
    end
    
    # Check faces perpendicular to y-axis (j fixed)
    for j in 0:1
        if ny[j + 1] > 0
            ny[j + 1] = 0
            for i in 0:1, k in 0:1, l in 0:1
                if n0[i + 1, j + 1, k + 1, l + 1] > 0
                    for ii in 0:1, kk in 0:1, ll in 0:1
                        fcube[ii + 1, kk + 1, ll + 1] = f0[ii + 1, j + 1, kk + 1, ll + 1]
                    end
                    
                    x1[1] = x0[1]
                    x1[2] = x0[2] + j * h0[2]
                    x1[3] = x0[3]
                    x1[4] = x0[4]
                    
                    if has_sign_change_3d(fcube)
                        check_dir = 0
                        break
                    end
                end
            end
        end
    end
    
    # Check faces perpendicular to z-axis (k fixed)
    for k in 0:1
        if nz[k + 1] > 0
            nz[k + 1] = 0
            for i in 0:1, j in 0:1, l in 0:1
                if n0[i + 1, j + 1, k + 1, l + 1] > 0
                    for ii in 0:1, jj in 0:1, ll in 0:1
                        fcube[ii + 1, jj + 1, ll + 1] = f0[ii + 1, jj + 1, k + 1, ll + 1]
                    end
                    
                    x1[1] = x0[1]
                    x1[2] = x0[2]
                    x1[3] = x0[3] + k * h0[3]
                    x1[4] = x0[4]
                    
                    if has_sign_change_3d(fcube)
                        check_dir = 0
                        break
                    end
                end
            end
        end
    end
    
    # Check faces perpendicular to w-axis (l fixed)
    for l in 0:1
        if nw[l + 1] > 0
            nw[l + 1] = 0
            for i in 0:1, j in 0:1, k in 0:1
                if n0[i + 1, j + 1, k + 1, l + 1] > 0
                    for ii in 0:1, jj in 0:1, kk in 0:1
                        fcube[ii + 1, jj + 1, kk + 1] = f0[ii + 1, jj + 1, kk + 1, l + 1]
                    end
                    
                    x1[1] = x0[1]
                    x1[2] = x0[2]
                    x1[3] = x0[3]
                    x1[4] = x0[4] + l * h0[4]
                    
                    if has_sign_change_3d(fcube)
                        check_dir = 0
                        break
                    end
                end
            end
        end
    end
    
    return check_dir
end

# Helper function to check if a 3D cube has sign changes
function has_sign_change_3d(fcube::Array{Float64, 3})
    np = 0
    nm = 0
    for i in 1:2, j in 1:2, k in 1:2
        if fcube[i, j, k] > 0
            np += 1
        elseif fcube[i, j, k] < 0
            nm += 1
        end
    end
    return np > 0 && nm > 0
end

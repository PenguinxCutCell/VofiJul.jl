pad_to_ndim(vec) = begin
    out = @MVector zeros(vofi_real, NDIM)
    for i in 1:min(length(vec), NDIM)
        out[i] = vofi_real(vec[i])
    end
    out
end

@inline function axis_index(dir)
    for i in 1:length(dir)
        if dir[i] != 0
            return i
        end
    end
    throw(ArgumentError("direction vector must have a non-zero component"))
end

@inline axis_length(dir, hvec) = hvec[axis_index(dir)]

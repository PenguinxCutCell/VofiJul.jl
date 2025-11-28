@inline function pad_to_ndim(vec)
    out = @MVector zeros(vofi_real, NDIM)
    @inbounds for i in 1:min(length(vec), NDIM)
        out[i] = vofi_real(vec[i])
    end
    out
end

@inline function axis_index(dir)
    @inbounds for i in 1:length(dir)
        if dir[i] != 0
            return i
        end
    end
    throw(ArgumentError("direction vector must have a non-zero component"))
end

@inline axis_length(dir, hvec) = @inbounds hvec[axis_index(dir)]

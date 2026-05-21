const TIMEKEEPERS_LOGO_PATH = normpath(joinpath(@__DIR__, "..", "images", "timekeepers-logo.svg"))
const _TIMEKEEPERS_ICON_CACHE = Ref{Union{Nothing, Vector{Matrix{NTuple{4, UInt8}}}}}(nothing)

_tk_byte(x::Real) = UInt8(clamp(round(Int, x), 0, 255))
_tk_rgba(r, g, b, a = 255) = (_tk_byte(r), _tk_byte(g), _tk_byte(b), _tk_byte(a))

function _tk_over(dst::NTuple{4, UInt8}, src::NTuple{4, UInt8}, alpha::Real)
    a = clamp(Float64(alpha) * Float64(src[4]) / 255.0, 0.0, 1.0)
    inv = 1.0 - a
    out_a = a + Float64(dst[4]) / 255.0 * inv
    out_a <= 0 && return _tk_rgba(0, 0, 0, 0)
    return _tk_rgba(
        (Float64(src[1]) * a + Float64(dst[1]) * inv) / out_a,
        (Float64(src[2]) * a + Float64(dst[2]) * inv) / out_a,
        (Float64(src[3]) * a + Float64(dst[3]) * inv) / out_a,
        255 * out_a,
    )
end

function _tk_stroke_alpha(distance::Real, width::Real)
    half = width / 2
    feather = width / 3
    return clamp((half + feather - abs(distance)) / feather, 0.0, 1.0)
end

function _tk_rounded_rect_alpha(x::Real, y::Real)
    cx = abs(x - 0.5)
    cy = abs(y - 0.5)
    half = 0.43
    radius = 0.20
    qx = max(cx - (half - radius), 0.0)
    qy = max(cy - (half - radius), 0.0)
    outside = sqrt(qx * qx + qy * qy) - radius
    return clamp((0.018 - outside) / 0.018, 0.0, 1.0)
end

function _timekeepers_icon(size::Integer)
    n = Int(size)
    n > 0 || error("Icon size must be positive")
    img = fill(_tk_rgba(0, 0, 0, 0), n, n)
    bg = _tk_rgba(32, 37, 43, 255)
    teal = _tk_rgba(46, 176, 199, 255)
    white = _tk_rgba(255, 255, 255, 255)
    coral = _tk_rgba(230, 90, 84, 255)
    @inbounds for row in 1:n, col in 1:n
        x = (col - 0.5) / n
        y = (row - 0.5) / n
        px = _tk_rgba(0, 0, 0, 0)

        bg_alpha = _tk_rounded_rect_alpha(x, y)
        px = _tk_over(px, bg, bg_alpha)

        ring = _tk_stroke_alpha(sqrt((x - 0.5)^2 + (y - 0.5)^2) - 0.242, 0.066)
        px = _tk_over(px, teal, ring)

        if 0.24 <= x <= 0.76
            t = (x - 0.24) / 0.52
            wave_y = 0.515 - 0.17 * sin(2π * t)
            wave = _tk_stroke_alpha(y - wave_y, 0.058)
            px = _tk_over(px, white, wave)
        end

        dot = 1.0 - clamp((sqrt((x - 0.76)^2 + (y - 0.515)^2) - 0.050) / 0.018, 0.0, 1.0)
        px = _tk_over(px, coral, dot)
        img[row, col] = px
    end
    return img
end

function _timekeepers_icons()
    cached = _TIMEKEEPERS_ICON_CACHE[]
    cached === nothing || return cached
    icons = [_timekeepers_icon(size) for size in (16, 32, 64, 128)]
    _TIMEKEEPERS_ICON_CACHE[] = icons
    return icons
end

function _tk_write_be(io::IO, x::UInt32)
    write(io, UInt8((x >> 24) & 0xff))
    write(io, UInt8((x >> 16) & 0xff))
    write(io, UInt8((x >> 8) & 0xff))
    write(io, UInt8(x & 0xff))
    return io
end

function _tk_crc32(bytes::Vector{UInt8})
    crc = UInt32(0xffffffff)
    for byte in bytes
        crc ⊻= UInt32(byte)
        for _ in 1:8
            crc = isodd(crc) ? (crc >> 1) ⊻ UInt32(0xedb88320) : crc >> 1
        end
    end
    return ~crc
end

function _tk_adler32(bytes::Vector{UInt8})
    a = UInt32(1)
    b = UInt32(0)
    for byte in bytes
        a = (a + UInt32(byte)) % UInt32(65521)
        b = (b + a) % UInt32(65521)
    end
    return (b << 16) | a
end

function _tk_write_chunk(io::IO, kind::String, data::Vector{UInt8})
    kind_bytes = Vector{UInt8}(codeunits(kind))
    _tk_write_be(io, UInt32(length(data)))
    write(io, kind_bytes)
    write(io, data)
    _tk_write_be(io, _tk_crc32([kind_bytes; data]))
    return io
end

function _tk_zlib_store(bytes::Vector{UInt8})
    out = IOBuffer()
    write(out, UInt8(0x78), UInt8(0x01))
    offset = 1
    while offset <= length(bytes)
        n = min(65535, length(bytes) - offset + 1)
        final = offset + n > length(bytes)
        write(out, final ? UInt8(0x01) : UInt8(0x00))
        len = UInt16(n)
        nlen = ~len
        write(out, UInt8(len & 0xff), UInt8((len >> 8) & 0xff))
        write(out, UInt8(nlen & 0xff), UInt8((nlen >> 8) & 0xff))
        write(out, view(bytes, offset:(offset + n - 1)))
        offset += n
    end
    _tk_write_be(out, _tk_adler32(bytes))
    return take!(out)
end

function _write_timekeepers_icon_png(path::AbstractString; icon_size::Integer = 256)
    icon = _timekeepers_icon(icon_size)
    n = size(icon, 1)
    raw = UInt8[]
    sizehint!(raw, n * (4n + 1))
    @inbounds for row in 1:n
        push!(raw, 0x00)
        for col in 1:n
            px = icon[row, col]
            push!(raw, px[1], px[2], px[3], px[4])
        end
    end

    ihdr = IOBuffer()
    _tk_write_be(ihdr, UInt32(n))
    _tk_write_be(ihdr, UInt32(n))
    write(ihdr, UInt8(8), UInt8(6), UInt8(0), UInt8(0), UInt8(0))

    open(path, "w") do io
        write(io, UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a])
        _tk_write_chunk(io, "IHDR", take!(ihdr))
        _tk_write_chunk(io, "IDAT", _tk_zlib_store(raw))
        _tk_write_chunk(io, "IEND", UInt8[])
    end
    return path
end

function _objc_class(name::String)
    return ccall((:objc_getClass, "/usr/lib/libobjc.A.dylib"), Ptr{Cvoid}, (Cstring,), name)
end

function _objc_sel(name::String)
    return ccall((:sel_registerName, "/usr/lib/libobjc.A.dylib"), Ptr{Cvoid}, (Cstring,), name)
end

function _objc_msg(receiver::Ptr{Cvoid}, selector::Ptr{Cvoid})
    return ccall((:objc_msgSend, "/usr/lib/libobjc.A.dylib"), Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}), receiver, selector)
end

function _objc_msg(receiver::Ptr{Cvoid}, selector::Ptr{Cvoid}, arg::Ptr{Cvoid})
    return ccall((:objc_msgSend, "/usr/lib/libobjc.A.dylib"), Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), receiver, selector, arg)
end

function _objc_msg(receiver::Ptr{Cvoid}, selector::Ptr{Cvoid}, arg::AbstractString)
    return ccall((:objc_msgSend, "/usr/lib/libobjc.A.dylib"), Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Cstring), receiver, selector, arg)
end

function _objc_msg_void(receiver::Ptr{Cvoid}, selector::Ptr{Cvoid}, arg::Ptr{Cvoid})
    return ccall((:objc_msgSend, "/usr/lib/libobjc.A.dylib"), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}), receiver, selector, arg)
end

function _nsstring(value::AbstractString)
    cls = _objc_class("NSString")
    return _objc_msg(cls, _objc_sel("stringWithUTF8String:"), String(value))
end

function _apply_macos_dock_icon!()
    Sys.isapple() || return false
    try
        icon_path = joinpath(tempdir(), "timekeepers-icon.png")
        _write_timekeepers_icon_png(icon_path; icon_size = 256)
        app = _objc_msg(_objc_class("NSApplication"), _objc_sel("sharedApplication"))
        image_alloc = _objc_msg(_objc_class("NSImage"), _objc_sel("alloc"))
        image = _objc_msg(image_alloc, _objc_sel("initWithContentsOfFile:"), _nsstring(icon_path))
        image == C_NULL && return false
        _objc_msg_void(app, _objc_sel("setApplicationIconImage:"), image)
        return true
    catch err
        @debug "Could not set macOS Dock icon" exception = err
        return false
    end
end

function _apply_timekeepers_icon!(screen)
    if Sys.isapple()
        _apply_macos_dock_icon!()
    else
        try
            GLMakie.GLFW.SetWindowIcon(screen.glscreen, _timekeepers_icons())
        catch err
            @debug "Could not set GLFW window icon" exception = err
        end
    end
    return screen
end

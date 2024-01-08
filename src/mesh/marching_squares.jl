struct FermiSurfaceSegment
    points::Vector{SVector{2,Float64}}
    isclosed::Bool
    length::Real
end

struct FermiSurface
    segments::Vector{FermiSurfaceSegment}
end

T, B, R, L = 0x01, 0x02, 0x04, 0x08
crossing_lookup = [L|B, B|R, L|R, T|R, 0x0, T|B, L|T, L|T, T|B, 0x0, T|R, L|R, B|R, L|B]

function get_cells(A::AbstractMatrix, level::Real = 0.0)
    x, y = axes(A)
    cells = Dict{Tuple{Int, Int}, UInt8}()

    @inbounds for i in first(x):last(x)-1
        for j in first(y):last(y)-1
            intersect = (A[i, j] > level) ? 0x01 : 0x00
            (A[i + 1, j] > level) && (intersect |= 0x02)
            (A[i + 1, j + 1] > level) && (intersect |= 0x04)
            (A[i, j + 1] > level) && (intersect |= 0x08)

            if !(intersect == 0x00) && !(intersect == 0x0f)
                cells[(i,j)] = crossing_lookup[intersect]
            end

        end
    end 

    return cells
end

shift = [(0,1), (0,-1), (1,0), (-1,0)]
new_edge = (B,T,L,R)
function get_next_cell(edge, index)
    i = trailing_zeros(edge) + 1
    index = index .+ shift[i]
   return new_edge[i], index
end

function find_contour(x, y, A::AbstractMatrix)
    fs = FermiSurface(FermiSurfaceSegment[])

    cells = get_cells(A)

    xax, yax = axes(A) # Valid indices for iterating over
    is = first(xax):last(xax)-1
    js = first(yax):last(yax)-1

    while Base.length(cells) > 0
        segment = Vector{SVector{2,Float64}}(undef, 0)
        start_index, case = first(cells)
        start_edge = 0x01 << trailing_zeros(case)

        end_index = follow_contour!(cells, segment, x, y, A, is, js, start_index, start_edge)

        # Check if the contour forms a loop
        isclosed = end_index == start_index ? true : false

        if !isclosed
            # Go back to the starting cell and walk the other direction
            edge, index = get_next_cell(start_edge, start_index)
            follow_contour!(cells, reverse!(segment), x, y, A, is, js, index, edge)
        end

        seg_length = 0.0
        for i in eachindex(segment)
            i == length(segment) && continue
            seg_length += norm(segment[i + 1] - segment[i])
        end
        isclosed && (seg_length += norm(last(segment) - first(segment)))

        push!(fs.segments, FermiSurfaceSegment(segment, isclosed, seg_length))
    end
    
    return fs
end

function follow_contour!(cells, contour, x, y, A, is, js, start_index, start_edge)
    index = start_index
    edge = start_edge

    push!(contour, get_crossing(x, y, A, index, edge))
    while true
        edge = pop!(cells, index) ⊻ edge
        
        push!(contour, get_crossing(x, y, A, index, edge))
        
        edge, index = get_next_cell(edge, index)

        (index == start_index || !(index[1] ∈ is) || !(index[2] ∈ js)) && break
    end

    return index
end

function get_crossing(x, y, A, index, edge)
    i, j = index
    if edge == 0x08 # Left
        xcoord = x[i]
        ycoord = y[j] - A[i,j] * (y[j + 1] - y[j]) / (A[i, j + 1] - A[i,j])
    elseif edge == 0x01 # Top
        xcoord = x[i] - A[i, j + 1] * (x[i + 1] - x[i]) / (A[i + 1, j + 1] - A[i, j + 1])
        ycoord = y[j + 1]
    elseif edge == 0x04 # Right
        xcoord = x[i + 1]
        ycoord = y[j] - A[i + 1, j] * (y[j + 1] - y[j]) / (A[i + 1, j + 1] - A[i + 1, j])
    elseif edge == 0x02 # Bottom
        xcoord = x[i] - A[i,j] * (x[i + 1] - x[i]) / (A[i + 1, j] - A[i,j])
        ycoord = y[j]
    end

    return SVector(xcoord, ycoord)
end

function uniform_fermi_surface(fs::FermiSurface, n::Int)
    perimeter = sum(map(x -> x.length, fs.segments))
    
    uniform_fs = FermiSurface(FermiSurfaceSegment[])
    for (i, segment) in enumerate(fs.segments)
        seg_n = round(Int, n * segment.length / perimeter)
        if segment.isclosed
            seg_range = LinRange(0.0, (seg_n - 1) * segment.length / seg_n, seg_n)
        else
            seg_range = LinRange(0.0, segment.length, seg_n + 1)
        end
        
        uniform_segment = [first(segment.points)]

        j = 1
        t = 0.0
        jmax = length(segment.points) 
        for s in seg_range
            s == 0.0 && continue
            while t < s && j < jmax
                t += norm(segment.points[j + 1] - segment.points[j])
                j += 1
            end
            # Go back one point after finding
            Δ = segment.points[j] - segment.points[j-1]
            t -= norm(Δ)
            j -= 1

            k_interp = segment.points[j - 1] + Δ * (s - t) / norm(Δ) 
            push!(uniform_segment, k_interp)
        end

        push!(uniform_fs.segments, FermiSurfaceSegment(uniform_segment, segment.isclosed, segment.length))
    end

    return uniform_fs
end




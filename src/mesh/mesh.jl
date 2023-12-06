using Contour
import ForwardDiff: gradient

export get_fermi_tube

function get_fermi_tube(f::Function, T::Real, n::Int = 1000)
    grid = LinRange(-0.5, 0.5, n) # Sample points in kx, ky in units of (2pi/a)
    energies = map(x -> f([x[1], x[2]]), collect(Iterators.product(grid, grid)))

    fs = find_contour(grid, grid, energies)

    for segment in fs.segments
        vel = gradient.(f, segment)
        for i in eachindex(segment)
            
        end

    end

    return fs
end
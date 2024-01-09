using FSGenerator
using CairoMakie
import StaticArrays: SVector
import LinearAlgebra: det

function main()
    grid_vals = Base._linspace(-0.5, 0.5, 2001)
    grid = map(x -> [x[1], x[2]], collect(Iterators.product(grid_vals, grid_vals)))
    es = hamiltonian.(grid)
    T = 12/Tf
    fd = f0.(es, T) .* (1 .- f0.(es, T))

    fs = FSGenerator.find_contour(grid_vals, grid_vals, es)
    @show first(fs.segments[1].points)
    uniform_fs = FSGenerator.uniform_fermi_surface(fs, 140)
    @time mesh, _ = temperature_broaden(uniform_fs, hamiltonian, T, 15)
    grid_momenta = vec(map(x -> x.momentum, mesh.patches))

    fig = Figure(size = (1000,1000))
    ax  = Axis(fig[1,1], aspect = 1.0)
    # heatmap!(ax, grid_vals, grid_vals, fd, marker = :rect, markersize= 5, colormap = :lisbon)
    
    # scatter!(ax, first.(vec(corners)), last.(vec(corners)), markersize = 2.0, color = :green)
    
    c = [map(y -> mesh.corners[y], x.corners) for x in mesh.patches]
    p = poly!(ax, c, color = map(x -> x.energy, mesh.patches), colormap = :berlin)
    Colorbar(fig[1,2], p)

    scatter!(ax, first.(grid_momenta), last.(grid_momenta), markersize = 1.0, color = :white)
    for segment in uniform_fs.segments
        lines!(ax, first.(segment.points), last.(segment.points), color = :green)
    end
    display(fig)

end

include(joinpath(@__DIR__, "bands", "Sr2RuO4", "gamma.jl"))
# include(joinpath(@__DIR__, "bands", "tbm.jl"))

main()
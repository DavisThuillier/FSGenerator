using FSGenerator
using CairoMakie
import StaticArrays: SVector
import LinearAlgebra: det

function main()
    grid_vals = Base._linspace(-0.5, 0.5, 1001)
    grid = map(x -> [x[1], x[2]], collect(Iterators.product(grid_vals, grid_vals)))
    es = hamiltonian.(grid)
    T = 15/Tf
    fd = f0.(es, T) .* (1 .- f0.(es, T))

    fs = FSGenerator.find_contour(grid_vals, grid_vals, es)
    uniform_fs = FSGenerator.uniform_fermi_surface(fs, 201)
    @time grid, corners = temperature_broaden(uniform_fs, hamiltonian, T, 5)
    grid_momenta = vec(map(x -> x.momentum, grid))

    fig = Figure()
    ax  = Axis(fig[1,1], aspect = 1.0)
    heatmap!(ax, grid_vals, grid_vals, fd, marker = :rect, markersize= 5, colormap = :lisbon)
    scatter!(ax, first.(grid_momenta), last.(grid_momenta), markersize = 2.0, color = vec(map(x -> det(x.jinv), grid)), colormap = :thermal)
    scatter!(ax, first.(vec(corners)), last.(vec(corners)), markersize = 2.0, color = :green)
    # for segment in uniform_fs.segments
    #     scatter!(ax, first.(segment.points), last.(segment.points), markersize = 4.0, color = :black)
    # end
    display(fig)

end

include(joinpath(@__DIR__, "bands", "Sr2RuO4", "gamma.jl"))
# include(joinpath(@__DIR__, "bands", "tbm.jl"))

main()
using FSGenerator
import ForwardDiff: gradient
import ProgressBars: ProgressBar

function main()
    grid_vals = Base._linspace(-0.5, 0.5, 1001)
    grid = map(x -> [x[1], x[2]], collect(Iterators.product(grid_vals, grid_vals)))
    es = hamiltonian.(grid)
    T = 15/Tf
    fd = f0.(es, T) .* (1 .- f0.(es, T))

    fs = FSGenerator.find_contour(grid_vals, grid_vals, es)
    uniform_fs = FSGenerator.uniform_fermi_surface(fs, 51)
    @time grid, corners = temperature_broaden(uniform_fs, hamiltonian, T, 5)
    grid[1]

    a = Vector{Float64}(undef, 6)
    @time for i in ProgressBar(eachindex(grid))
        ki = grid[i].momentum
        Ji = grid[i].jinv
        for j in eachindex(grid)
            kj = grid[j].momentum
            Jj = grid[j].jinv
            for ℓ in eachindex(grid)
                kℓ = grid[ℓ].momentum
                Jℓ = grid[ℓ].jinv
                a[1], a[2] = gradient(x -> hamiltonian(x + kj - kℓ), ki)' * Ji
                a[3], a[4] = gradient(x -> hamiltonian(ki + x - kℓ), kj)' * Jj
                a[5], a[6] = gradient(x -> hamiltonian(ki + kj - x), kℓ)' * Jℓ
                hamiltonian(ki + kj - kℓ)
            end
        end
    end

end

include(joinpath(@__DIR__, "bands", "Sr2RuO4", "gamma.jl"))

main()
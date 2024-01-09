
export temperature_broaden

function get_k_bound(hamiltonian::Function, e_bound::Float64, fs_k::SVector{2, Float64}, velocity::SVector{2,Float64}; max_iterations = 10000, tolerance = 0.0001)
    e_bound  == 0.0 && return fs_k

    n = velocity / norm(velocity)
    step = (e_bound / norm(velocity))

    i::Int = 0
    i_limit::Int = abs(div(0.5, step)) # Number of step corresponding to half-width of Brillouin zone

    endpoint = fs_k
    startpoint = fs_k
    
    if e_bound > 0
        while hamiltonian(endpoint) < e_bound && i < i_limit
            endpoint += step * n
            (abs(endpoint[1]) > 0.5 || abs(endpoint[2]) > 0.5) && break
            i += 1
        end
        endpoint += step * n
    else
        while hamiltonian(startpoint) > e_bound && i < i_limit
            startpoint += step * n
            i += 1
        end
        startpoint += step * n 
    end

    j::Int = 0
    while j < max_iterations
        midpoint = (startpoint + endpoint) / 2
        delta_E = hamiltonian(midpoint) - e_bound
        norm(startpoint - endpoint) < tolerance && break
        
        if sign(delta_E) == sign(hamiltonian(startpoint) - e_bound)
            startpoint = midpoint
        else
            endpoint   = midpoint
        end
        j += 1
    end

    k_bound = (startpoint + endpoint) / 2

    # Check if k_bound lies outside the Brillouin zone
    if abs(k_bound[1]) > 0.5
        k_bound = fs_k + ( (0.5 * sign(k_bound[1]) - fs_k[1]) / n[1]) * n
    elseif abs(k_bound[2]) > 0.5
        k_bound = fs_k + ( (0.5 * sign(k_bound[2]) - fs_k[2]) / n[2]) * n
    end

    return k_bound
end

# Integration variable containing the momentum, energy, and size of a mesh patch
struct Patch
    momentum::SVector{2,Float64}
    energy::Float64
    dV::Float64
    jinv::Matrix{Float64} # Jacobian of transformation from (kx, ky) --> (E, s)
    corners::Vector{Int}
end

Patch(momentum, f::Function, dV, corners) = Patch(momentum, f(momentum), dV, get_jacobian_inverse(momentum, f), corners)

struct Mesh
    patches::Vector{Patch}
    corners::Vector{SVector{2, Float64}}
end

function get_jacobian_inverse(k::SVector{2, Float64}, f::Function)
    J = Matrix{Float64}(undef, 2, 2)
    v = gradient(f, k)
    J[1,1] = v[1]
    J[1,2] = v[2]
    J[2,1] = - v[2] / norm(v)
    J[2,2] = v[1] / norm(v)
    return inv(J)
end

function temperature_broaden(fs::FermiSurface, f::Function, T::Real, n::Int, threshold::Real = 0.001)
    e_max::Float64 = 2 * T * acosh(1 / (2 * sqrt(threshold)))
    iseven(n) && (n += 1) # Enforce that there is a center point on the Fermi surface

    energies = LinRange(-e_max, e_max, n + 1)
    Δε = 2 * e_max / n
    mesh = Matrix{Patch}(undef, n, 0)
    mesh_corners = Matrix{SVector{2,Float64}}(undef, n + 1, 0)

    for segment in fs.segments
        ℓ = length(segment.points)
        vel = gradient.(f, segment.points)

        corners = Matrix{SVector{2, Float64}}(undef, n + 1, ℓ)
        for j in eachindex(segment.points)
            for i in eachindex(energies)
                corners[i, j] = get_k_bound(f, energies[i], segment.points[j], vel[j])
            end
        end
        
        m = segment.isclosed ? ℓ : ℓ - 1
        segment_mesh = Matrix{Patch}(undef, n, m)
        for i in 1:n
            for j in 1:ℓ-1
                patch_k = (corners[i, j] + corners[i + 1, j] + corners[i + 1, j + 1] + corners[i, j + 1]) / 4
                segment_mesh[i, j] = Patch(
                    patch_k, 
                    f,
                    get_patch_area(corners, i, j),
                    # [(i-1)*m + j, i*m + j, i*m + (j+1), (i-1)*m + (j+1)]
                    [(j-1)*(n+1) + i, (j-1)*(n+1) + i+1, j*(n+1) + i+1, j*(n+1) + i]
                    )
            end

            # Extra mesh point to close Fermi surface if closed
            if segment.isclosed
                patch_k = (corners[i, ℓ] + corners[i + 1, ℓ] + corners[i + 1, 1] + corners[i, 1]) / 4
                segment_mesh[i, ℓ] = Patch(
                    patch_k,
                    f,
                    get_patch_area(corners, i, ℓ),
                    [(ℓ-1)*(n+1) + i, (ℓ-1)*(n+1) + i+1, 1*(n+1) + i+1, 1*(n+1) + i]
                    )
            end
        end

        mesh = hcat(mesh, segment_mesh)
        mesh_corners = hcat(mesh_corners, corners)
    end

    return Mesh(vec(mesh), vec(mesh_corners)), Δε
end

function get_patch_area(A::Matrix, i::Int, j::Int)
    m = size(A)[2]
    if j != m
        α = A[i, j + 1] - A[i, j]
        β = A[i + 1, j + 1] - A[i, j + 1]
        γ = A[i + 1, j] - A[i + 1, j + 1]
    else
        α = A[i, 1] - A[i, j]
        β = A[i + 1, 1] - A[i, 1]
        γ = A[i + 1, j] - A[i + 1, 1]
    end
    δ = A[i, j] - A[i + 1, j]

    tri_area_1 = abs(α[1]*δ[2] - α[2]*δ[1]) / 2
    tri_area_2 = abs(β[1]*γ[2] - β[2]*γ[1]) / 2

    return tri_area_1 + tri_area_2
end
module FSGenerator

    export f0

    """
        f0(E, T)
    
    Return the value of the Fermi-Dirac distribution for energy E and temperature T.
    """
    f0(E::Float64, T::Float64) = 1 / (exp(E/T) + 1)

    include("./mesh/mesh.jl")

    include("./mesh/marching_squares.jl")

end # module FSGenerator

### Gamma Band of Sr2RuO4###

# In units of t0 = 119 meV
const mu::Float64 = 1.48 
const tx::Float64 = 1.0 
const ty::Float64 = 1.0 
const tp::Float64 = 0.392 
const ef::Float64 = 0.65 * (2.0*tx + 2.0*ty + 4*tp + mu) * 0.119
const Tf = (ef / 8.617333e-5) # K
# const λchar =  2pi / (6.582119569e−16) * 0.65 * ef * 0.119 * 1e-12 # Characteristic rate in ps^-1

efn = (2.0*tx + 2.0*ty + 4*tp + mu)
hamiltonian(k::AbstractVector) = (- 2.0 * tx * cos(2pi*k[1]) - 2.0 * ty * cos(2pi * k[2]) - 4 * tp * cos(2pi * k[1]) * cos(2pi * k[2]) - mu) / efn

### Deformation Potentials ###
const alph::Float64 = 7.604
const alph_p::Float64 = 7.604

# function dxx(k::SVector{2,Float64}, v::SVector{2,Float64})
#     return 2 * alph * (tx/ef) * cos(2pi*k[1]) + 2 * alph_p * (tp/ef) * cos(2pi*k[1]) * cos(2pi*k[2]) - k[1] * v[1] / ef
# end

# function dyy(k::SVector{2,Float64}, v::SVector{2,Float64})
#     return 2 * alph * (ty/ef) * cos(2pi*k[2]) + 2 * alph_p * (tp/ef) * cos(2pi*k[1]) * cos(2pi*k[2]) - k[2] * v[2] / ef
# end

# function dxy(k::SVector{2,Float64}, v::SVector{2,Float64})
#     return - 2 * alph_p * (tp/ef) * sin(2pi*k[1]) * sin(2pi*k[2]) - k[1] * v[2] / ef
# end
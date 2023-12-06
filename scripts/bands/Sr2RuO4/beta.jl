### Beta Band of Sr2RuO4###

# In units of t2 = 165 meV
const μ::Float64  = 1.08 
const t2::Float64 = 1.00
const t3::Float64 = 0.08 
const t5::Float64 = 0.13
const ef::Float64 = 3.24 * 0.165
efn               = 3.24
# const ef::Float64 = 2.0*t2 + 2.0*ty + 4*tp + mu # Without cyclotron mass correction
const Tf = (ef * 0.60 / 8.617333e-5)# K

exz(k) = (-2.0 * t2 * cos(2pi*k[1]) - 2.0 * t3 * cos(2pi*k[2]) - μ) #/ ef
eyz(k) = (-2.0 * t3 * cos(2pi*k[1]) - 2.0 * t2 * cos(2pi*k[2]) - μ) #/ ef
V(k)   = (4.0 * t5 * sin(2pi*k[1]) * sin(2pi*k[2])) / ef
hamiltonian(k) = 0.5 * ( (exz(k) + eyz(k)) + sqrt((exz(k) + eyz(k))^2 - 4 * (exz(k)*eyz(k) -  V(k)^2)) ) / efn


### Deformation Potentials ###
# const alph::Float64 = 7.604
# const alph_p::Float64 = 7.604

# function dxx(k::SVector{2,Float64}, v::SVector{2,Float64})
#     return 2 * alph * (tx/ef) * cos(2pi*k[1]) + 2 * alph_p * (tp/ef) * cos(2pi*k[1]) * cos(2pi*k[2]) - k[1] * v[1] / ef
# end

# function dyy(k::SVector{2,Float64}, v::SVector{2,Float64})
#     return 2 * alph * (ty/ef) * cos(2pi*k[2]) + 2 * alph_p * (tp/ef) * cos(2pi*k[1]) * cos(2pi*k[2]) - k[2] * v[2] / ef
# end

# function dxy(k::SVector{2,Float64}, v::SVector{2,Float64})
#     return - 2 * alph_p * (tp/ef) * sin(2pi*k[1]) * sin(2pi*k[2]) - k[1] * v[2] / ef
# end
### Example Band ###
# Near diamond FS 

const mu::Float64 = 0.1
const tx::Float64 = 1.0
const ty::Float64 = 1.0
const ef::Float64 = - 2.0 * tx - 2.0 * ty - mu
const Tf::Float64 = 4000

hamiltonian(k) = (- 2.0 * tx * cos(k[1]*2pi) - 2.0 * ty * cos(k[2]*2pi) - mu) / ef
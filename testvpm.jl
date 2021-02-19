using HDF5
using GeometricTools
using LinearAlgebra
using FLOWVPM
using LinearAlgebra

vpm=FLOWVPM;
gt = GeometricTools;

UJ = vpm.UJ_direct
zeta = vpm.zeta_direct

h = 10.0                  # Length to probe
nprobes = 20              # Number of probes per direction
D = rand(3)               # Random direction
D ./= norm(D)
X = rand(3)               # Position of particle
Gamma = (1/rand())*cross([1,0,0], D) # Vortex strength
sigma = range(0.05, 0.15, length=10)[ceil(Int, rand()*10)]*h        # Smoothing radius

# Create fluid domain grid
fdom = gt.Grid([-1,-1,-1],[1,1,1],[6,6,6])

# Initialize field
pfield = vpm.ParticleField(fdom.nnodes + 1)
# Add probes
for ni in 1:fdom.nnodes
    Xprobe = gt.get_node(fdom, ni)
    vpm.add_particle(pfield, Xprobe, 1e-10*ones(3), h/100)
end

    # Add vortex
vpm.add_particle(pfield, X, Gamma, sigma)

    # Calculate 𝝎 = ∇×𝐮
UJ(pfield)
Us = [vpm.get_U(P) for P in vpm.iterate(pfield)][1:end-1]
Ws = [vpm.get_W(P) for P in vpm.iterate(pfield)][1:end-1]

    # Calculate 𝝎 = ∑𝚪𝑝 𝜁𝜎(𝐱−𝐱𝑝)
zeta(pfield)

    # Save fluid domain
gt.add_field(fdom, "U", "vector", Us, "node")
gt.add_field(fdom, "W", "vector", Ws, "node")
gt.save(fdom,"/media/flowlab/Storage/gdevenport/simulations/pfield/testVTK")    
vpm.save(pfield,"/media/flowlab/Storage/gdevenport/simulations/pfield/pfield_save_test")
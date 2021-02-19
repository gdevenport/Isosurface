using HDF5
using GeometricTools
using LinearAlgebra
using FLOWVPM
using LinearAlgebra

vpm=FLOWVPM;
gt = GeometricTools;

UJ = vpm.UJ_direct
zeta = vpm.zeta_direct


# Create fluid domain grid
fdom = gt.Grid([-5,-5,-5],[5,5,5],[11,11,11])

# Read in the particle field
file = h5open("/media/flowlab/Storage/gdevenport/simulations/turbine_validation20/sim_pfield.5.h5","r");

# We now extract the various groups, Gamma, X, sigma, etc.
group_gamma = file["Gamma"];
group_X = file["X"];
group_Sigma = file["sigma"];

# We now extract the data from the various groups.
Gamma = read(group_gamma);
X = read(group_X);
Sigma = read(group_Sigma);

# Calculate the length of the data set
lengthGamma = size(Gamma)[2];
lengthX = size(X)[2];
lengthSigma = size(Sigma)[1];
println(lengthGamma,"\t\t", lengthX, "\t\t", lengthSigma)

Xpoints = X[1,:];
Ypoints = X[2,:];
Zpoints = X[3,:];
# Initialize field
pfield = vpm.ParticleField(fdom.nnodes + 1 + lengthX)

# Add probes to the particle field at the nodes of the grid.
for i in 1:fdom.nnodes
    Xprobe = gt.get_node(fdom, i)
end

# Now add the particles to the particle field with the probes, based on the data read in.
for i in 1:lengthX
    vpm.add_particle(pfield, [Xpoints[i],Ypoints[i],Zpoints[i]], Gamma[i], Sigma[i])
end

    # Calculate 𝝎 = ∇×𝐮
UJ(pfield)
Us = [vpm.get_U(P) for P in vpm.iterate(pfield)][1:fdom.nnodes]
Ws = [vpm.get_W(P) for P in vpm.iterate(pfield)][1:fdom.nnodes]

zeta(pfield)
omegaapproxs = [[P.Jexa[1], P.Jexa[2], P.Jexa[3]] for P in vpm.iterate(pfield)][1:fdom.nnodes]


    # Calculate 𝝎 = ∑𝚪𝑝 𝜁𝜎(𝐱−𝐱𝑝)

    # Save fluid domain
    
gt.add_field(fdom, "U", "vector", Us, "node")
gt.add_field(fdom, "W", "vector", Ws, "node")
gt.add_field(fdom, "omegaapproxs", "vector", omegaapproxs, "node")

gt.save(fdom,"/media/flowlab/Storage/gdevenport/simulations/pfield/testVTK")
vpm.save(pfield,"/media/flowlab/Storage/gdevenport/simulations/pfield/pfield_save_test")

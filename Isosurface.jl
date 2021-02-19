# Isosurface Creation
# Author: Greg Devenport
# Date: February 19, 2021
# This code creates the file necessary for isosurface creation in Para View.


using HDF5
using GeometricTools
using LinearAlgebra
using FLOWVPM

vpm=FLOWVPM;
gt = GeometricTools;
UJ = vpm.UJ_direct
zeta = vpm.zeta_direct

# Define the path and file name of the desired h5 file. 
path = "/media/flowlab/Storage/gdevenport/simulations/turbine_validation20/";
name = "sim_pfield.5.h5";

# Create fluid domain grid. Grid([x,y,z lower bounds],[x,y,z upper bounds],[nx,ny,nz number of divisions for each coordinate])
# The number of nodes for this grid will be (nx+1)*(ny+1)*(nz+1)
fdom = gt.Grid([-5,-5,-5],[5,5,5],[11,11,11])

"""
    readh5(file_name, file_path)

    Reads in an h5 file and extracts the Gamma, X (x,y,z position), and Sigma data. 
    Returns the extracted data and the number of particles in the h5 file. 

    return X, Gamma, Sigma, lengthX
"""
function readh5(file_name, file_path)
    file = h5open("$file_path$file_name","r");

    # We now extract the various groups, Gamma, X, sigma and extract the data from the various groups.
    # file["Gamma"] extracts the group Gamma. The read command reads in the data as an array. 
    Gamma = read(file["Gamma"]);
    X = read(file["X"]);
    Sigma = read(file["sigma"]);

    # Calculate the length of the data set 
    lengthX = size(X)[2];

    return X, Gamma, Sigma, lengthX
end

# Initialize field
pfield = vpm.ParticleField(fdom.nnodes + lengthX + 1)

# Add probes to the particle field at the nodes of the grid.
# add_particle(pfield, [x,y,z], Gamma, Sigma)
for i in 1:fdom.nnodes
    Xprobe = gt.get_node(fdom, i)
    vpm.add_particle(pfield, Xprobe, 1e-10*ones(3), 0.01)
end

# Now add the particles to the particle field with the probes, based on the data read in from the h5 file. 
# It should be noted that the [x,y,z] data are arranged with x,y,z as rows and the various particles as columns.
for i in 1:lengthX
    vpm.add_particle(pfield, X[:,i], Gamma[i], Sigma[i])
end

    # Calculate 𝝎 = ∇×𝐮
UJ(pfield)
Us = [vpm.get_U(P) for P in vpm.iterate(pfield)][1:fdom.nnodes]
Ws = [vpm.get_W(P) for P in vpm.iterate(pfield)][1:fdom.nnodes]

zeta(pfield)
omegaapproxs = [[P.Jexa[1], P.Jexa[2], P.Jexa[3]] for P in vpm.iterate(pfield)][1:fdom.nnodes]

# Add the velocity (Us) and vorticity (Ws) data to the fluid domain. 
gt.add_field(fdom, "U", "vector", Us, "node")
gt.add_field(fdom, "W", "vector", Ws, "node")
gt.add_field(fdom, "omegaapproxs", "vector", omegaapproxs, "node")

# Save the grid as a VTK file. 
gt.save(fdom,"/media/flowlab/Storage/gdevenport/simulations/pfield/testVTK")
vpm.save(pfield,"/media/flowlab/Storage/gdevenport/simulations/pfield/pfield_save_test")

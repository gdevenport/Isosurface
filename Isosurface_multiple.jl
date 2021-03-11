# Isosurface Creation
# Author: Greg Devenport
# Date: February 19, 2021
# This code creates the file necessary for isosurface creation in Para View.
using HDF5
using GeometricTools
using LinearAlgebra
using FLOWVPM

vpm = FLOWVPM;
gt = GeometricTools;
UJ = vpm.UJ_direct
zeta = vpm.zeta_direct

# Define the path and file name of the desired h5 file.
# path = "/media/flowlab/Storage/gdevenport/simulations/turbine_validation20_copy/";
# path = "/media/flowlab/Storage/jmehr/simulations/straight_path_test/";
path = "/media/flowlab/Storage/gdevenport/simulations/single_rotor_output/";
# File name is provided by the input for the function.

# Define the pfield path and file name for the saved file.
p_save_path = "/media/flowlab/Storage/gdevenport/simulations/pfield/";
p_save_name = "sim_pfield_test"

# Define the vtk path and file name for the saved file.
vtk_save_path = "/media/flowlab/Storage/gdevenport/simulations/vtk/";
vtk_save_name = "straight_path";


"""
iso_vtk(name::String, iteration::Int)
Inputs are 'name' which is the name of the h5 file, and 'iteration' which is an integer
representing the number associated with the h5 file. See get_file_number(). 
"""
function iso_vtk(name::String, iteration::Int)

# Create fluid domain grid. Grid([x,y,z lower bounds],[x,y,z upper bounds],[nx,ny,nz number of divisions for each coordinate])
# The number of nodes for this grid will be (nx+1)*(ny+1)*(nz+1)
fdom = gt.Grid([-5,-5,-5],[0.5,5,5],[11,11,11])

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
    gamma = read(file["Gamma"]);
    position = read(file["X"]);
    sigma = read(file["sigma"]);

    # Calculate the length of the data set 
    len_data = size(position)[2];

    return position, gamma, sigma, len_data
end

X, Gamma, Sigma , lengthX = readh5(name, path);

# Initialize particle field. We will add the test probes and the particles from the h5 file.
pfield = vpm.ParticleField(fdom.nnodes + lengthX + 1);

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
gt.save(fdom,"$vtk_save_path.$vtk_save_name$iteration")

# The variable 'iteration' pairs the original h5 file number with this output vtk file number so that
# the output files are kept in the same order as the input files. 
vpm.save(pfield, "$p_save_name.$iteration"; path = p_save_path)

end;


"""
get_file_number(testfile::String)

Input is a file name.
Coded for an h5 file name input. Extracts the h5 file number. 
Using sim_pfield.0012.h5 as the file name would return 0012
"""
function get_file_number(testfile::String)
file_number = 0;

# Iteration starts at end-3, this skips the '.h5' portion of the filename. 
i=3; 
while testfile[end-i] != '.'
    # Since i starts at 3, we want j to start at 0. 
    j = i - 3;
    file_number = file_number + 10^j * parse(Int, testfile[end-i]);
    i += 1;
end
return file_number
end;


function whole_folder()
filelist = readdir(path);

for file in filelist
    # Skip all files not ending in '5' this is designed to only use '.h5' files. 
    if file[end]=='5'
        iso_vtk(file, get_file_number(file))
        number = get_file_number(file);
        println("File number $number")
    end
end
end

function part_folder()
print("Enter file start: ");
file_start = parse(Int, readline());
print("Enter file end: ");
file_end = parse(Int,readline());
# file_name = "sim_pfield"
file_name = "tutorial_pfield"


for i in file_start:file_end
    iso_vtk("$file_name.$i.h5",i)
    println("Creating Isosurface for file $i")
end
end

part_folder()

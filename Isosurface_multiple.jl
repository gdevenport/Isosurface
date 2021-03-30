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


"""
    create_iso( h5file_name, iteration, bounds, freestream, data_path, save_path, vtk_save_name, pfield_save_name)
    Inputs are
    'h5file_name' which is the name of the h5 file, 
    'iteration' which is an integer representing the number associated with the h5 file. See get_file_number(), 
    'bounds' which is an array containing the bounds for the fluid domain grid,
    'freestream' which is a vector containing the freestream components used in the simulation,
    'data_path' which is a string where the h5 files are contained,
    'save_path' which is a string where you want to store the pfield and/or vtk files,
    'vtk_save_name' which is a string of the vtk file save names,
    'pfield_save_name' which is a string of the pfield file save names,
    'verbose' which is a bool, setting to true will cause many lines of text to be printed as you monitor the code progress. 
"""
function create_iso(
    h5file_name, 
    iteration,
    bounds,
    freestream,
    data_path,
    save_path,
    vtk_save_name,
    pfield_save_name,
    verbose
    )

# --------------------------------------------------------Create Fluid Domain------------------------------------------------------
    # Create fluid domain grid. Grid([x,y,z lower bounds],[x,y,z upper bounds],[nx,ny,nz number of divisions for each coordinate])
    # The number of nodes for this grid will be (nx+1)*(ny+1)*(nz+1)
    fdom = gt.Grid(bounds[1],bounds[2],convert(Array{Int64,1},bounds[3]))

    # Print file number code is running on. 
    if verbose
        println("Creating Isosurface for file $iteration")
    end

    X, Gamma, Sigma , lengthX = readh5(h5file_name, data_path);

    if verbose
        println("Building particle field...number of particles: $lengthX")
    end

# --------------------------------------------------------Create Particle Field----------------------------------------------------
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
    # The same is true of the Gamma data. 
    for i in 1:lengthX
        vpm.add_particle(pfield, X[:,i], Gamma[:,i], Sigma[i])
    end
 
# --------------------------------------------------------Calculate Vorticity and Velocity------------------------------------------
       #The pfield must be reset each iteration so that the velocities do not continue to add on top of eachother. 
    if verbose
        println("Calculating vorticity and velocity...")
        println("\t Resetting particle field...")
    end
        

    vpm._reset_particles(pfield)

    if verbose
        println("\t Calculating particle on particle interations...")
    end

    UJ(pfield)

    if verbose
        println("\t Calculating velocity...")
    end

    Us = [vpm.get_U(P)+freestream for P in vpm.iterate(pfield)][1:fdom.nnodes]

    if verbose
        println("\t Calculating vorticity...\n")
    end

    Ws = [vpm.get_W(P) for P in vpm.iterate(pfield)][1:fdom.nnodes]
    
    #zeta(pfield)
    #omegaapproxs = [[P.Jexa[1], P.Jexa[2], P.Jexa[3]] for P in vpm.iterate(pfield)][1:fdom.nnodes]
    

# --------------------------------------------------------Add Solutions to Fluid Domain--------------------------------------------- 
    # Add the velocity (Us) and vorticity (Ws) data to the fluid domain.
    gt.add_field(fdom, "U", "vector", Us, "node")
    gt.add_field(fdom, "W", "vector", Ws, "node")
    #gt.add_field(fdom, "omegaapproxs", "vector", omegaapproxs, "node")
    
    # Generate the file number to match the input h5 file. Output in .%4d format (0001, 0010, 0100, 1000).
    if iteration < 10
        file_number = "000$iteration";
    elseif iteration < 100
        file_number = "00$iteration";
    elseif iteration < 1000
        file_number = "0$iteration";
    else
        file_number = "$iteration";
    end


# --------------------------------------------------------Save VTK Files-----------------------------------------------------------
    # Save the grid as a VTK file.
    gt.save(fdom,"$save_path$vtk_save_name";num=iteration)
    
    # The variable 'iteration' pairs the original h5 file number with this output vtk file number so that
    # the output files are kept in the same order as the input files.
    # vpm.save(pfield, "$pfield_save_name"; path = save_path,num=iteration) #Turn on if you want pfield to save. 
    
end;


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


function iterate_iso(;
    file_start=file_start, 
    file_end=file_end,
    bounds=bounds,
    freestream=freestream,
    data_path=data_path,
    pfield_file_name=pfield_file_name,
    save_path=save_path,
    vtk_save_name=vtk_save_name,
    pfield_save_name=pfield_save_name,
    verbose = verbose
    )

    for i in file_start:file_end
        create_iso("$pfield_file_name.$i.h5",i,bounds,freestream,data_path,save_path,vtk_save_name,pfield_save_name,verbose)
    end

end

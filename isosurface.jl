# Isosurface Creation
# Author: Greg Devenport
# Date: February 19, 2021
# This code creates the files necessary for isosurface creation in Para View.
using HDF5
using GeometricTools
using LinearAlgebra
using FLOWVPM

vpm = FLOWVPM;
gt = GeometricTools;
UJ = vpm.UJ_direct


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
    'verbose' which is a bool, setting to true will cause many lines of text to be printed as you monitor the code progress. 
"""
function create_iso(;
    file_start=file_start,
    file_end=file_end,
    freestream=freestream,
    data_path=data_path,
    pfield_file_name=pfield_file_name,
    save_path=save_path,
    vtk_save_name=vtk_save_name,
    verbose=verbose,
    circular=circular,
    center=center,
    dimensions=dimensions,
    v=v,
    divisions=divisions,
    t_total=t_total,
    rotation_center=rotation_center
    )

    x_length = dimensions[1];
    y_length = dimensions[2];
    z_length = dimensions[3];

    # Change this if circlular path is in a different plane. 
    r = sqrt((center[2]-rotation_center[2])^2 + (center[3]-rotation_center[3])^2);

    n_steps = file_end - file_start;


    @time begin
        for i in file_start:file_end

            # --------------------------------------------------------Create Fluid Domain------------------------------------------------------
            # Create fluid domain grid. Grid([x,y,z lower bounds],[x,y,z upper bounds],[nx,ny,nz number of divisions for each coordinate])
            # The number of nodes for this grid will be (nx+1)*(ny+1)*(nz+1)
            
            if circular

                t = i*(t_total/n_steps) # This is the real time. 

                # So far this only works for a circular path in a plane (x/y, x/z, y/z) and not in three dimensions. 
                # Change x1, x2, y1, y2 to the appropriate variables so the circular path is in the desired plane. 
                z1 = rotation_center[3] + r*cos(t*v/r) - z_length/2;
                z2 = rotation_center[3] + r*cos(t*v/r) + z_length/2;
                y1 = rotation_center[2] + r*sin(t*v/r) + y_length/2;
                y2 = rotation_center[2] + r*sin(t*v/r) - y_length/2;
                x1 = center[1]-x_length/2;
                x2 = center[1]+x_length/2;

                circle_path_coordinates = [[min(x1,x2),min(y1,y2),min(z1,z2)],[max(x1,x2),max(y1,y2),max(z1,z2)]]
                
                fdom = gt.Grid(circle_path_coordinates[1],circle_path_coordinates[2],convert(Array{Int64,1}, divisions))
            end

            if ! circular

                x1 = center[1]-x_length/2;
                x2 = center[1]+x_length/2;
                y1 = center[2]-y_length/2;
                y2 = center[2]+y_length/2;
                z1 = center[3]+z_length/2;
                z2 = center[3]-z_length/2;

                fdom = gt.Grid([min(x1,x2),min(y1,y2),min(z1,z2)],[max(x1,x2),max(y1,y2),max(z1,z2)],convert(Array{Int64,1}, divisions))
            end
            # Print file number code is running on. 
            if verbose println("Creating Isosurface for file $i") end

            X, Gamma, Sigma , lengthX = readh5("$pfield_file_name.$i.h5",data_path);

            if verbose println("Building particle field...number of particles: $lengthX") end

            # --------------------------------------------------------Create Particle Field----------------------------------------------------
            # Initialize particle field. We will add the test probes and the particles from the h5 file.
            pfield_from_h5_file = vpm.ParticleField(lengthX);
            pfield_for_fluid_domain = vpm.ParticleField(fdom.nnodes + 1);

            # Add probes to the particle field at the nodes of the grid.
            # add_particle(pfield, [x,y,z], Gamma, Sigma)
            for i in 1:fdom.nnodes
                Xprobe = gt.get_node(fdom, i)
                vpm.add_particle(pfield_for_fluid_domain, Xprobe, 1e-10*ones(3), 0.01)
            end

            # Now add the particles to the particle field with the probes, based on the data read in from the h5 file.
            # It should be noted that the [x,y,z] data are arranged with x,y,z as rows and the various particles as columns.
            # The same is true of the Gamma data. 
            for i in 1:lengthX
                vpm.add_particle(pfield_from_h5_file, X[:,i], Gamma[:,i], Sigma[i])
            end
        
            # --------------------------------------------------------Calculate Vorticity and Velocity------------------------------------------
            #The pfield must be reset each iteration so that the velocities do not continue to add on top of eachother. 
            if verbose println("Calculating vorticity and velocity..."); println("\t Resetting particle field...") end
                

            vpm._reset_particles(pfield_for_fluid_domain)
            vpm._reset_particles(pfield_from_h5_file)


            if verbose println("\t Calculating particle on particle interations...") end

            UJ(pfield_from_h5_file,pfield_for_fluid_domain)

            if verbose println("\t Calculating velocity...") end

            Us = [vpm.get_U(P)+freestream for P in vpm.iterate(pfield_for_fluid_domain)]

            if verbose println("\t Calculating vorticity...\n") end

            Ws = [vpm.get_W(P) for P in vpm.iterate(pfield_for_fluid_domain)]
            
            #zeta(pfield)
            #omegaapproxs = [[P.Jexa[1], P.Jexa[2], P.Jexa[3]] for P in vpm.iterate(pfield)][1:fdom.nnodes]
            

            # --------------------------------------------------------Add Solutions to Fluid Domain--------------------------------------------- 
            # Add the velocity (Us) and vorticity (Ws) data to the fluid domain.
            gt.add_field(fdom, "U", "vector", Us, "node")
            gt.add_field(fdom, "W", "vector", Ws, "node")
            #gt.add_field(fdom, "omegaapproxs", "vector", omegaapproxs, "node")
            
            # Generate the file number to match the input h5 file. Output in .%4d format (0001, 0010, 0100, 1000).
            if i < 10
                file_number = "000$i";
            elseif i < 100
                file_number = "00$i";
            elseif i < 1000
                file_number = "0$i";
            else
                file_number = "$i";
            end


            # --------------------------------------------------------Save VTK Files-----------------------------------------------------------
            # Save the grid as a VTK file.
            gt.save(fdom,"$save_path$vtk_save_name";num=i)
            
            # The variable 'iteration' pairs the original h5 file number with this output vtk file number so that
            # the output files are kept in the same order as the input files.
        end
    end
end


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


function local_isosurface()
    # File start and stop
    file_start=9;
    file_end=12;

    # Define bounds

    # Define freestream
    freestream = [-20.0,0.0,0.0];

    # Data read
    data_path="/media/sf_Virtual_Drive/Winter 2021/Wake Visualization/Simulation Visualization Files/";
    pfield_file_name = "greg_test_pfield"

    # Data write
    save_path = "/media/sf_Virtual_Drive/Winter 2021/Wake Visualization/";
    vtk_save_name = "delete_me";
    
    verbose = true;
    circular = true;

    # Length of each side of fluid grid
    dimensions = [2,2,2];

    # Number of divisions per dimension. [4,4,4] would be (4+1)=5 divisions per dimension, and (4+1)^3 probes. 
    divisions = [1,1,1];

    # Center of fluid domain grid. 
    center = [4.0,0.0,0.0];

    # Center about which craft flies around. 
    rotation_center = [0.0,0.0,0.0];

    # Velocity of windcraft flying in circle. 
    v = 10;

    # Total time simulation simulates. 
    t_total = 1.0;
        
    create_iso(;
        file_start=file_start,
        file_end=file_end,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        vtk_save_name=vtk_save_name,
        verbose = verbose,
        circular = circular,
        center=center,
        dimensions=dimensions,
        v=v,
        divisions=divisions,
        t_total=t_total,
        rotation_center=rotation_center
        ) 
end



function windcraft_pylon_circle()
    # File start and stop
    file_start=1;
    file_end=100;

    # Define freestream
    freestream = [0.0,0.0,0.0];

    # Data read
    data_path="/media/flowlab/Storage/jmehr/simulations/4Rotor_straight_path_test_withpylons/";
    pfield_file_name = "sim_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/isosurface/windcraft_pylons_circle/";
    vtk_save_name = "wind_4_pylon_circle";

    verbose = true;
    circular = true;

    # Length of each side of fluid grid
    dimensions = [4,4,4];

    # Number of divisions per dimension. [4,4,4] would be (4+1)=5 divisions per dimension, and (4+1)^3 probes. 
    divisions = [1,1,1];

    # Center of fluid domain grid at start of simulation. 
    center = [0.0, 11.0, 0.0];
    
 
    # Center about which craft flies around. 
    rotation_center = [0.0,0.0,-60.0];

    # Velocity of windcraft flying in circle. 
    v = 1;

    # Total time simulation simulates. 
    t_total = 8.0;

    create_iso(
        file_start=file_start,
        file_end=file_end,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        vtk_save_name=vtk_save_name,
        verbose = verbose,
        circular = circular,
        center=center,
        dimensions=dimensions,
        v=v,
        divisions=divisions,
        t_total=t_total,
        rotation_center=rotation_center
        )

end

windcraft_pylon_circle()

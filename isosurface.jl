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
    create_iso_circular(file_start, file_end, freestream, data_path, pfield_file_name, save_path, vtk_save_name, verbose, 
    circular, center, dimensions, v, divisions, t_total, rotation_center)

    Inputs are

    'file_start' an int representing the file number start.
    'file_end' an int representing the file number end.
    'freestream' which is a vector containing the freestream components used in the simulation,
    'data_path' which is a string where the h5 files are contained,
    'pfield_file_name' a string specifying the pfield name, usually "sim_pfield"
    'save_path' which is a string where you want to store the pfield and/or vtk files,
    'vtk_save_name' which is a string of the vtk file save names,
    'verbose' which is a bool, setting to true will cause many lines of text to be printed as you monitor the code progress.
    'circular' a bool set to true if the simulation involves a circular path.
    'center' a vector of the origin of the fluid domain.
    'dimensions' a vector specifying the x,y,z dimensions of the fluid domain.
    'v' is the velocity of the craft.
    'divisions' a vector specifying the number of divisions in the x,y,z planes. 
    't_total' the total time the simulation ran for.
    'rotation_center' a vector specifying the point around which the vehicle moves about. 
"""
function create_iso_circular(;
    file_start=file_start,
    file_end=file_end,
    freestream=freestream,
    data_path=data_path,
    pfield_file_name=pfield_file_name,
    save_path=save_path,
    vtk_save_name=vtk_save_name,
    verbose=verbose,
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
            
          

            t = i*(t_total/n_steps) # This is the real time. 

            # So far this only works for a circular path in a plane (x/y, x/z, y/z) and not in three dimensions. 
            # Change x1, x2, y1, y2 to the appropriate variables so the circular path is in the desired plane. 
            z1 = rotation_center[3] + r*cos(t*v/r) - z_length/2 + center[3];
            z2 = rotation_center[3] + r*cos(t*v/r) + z_length/2 + center[3];
            y1 = rotation_center[2] + r*sin(t*v/r) + y_length/2 + center[2];
            y2 = rotation_center[2] + r*sin(t*v/r) - y_length/2 + center[2];
            x1 = rotation_center[1] - x_length/2 + center[1];
            x2 = rotation_center[1] + x_length/2 + center[1];

            circle_path_coordinates = [[min(x1,x2),min(y1,y2),min(z1,z2)],[max(x1,x2),max(y1,y2),max(z1,z2)]]
                
            fdom = gt.Grid(circle_path_coordinates[1],circle_path_coordinates[2],convert(Array{Int64,1}, divisions))

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
    create_iso_stationary(file_start, file_end, freestream, data_path, pfield_file_name, save_path, vtk_save_name, verbose, 
    circular, center, dimensions, v, divisions, t_total, rotation_center)

    Inputs are
    
    'file_start' an int representing the file number start.
    'file_end' an int representing the file number end.
    'freestream' which is a vector containing the freestream components used in the simulation,
    'data_path' which is a string where the h5 files are contained,
    'pfield_file_name' a string specifying the pfield name, usually "sim_pfield"
    'save_path' which is a string where you want to store the pfield and/or vtk files,
    'vtk_save_name' which is a string of the vtk file save names,
    'verbose' which is a bool, setting to true will cause many lines of text to be printed as you monitor the code progress.
    'center' a vector of the origin of the fluid domain.
    'dimensions' a vector specifying the x,y,z dimensions of the fluid domain.
    'divisions' a vector specifying the number of divisions in the x,y,z planes. 
"""
function create_iso_stationary(;
    file_start=file_start,
    file_end=file_end,
    freestream=freestream,
    data_path=data_path,
    pfield_file_name=pfield_file_name,
    save_path=save_path,
    vtk_save_name=vtk_save_name,
    verbose=verbose,
    center=center,
    dimensions=dimensions,
    divisions=divisions,
    )

    x_length = dimensions[1];
    y_length = dimensions[2];
    z_length = dimensions[3];


    @time begin
        for i in file_start:file_end

            # --------------------------------------------------------Create Fluid Domain------------------------------------------------------
            # Create fluid domain grid. Grid([x,y,z lower bounds],[x,y,z upper bounds],[nx,ny,nz number of divisions for each coordinate])
            # The number of nodes for this grid will be (nx+1)*(ny+1)*(nz+1)
            

            x1 = center[1] - x_length/2;
            x2 = center[1] + x_length/2;
            y1 = center[2] - y_length/2;
            y2 = center[2] + y_length/2;
            z1 = center[3] + z_length/2;
            z2 = center[3] - z_length/2;

            fdom = gt.Grid([min(x1,x2),min(y1,y2),min(z1,z2)],[max(x1,x2),max(y1,y2),max(z1,z2)],convert(Array{Int64,1}, divisions))
        

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
include("Isosurface_multiple.jl")

function local_isosurface()
    # File start and stop
    file_start=9;
    file_end=12;

    # Define bounds
    bounds=[[1.0,1.0,1.0],[2.0,2.0,2.0],[3,3,3]];

    # Define freestream
    freestream = [-20.0,0.0,0.0];

    # Data read
    data_path="/media/sf_Virtual_Drive/Winter 2021/Wake Visualization/Simulation Visualization Files/";
    pfield_file_name = "greg_test_pfield"

    # Data write
    save_path = "/media/sf_Virtual_Drive/Winter 2021/Wake Visualization/";
    pfield_save_name = "delete_me";
    vtk_save_name = "delete_me";
        
    create_isosurface(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name
        ) 

end

function turbine_validation_isosurface()

    # File start and stop
    file_start=1;
    file_end=10;

    # Define bounds
    
    # bounds=[[-10.0,-0.5,-5.0],[0.5,0.5,5.0],[21,3,21]];

    bounds=[[-10.0,-5.0,-5.0],[0.5,5.0,5.0],[17,17,17]];


    # Define freestream
    freestream = [-20.0,0.0,0.0];

    # Data read
    data_path="/media/flowlab/Storage/gdevenport/simulations/turbine_validation20_copy/";
    pfield_file_name = "sim_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/vtk/turbine_output/";
    pfield_save_name = "turbine_validation20";
    vtk_save_name = "turbine_validation20";
        
    create_isosurface(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name
        ) 

end

function straight_path_isosurface()

    # File start and stop
    file_start=1;
    file_end=10;

    # Define bounds
    bounds=[[-1.0,-1.0,-1.0],[2.0,1.0,1.0],[20,5,20]];

    # Define freestream
    freestream = [-20.0,0.0,0.0];

    # Data read
    data_path="/media/flowlab/Storage/jmehr/simulations/straight_path_test/";
    pfield_file_name = "sim_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/straight_path_output/";
    pfield_save_name = "straight_path";
    vtk_save_name = "straight_path";
        
    create_isosurface(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name
        ) 
end

function single_rotor_isosurface()

    # File start and stop
    file_start=1;
    file_end=80;

    # Define bounds
    bounds=[[-0.02,-0.001,-0.15],[0.2,0.001,0.15],[20,5,20]];
    #bounds=[[-0.02,-0.15,-0.15],[0.2,0.15,0.15],[17,17,17]];

    # Define freestream
    freestream = [0.0,0.0,0.0];

    # Data read
    data_path="/media/flowlab/Storage/gdevenport/simulations/single_rotor/";
    pfield_file_name = "singlerotor_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/single_rotor_output/";
    pfield_save_name = "single_rotor";
    vtk_save_name = "single_rotor";
        
    create_isosurface(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name
        ) 
end

function phantom_hover_isosurface()
    # File start and stop
    file_start=1;
    file_end=80;

    # Define bounds
    bounds=[[-0.1,-0.2,-0.2],[0.25,0.2,0.2],[20,5,20]];

    # Define freestream
    freestream = [0.0,0.0,0.0];

    # Data read
    data_path="/media/flowlab/Storage/gdevenport/simulations/phantom_hover/singlerotor_hover_test00/";
    pfield_file_name = "singlerotor_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/Vorticity_validation/";
    pfield_save_name = "phantom";
    vtk_save_name = "phantom";
        
    create_isosurface(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name
        ) 

end

function phantom_flight_isosurface()
    # File start and stop
    file_start=1;
    file_end=80;

    # Define bounds
    bounds=[[-0.1,-0.2,-0.2],[0.25,0.2,0.2],[20,5,20]];
    # Define freestream
    freestream = [0.0,0.0,0.0];

    # Data read
    data_path="/media/flowlab/Storage/gdevenport/simulations/data/phantom/singlerotor_fflight_test01/";
    pfield_file_name = "singlerotor_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/isosurface/phantom flight/";
    pfield_save_name = "phantom";
    vtk_save_name = "phantom";
        
    create_isosurface(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name
        ) 

end

phantom_flight_isosurface()
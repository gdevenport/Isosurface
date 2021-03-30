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
    
    verbose = true;
        
    iterate_iso(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name,
        verbose = verbose
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
    data_path="/media/flowlab/Storage/gdevenport/simulations/data/turbine_validation20/";
    pfield_file_name = "sim_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/isosurface/turbine validation/";
    pfield_save_name = "turbine_validation20";
    vtk_save_name = "turbine_validation20";

    verbose = true;
        
    iterate_iso(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name,
        verbose = true
        ) 

end

function straight_path_isosurface()

    # File start and stop
    file_start=1;
    file_end=100;

    # Define bounds
    bounds=[[-1.0,-1.0,-1.0],[2.0,1.0,1.0],[20,5,20]];

    # Define freestream
    freestream = [-20.0,0.0,0.0];

    # Data read
    data_path="/media/flowlab/Storage/jmehr/simulations/straight_path_test/";
    pfield_file_name = "sim_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/isosurface/straight path/";
    pfield_save_name = "straight_path";
    vtk_save_name = "straight_path";

    verbose = true;
        
    iterate_iso(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name,
        verbose = verbose
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
    data_path="/media/flowlab/Storage/gdevenport/simulations/data/single_rotor/";
    pfield_file_name = "singlerotor_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/isosurface/single rotor/";
    pfield_save_name = "single_rotor";
    vtk_save_name = "single_rotor";

    verbose = true;
        
    iterate_iso(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name,
        verbose = verbose
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
    data_path="/media/flowlab/Storage/gdevenport/simulations/data/phantom/singlerotor_hover_test00/";
    pfield_file_name = "singlerotor_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/isosurface/phantom hover/";
    pfield_save_name = "phantom";
    vtk_save_name = "phantom";

    verbose = true;
        
    iterate_iso(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name,
        verbose = verbose
        ) 

end

function phantom_flight_isosurface()
    # File start and stop
    file_start=1;
    file_end=80;

    # Define bounds
    bounds=[[-0.1,-0.01,-0.2],[0.25,0.01,0.2],[20,5,20]];
    # Define freestream
    freestream = [0.0,0.0,0.0];

    # Data read
    data_path="/media/flowlab/Storage/gdevenport/simulations/data/phantom/singlerotor_fflight_test01/";
    pfield_file_name = "singlerotor_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/isosurface/phantom flight/";
    pfield_save_name = "phantom";
    vtk_save_name = "phantom";
        
    verbose = true;

    iterate_iso(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name,
        verbose = verbose
        ) 

end

function turbine_24()
    # File start and stop
    file_start=0;
    file_end=72;

    # Define bounds
    bounds=[[-20,-7,-7],[1,7,7],[30,30,30]];
    # Define freestream
    freestream = [0.0,0.0,0.0];

    # Data read
    data_path="/media/flowlab/Storage/jmehr/simulations/oldform_turbinevalidation_FINAL_GOOD_DONOTDELETE_24/";
    pfield_file_name = "sim_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/isosurface/turbine 24/";
    pfield_save_name = "turbine_24";
    vtk_save_name = "turbine_24";
    
    verbose = true;

    iterate_iso(
        file_start=file_start,
        file_end=file_end,
        bounds=bounds,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        pfield_save_name=pfield_save_name,
        vtk_save_name=vtk_save_name,
        verbose = verbose
        ) 

end

turbine_24()

include("isosurface.jl")
include("read_h5.jl")

function windcraft_pylon()
    # File start and stop
    file_start=0;
    file_end=100;

    # Define freestream
    freestream = [0.0,0.0,0.0];

    # Data read
    data_path="/media/flowlab/Storage/jmehr/simulations/4Rotor_straight_path_test_withpylons/";
    pfield_file_name = "sim_pfield"

    # Data write
    save_path = "/media/flowlab/Storage/gdevenport/simulations/isosurface/windcraft_pylons_circle/";
    vtk_save_name = "wind_4_pylon";

    verbose = true;

    dimensions = [4.0,4.0,5.5]

    divisions = [10,10,10]

    center = [0.0,0.0,0.0]


    create_iso_stationary(
        file_start=file_start,
        file_end=file_end,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        vtk_save_name=vtk_save_name,
        verbose = verbose,
        center=center,
        dimensions=dimensions,
        divisions=divisions
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
    dimensions = [5,5,5];

    # Number of divisions per dimension. [4,4,4] would be (4+1)=5 divisions per dimension, and (4+1)^3 probes. 
    divisions = [10,10,10];

    # Center of fluid domain grid at start of simulation. 
    center = [0.0, 1.0, 0.0];
    
 
    # Center about which craft flies around. 
    rotation_center = [0.0,0.0,-60.0];

    # Velocity of windcraft flying in circle. 
    v_vehicle = 1;

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
        v_vehicle=v_vehicle,
        divisions=divisions,
        t_total=t_total,
        rotation_center=rotation_center
        )

end

function local_isosurface()
    # File start and stop
    file_start=9;
    file_end=12;

    # Define freestream
    freestream = [-20.0,0.0,0.0];

    # Data read
    data_path="/media/sf_Virtual_Drive/Winter 2021/Wake Visualization/Simulation Visualization Files/";
    pfield_file_name = "greg_test_pfield"

    # Data write
    save_path = "/media/sf_Virtual_Drive/Winter 2021/Wake Visualization/";
    vtk_save_name = "delete_me";

    verbose = true;

    # Length of each side of fluid grid
    dimensions = [2,2,2];

    # Number of divisions per dimension. [4,4,4] would be (4+1)=5 divisions per dimension, and (4+1)^3 probes. 
    divisions = [1,1,1];

    # Center of fluid domain grid. 
    center = [4.0,0.0,0.0];

    create_iso_stationary(;
        file_start=file_start,
        file_end=file_end,
        freestream=freestream,
        data_path=data_path,
        pfield_file_name=pfield_file_name,
        save_path=save_path,
        vtk_save_name=vtk_save_name,
        verbose = verbose,
        center=center,
        dimensions=dimensions,
        divisions=divisions,
        ) 
end

local_isosurface()
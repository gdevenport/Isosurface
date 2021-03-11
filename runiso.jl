include("Isosurface_multiple.jl")

# File start and stop
file_start=1;
file_end=10;

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
    data_path=data_path,
    pfield_file_name=pfield_file_name,
    save_path=save_path,
    pfield_save_name=pfield_save_name,
    vtk_save_name=vtk_save_name
    ) 

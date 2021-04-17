# Isosurface
This is my code that will read in h5 files and create vtk files to be read by paraview. This code calculates velocity and vorticity. Old pieces of code are found in **old code**. A recent update added the ability for the location of isosurface creation to change in time. This is helpful for the case of a moving aircraft, such as viewing the velocity field around the propellers. 

Here are a few samples of what the code can do. 
**Velocity slice showing velocity profile of wake behind turbine**
![alt-text](https://github.com/gdevenport/Isosurface/blob/main/media/front_side.gif)

**Velocity behind windcraft turbine blades**
![alt-text](https://github.com/gdevenport/Isosurface/blob/main/media/velocity_front.gif)

**Vorticity isosurfaces of a stationary propeller**
![alt-text](https://github.com/gdevenport/Isosurface/blob/main/media/tips_tubes_vorticity.gif)
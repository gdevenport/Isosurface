using FLOWUnsteady
uns = FLOWUnsteady
vlm = uns.vlm

span = 1.0              #wing span
aspectratio = 10.0      #wing aspect ratio
taperratio = 0.5        #wing taper ratio
wingtwist = 0.0         #wing twist
wingsweep = 10.0        #wing sweep in degrees
wingdihedral = 7.0      #wing dihedral in degrees

mainwing = vlm.simpleWing(span,aspectratio,taperratio,wingtwist,wingsweep,wingdihedral)

system = vlm.WingSystem()

vlm.addwing(system,"mainwing",mainwing)

Vinf(x,t) = [1,0,0]         #non-dimensional function defining free stream velocity
vlm.setVinf(system, Vinf)   #set freestream velocity for the system

run_name = "tutorial"           #define identifier at beginning of file names
#save_path = "/media/sf_Virtual_Drive/Winter 2021/Isosurface/single rotor output2.0"     #define directory where files will be saved
save_path = "/media/flowlab/Storage/gdevenport/simulations/single_rotor_output"     #define directory where files will be saved

run(`rm -rf $save_path`)        #clear out directory where files will be saved
run(`mkdir $save_path`)         #re-create directory fresh

vlm.save(system, run_name; path=save_path)  #save geometry in a .vtk file format

rotor_file = "apc10x7.csv"
data_path = uns.def_data_path

R, B = uns.read_rotor(rotor_file; data_path=data_path)[[1,3]] #get the radius for later

rotor = uns.generate_rotor(rotor_file; pitch=0.0,
                                            n=10, CW=true, ReD=1.5e6,
                                            verbose=false, xfoil=false, # changed xfoil here
                                            data_path=data_path,
                                            plot_disc=false);

rotors = vlm.Rotor[rotor];           

vehicleorigin = [0.0; 0.0; 0.0]
vehicleaxis = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]

rotororigin = [-0.1; 0.0; 0.0]

for rotor in rotors
    vlm.setcoordsystem(rotor, rotororigin, vehicleaxis; user=true)
end

rotor_systems = (rotors,);

for rotor in rotors; vlm.addwing(system, run_name, rotor); end;

RPMref = 60       #reference RPM
for rotor in rotors; vlm.setRPM(rotor, RPMref); end;

run(`rm -rf $save_path`)
run(`mkdir $save_path`)

vlm.setVinf(system, Vinf)   #set freestream velocity for the system

vlm.save(system, run_name; path=save_path)

vlm_system = vlm.WingSystem()

vlm.addwing(vlm_system, "mainwing", mainwing)

vlm.setVinf(vlm_system, Vinf)   #set freestream velocity for the system


wake_system = vlm.WingSystem()

vlm.addwing(wake_system, "SolveVLM", vlm_system)

vlm.setVinf(wake_system, Vinf)   #set freestream velocity for the system


for rotor in rotors; vlm.addwing(wake_system, run_name, rotor); end;

tilting_systems = ();

Vvehicle(t) = [-1.0,0.0,0.0]

anglevehicle(t) = zeros(3)

angle = ();

RPM_fun(t) = 1.0

RPM = (RPM_fun, );

maneuver = uns.KinematicManeuver(angle, RPM, Vvehicle, anglevehicle)

vehicle = uns.VLMVehicle(   system;
                            tilting_systems = tilting_systems,
                            rotor_systems   = rotor_systems,
                            vlm_system      = vlm_system,
                            wake_system     = wake_system,
                        );
Vref = 10.0         #define a reference velocity for the vehicle
ttot = 1.0          #define a total simulation time, in seconds
nsteps = 30        #define the number of steps the simulation will take
                        
                        #initial conditions
tinit = 0.0                                  #initial time
Vinit = Vref*maneuver.Vvehicle(tinit/ttot)   #initial linear velocity
Winit = zeros(3)                             #initial angular velocity

simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot; Vinit=Vinit, Winit=Winit, t=tinit);



nullfunc(args...) = false
    
pfield = uns.run_simulation(simulation, nsteps;
                        surf_sigma=R/10,
                        Vinf=Vinf,
                        save_path=save_path,
                        run_name=run_name,
                        prompt=true,
                        verbose=true,
                        extra_runtime_function=nullfunc
                                                )
                                
run(`paraview --data="$save_path/$(files);tutorial_pfield...vtk"`)
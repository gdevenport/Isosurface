# Load simulation engine
# import FLOWUnsteady
using FLOWUnsteady
uns = FLOWUnsteady
vlm = uns.vlm

import GeometricTools
gt = GeometricTools

# ------------ GLOBAL VARIABLES ------------------------------------------------
# Default path where to save data
extdrive_path = "/media/flowlab/Storage/gdevenport/simulations/single_rotor";
# extdrive_path = "temps/"



# ------------ DRIVERS ---------------------------------------------------------
function run_singlerotor_hover(; xfoil=true, prompt=true)

    J = 0.00                # Advance ratio Vinf/(nD)
    angle = 0.0             # (deg) angle of freestream (0 == climb, 90==forward flight)

    singlerotor(;   xfoil=xfoil,
                    VehicleType=uns.VLMVehicle,
                    J=J,
                    DVinf=[cos(pi/180*angle), sin(pi/180*angle), 0],
                    save_path=extdrive_path,
                    prompt=prompt)
end

# ------------------------------------------------------------------------------

function singlerotor(;  xfoil       = false,             # Whether to run XFOIL
                        VehicleType = uns.VLMVehicle,   # Vehicle type
                        J           = 0.0,              # Advance ratio
                        DVinf       = [1.0, 0, 0],      # Freestream direction
                        nrevs       = 6,                # Number of revolutions
                        nsteps_per_rev = 72,            # Time steps per revolution
                        # OUTPUT OPTIONS
                        save_path   = nothing,
                        run_name    = "singlerotor",
                        prompt      = true,
                        verbose     = true,
                        v_lvl       = 0)

    # TODO: Wake removal ?

    # ------------ PARAMETERS --------------------------------------------------

    # Rotor geometry
    rotor_file = "DJI-II.csv"           # Rotor geometry
    data_path = uns.def_data_path       # Path to rotor database
    pitch = 0.0                         # (deg) collective pitch of blades
    # n = 50                              # Number of blade elements
    n = 10
    CW = false                          # Clock-wise rotation
    # xfoil = false                     # Whether to run XFOIL

    # Read radius of this rotor and number of blades
    R, B = uns.read_rotor(rotor_file; data_path=data_path)[[1,3]]

    # Simulation parameters
    RPM = 81*60                         # RPM
    # J = 0.00001                       # Advance ratio Vinf/(nD)
    rho = 1.225                         # (kg/m^3) air density
    mu = 1.81e-5                        # (kg/ms) air dynamic viscosity
    ReD = 2*pi*RPM/60*R * rho/mu * 2*R  # Diameter-based Reynolds number

    magVinf = J*RPM/60*(2*R)
    Vinf(X,t) = magVinf*DVinf           # (m/s) freestream velocity

    # Solver parameters
    # nrevs = 6                         # Number of revolutions in simulation
    # nsteps_per_rev = 72                 # Time steps per revolution
    p_per_step = 2                      # Sheds per time step
    ttot = nrevs/(RPM/60)               # (s) total simulation time
    nsteps = nrevs*nsteps_per_rev       # Number of time steps
    lambda = 2.125                      # Core overlap
    overwrite_sigma = lambda * 2*pi*R/(nsteps_per_rev*p_per_step) # Smoothing core size
    surf_sigma = R/10                   # Smoothing radius of lifting surface
    vlm_sigma = surf_sigma              # Smoothing radius of VLM
    shed_unsteady = true                # Shed particles from unsteady loading

    max_particles = ((2*n+1)*B)*nrevs*nsteps_per_rev*p_per_step # Max particles for memory pre-allocation
    plot_disc = false                    # Plot blade discretization for debugging


    # ------------ SIMULATION SETUP --------------------------------------------
    # Generate rotor
    rotor = uns.generate_rotor(rotor_file; pitch=pitch,
                                            n=n, CW=CW, ReD=ReD,
                                            verbose=verbose, xfoil=xfoil,
                                            data_path=data_path,
                                            plot_disc=plot_disc)
    # ----- VEHICLE DEFINITION
    # System of all FLOWVLM objects
    system = vlm.WingSystem()
    vlm.addwing(system, run_name, rotor)

    # Systems of rotors
    rotors = vlm.Rotor[rotor]   # Defining this rotor as its own system
    rotor_systems = (rotors,)

    # Wake-shedding system (doesn't include the rotor if quasi-steady vehicle)
    wake_system = vlm.WingSystem()

    if VehicleType != uns.QVLMVehicle
        vlm.addwing(wake_system, run_name, rotor)
    else
        # Mute colinear warnings. This is needed since the quasi-steady solver
        #   will probe induced velocities at the lifting line of the blade
        uns.vlm.VLMSolver._mute_warning(true)
    end

    # FVS's Vehicle object
    vehicle = VehicleType(   system;
                                rotor_systems=rotor_systems,
                                wake_system=wake_system
                             )

    # ----- MANEUVER DEFINITION
    RPM_fun(t) = 1.0                # RPM (normalized by reference RPM) as a
                                    # function of normalized time

    angle = ()                      # Angle of each tilting system (none in this case)
    sysRPM = (RPM_fun, )              # RPM of each rotor system
    Vvehicle(t) = zeros(3)          # Translational velocity of vehicle over Vcruise
    anglevehicle(t) = zeros(3)      # (deg) angle of the vehicle

    # FVS's Maneuver object
    maneuver = uns.KinematicManeuver(angle, sysRPM, Vvehicle, anglevehicle)

    # Plot maneuver path and controls
    # uns.plot_maneuver(maneuver; vis_nsteps=nsteps)


    # ----- SIMULATION DEFINITION
    RPMref = RPM
    Vref = 0.0
    simulation = uns.Simulation(vehicle, maneuver, Vref, RPMref, ttot)

    monitor(args...) = false


    # ------------ RUN SIMULATION ----------------------------------------------
    pfield = uns.run_simulation(simulation, nsteps;
                                      # SIMULATION OPTIONS
                                      Vinf=Vinf,
                                      # SOLVERS OPTIONS
                                      p_per_step=p_per_step,
                                      overwrite_sigma=overwrite_sigma,
                                      vlm_sigma=vlm_sigma,
                                      surf_sigma=surf_sigma,
                                      max_particles=max_particles,
                                      shed_unsteady=shed_unsteady,
                                      extra_runtime_function=monitor,
                                      # OUTPUT OPTIONS
                                      save_path=save_path,
                                      run_name=run_name,
                                      prompt=prompt,
                                      verbose=verbose, v_lvl=v_lvl,
                                      )
    return pfield, rotor
end

singlerotor()

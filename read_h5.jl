# This code is to read in h5 files and view the data. 
using HDF5
using GeometricTools
using LinearAlgebra
using Plots
gt = GeometricTools;

# We first read in the h5 file. 
function readh5(filename)
    global file = h5open("/media/sf_Virtual_Drive/Winter 2021/Isosurface/Pfield_from_Judd/$filename","r");

    # We now extract the various groups, Gamma, X, sigma, etc. 
    group_gamma = file["Gamma"];
    group_X = file["X"];

    # We now extract the data from the various groups.
    global Gamma = read(group_gamma);
    global X = read(group_X);

    # Calculate the length of the data set 

    global lengthGamma = size(Gamma)[2];
    global lengthX = size(X)[2];

end

function printGamma() # I believe that Gamma should be a scalar value, so I'm not sure why the gamma data has three components...
    for i=1:lengthGamma
        if i==1
            println("\tGamma-X\t\t\tGamma-Y\t\t\tGamma-Z")
            println("---------------------------------------------------------------------")
        end
        println(Gamma[1,i],"\t", Gamma[2,i],"\t",Gamma[3,i])
    end
end

function printX()
    for i=1:lengthX
        if i==1
            println("\tX\t\t\tY\t\t\tZ")
            println("---------------------------------------------------------------------")
        end
        println(X[1,i],"\t", X[2,i],"\t",X[3,i])
    end
end

function writeX()
    io = open("writefile.txt","w")
    write(io, "This file is read only")
    close(io)
end


function normGamma() # I believe that Gamma should be a scalar value, so I'm not sure why the gamma data has three components...
    global normGammaList = [];
    for i=1:lengthGamma
        append!(normGammaList,norm(Gamma[:,i]))
    end
    for i=1:lengthGamma
        println(normGammaList[i])
    end
end

function makeGrid()
    fdom = gt.Grid([-1,-1,-1],[1,1,1],[6,6,6])
    newArray = normGammaList[1:343]
    gt.add_field(fdom,"Test","scalar",newArray,"node")
    gt.save(fdom,"testVTK")

end

function buildArray(n)
    array = [];
    for i=1:n
        append!(array,i)
    end
    return array;
end

file = h5open("/media/sf_Virtual_Drive/Winter 2021/Wake Visualization/simulations/greg_test_pfield.6.h5","r");

# We now extract the various groups, Gamma, X, sigma, etc. 
group_gamma = file["Gamma"];
group_X = file["X"];
group_sigma = file["sigma"];

# We now extract the data from the various groups.
Gamma = read(group_gamma);
X = read(group_X);
Sigma = read(group_sigma);

# Calculate the length of the data set 

lengthGamma = size(Gamma)[2];
lengthX = size(X)[2];
lengthSigma = size(Sigma)[1];
println(lengthGamma,"\t\t", lengthX, "\t\t", lengthSigma)

function plotstuff()
    readh5("sim_pfield.2.h5")
    xpoints=X[1,:];
    ypoints=X[2,:];
    zpoints=X[3,:];
    scatter(ypoints,zpoints,xpoints)
end


# Isosurface Creation
# Author: Greg Devenport
# Date: February 19, 2021
# This code extracts the needed data from an h5 file for use with isosurface.jl
using HDF5

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
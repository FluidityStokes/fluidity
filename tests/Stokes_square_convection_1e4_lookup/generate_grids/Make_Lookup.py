#!/usr/bin/python

# Written by Rhodri Davies, October 2011.
# Program to transform lookup table into a grid structure, thus allowing for ease of navigation:
# Note that this currently works when lookup table increments in pressure first then temperature.

# Import required functionality:

import sys
import numpy

# Check script input:
def usage() :
    msg="Usage: %prog lookup_file_1 tmin tmax tinc pmin (GPa) pmax (GPa) pinc (GPa)"
    return msg

if len(sys.argv) < 7:
    print usage()
    sys.exit(1)

# Command line arguments:
main_datafile        = sys.argv[1]
tmin                 = float(sys.argv[2])
tmax                 = float(sys.argv[3])
tinc                 = float(sys.argv[4])
pmin                 = float(sys.argv[5])
pmax                 = float(sys.argv[6])
pinc                 = float(sys.argv[7])

# Load data:
main_data=numpy.loadtxt(main_datafile,usecols=(0,1,2,3,4,5,6,7,8,9,10,11)) # Read in but skip final column for now

# Store columns of interest and create new text file:
final_data=numpy.zeros((main_data.shape[0])*4).reshape(main_data.shape[0],4)
final_data[:,0] = main_data[:,0]                # Pressure (GPa)
final_data[:,1] = main_data[:,2]                # Temperature (k)
final_data[:,2] = main_data[:,3]                # Density (g/cm^3) (1 g/cm^3 = 1000 kg/m^3)
final_data[:,3] = main_data[:,9]                # H: Enthalpy (kJ/Kg)

############### Now transform data into grid structure:   ########################
tpoints     = ((tmax - tmin) / tinc) + 1
ppoints     = ((pmax - pmin) / pinc) + 1
nattributes = 2

pressure_values    = numpy.linspace(pmin,pmax,ppoints)
temperature_values = numpy.linspace(tmin,tmax,tpoints)

grid_data   = numpy.zeros(ppoints*tpoints*nattributes).reshape(tpoints,ppoints,nattributes)

for i in range(grid_data.shape[0]):
    for j in range(grid_data.shape[1]):
        lookup_point = (ppoints*i) + j
        grid_data[i,j,0] = final_data[lookup_point,2] # Density
        grid_data[i,j,1] = final_data[lookup_point,3] # Enthalpy

numpy.save("Grid_Data",grid_data)



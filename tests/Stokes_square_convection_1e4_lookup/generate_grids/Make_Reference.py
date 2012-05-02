#!/usr/bin/python

import sys
import numpy

def usage() :
    msg="Usage: %prog reference_lookup_file"
    return msg

if len(sys.argv) < 1:
    print usage()
    sys.exit(1)

# Command line arguments:
main_datafile = sys.argv[1]

# Load data:
main_data = numpy.loadtxt(main_datafile)

# Store columns of interest and create new text file:
final_data=numpy.zeros((main_data.shape[0])*4).reshape(main_data.shape[0],4)
final_data[:,0] = main_data[:,1]*1000  # Depth (m)
final_data[:,1] = main_data[:,0]       # Reference Pressure (GPa)
final_data[:,2] = main_data[:,2]       # Reference Temperature (k)
final_data[:,3] = main_data[:,3]       # Reference Density (g/cm^3)

# Save text file:
numpy.savetxt("Lookup_Reference_Final.dat",final_data,fmt="%f")
numpy.save("Lookup_Reference",final_data)





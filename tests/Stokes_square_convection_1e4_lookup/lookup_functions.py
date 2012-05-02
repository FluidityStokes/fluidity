#!/usr/bin/env python   

import numpy
import math

# Start by specifying key parameters of lookup table:
# First specify temperature info:
tmin = 300.     # Minimum T (k)
tmax = 4500.    # Maximum T (k)
tinc = 10.      # Temperature increment (k)
# And pressure info:
pmin = 0.       # Minimum pressure (GPa)
pmax = 140.      # Maximum pressure (GPa)
pinc = 0.1      # Pressure increment (GPa)

# Use this information to form arrays of lookup temperature and pressure fields:
Lookup_Ts = numpy.arange(tmin,tmax+tinc,tinc)
Lookup_Ps = numpy.arange(pmin,pmax+pinc,pinc)

# And for the reference state (used for reference temperauture (K)):
adiabat_footing_temperature = 1600.

def get_attribute_column(field_name):
  if(field_name == "FullLookupDensity"):
    att_column = 0
  if(field_name == "Enthalpy"):
    att_column = 1
  if(field_name == "IsobaricSpecificHeatCapacity"):
    att_column = 2
  if(field_name == "IsobaricThermalExpansivity"):
    att_column = 3
  return att_column

def get_reference_attribute_column(field_name):
  if(field_name == "Depth"):
    att_column = 0
  if(field_name == "Reference_Pressure"):
    att_column = 1
  if(field_name == "CompressibleReferenceTemperature"):
    att_column = 2
  if(field_name == "CompressibleReferenceDensity"):
    att_column = 3
  return att_column

def find_nearest_index(array,value):
  # Function to find index of nearest value
  idx=array.searchsorted(value)
  return idx-1

def central_difference(H_up, H_down, increment):
  # Central difference H_up and H_down:
  val = (H_up-H_down) / increment
  return val

def get_material_property(lookup_data, T, P, field_name):
  # Use field name to determine which column of lookup table is appropriate for this field:
  att_column = get_attribute_column(field_name)
  # This function determines a value
  iP = int(numpy.floor((P-Lookup_Ps[0]) / pinc)) 
  iT = int(numpy.floor((T-Lookup_Ts[0]) / tinc))
  # Now interpolate data from surrounding (T,P) space to node:
  # First normalize values - begin with pressure:
  norm_pres = 1. - ((P - Lookup_Ps[iP]) / pinc)
  assert(norm_pres >= 0. and norm_pres <= 1.)
  # Next normalize temperature:
  norm_temp = 1. - ((T - Lookup_Ts[iT]) / tinc)
  assert(norm_temp >= 0. and norm_temp <= 1.)
  # Calculate contributions from each vertex:
  v1n = norm_temp     * norm_pres      * lookup_data[iT   ,iP   , att_column]
  v2n = (1.-norm_temp)* norm_pres      * lookup_data[iT+1 ,iP   , att_column]
  v3n = (1.-norm_temp)* (1.-norm_pres) * lookup_data[iT+1 ,iP+1 , att_column]
  v4n = norm_temp     * (1.-norm_pres) * lookup_data[iT   ,iP+1 , att_column]
  # Sum all contributions to get final value at this node:
  nodal_value = v1n+v2n+v3n+v4n
  return nodal_value

def get_reference_property(reference_data,X,field_name):
  att_column = get_reference_attribute_column(field_name)
  depth = -X[1]
  if (depth < 1e-50) :
    depth = 1e-50
  nearest_index = find_nearest_index(reference_data[:,0],depth)
  # Now interpolate data from surrounding depth points:
  # First normalize distances:
  depth_range = reference_data[nearest_index+1,0] - reference_data[nearest_index,0]
  d1 = 1. - (depth - reference_data[nearest_index,0]) / depth_range
  d2 = 1. - (reference_data[nearest_index+1,0] - depth) / depth_range
  assert(d1 >= 0. and d1 <= 1.)
  assert(d2 >= 0. and d2 <= 1.)
  reference_value = (d1*reference_data[nearest_index,att_column]+d2*reference_data[nearest_index+1,att_column])
  if(field_name == "CompressibleReferenceTemperature"):
    reference_value = reference_value - adiabat_footing_temperature
  return reference_value

def calculate_effective_cp(lookup_data, T, P):
  # Use field name to determine which column of lookup table is the enthalpy:
  enth_column = get_attribute_column("Enthalpy")
  # Find nearest T & P values in lookup tables (lower bounds):
  iP = int(numpy.floor((P-Lookup_Ps[0]) / pinc)) 
  iT = int(numpy.floor((T-Lookup_Ts[0]) / tinc))
  # Now interpolate data from surrounding (T,P) space to node:
  # First normalize values - begin with pressure:
  norm_pres = 1. - ((P - Lookup_Ps[iP]) / pinc)
  assert(norm_pres >= 0. and norm_pres <= 1.)
  # Next normalize temperature:
  norm_temp = 1. - ((T - Lookup_Ts[iT]) / tinc)
  assert(norm_temp >= 0. and norm_temp <= 1.)
  # Calculate Cp eff (dH/dT)p on each surrounding point (11,21,22,12) using a simple central difference scheme. Note that
  # at the surface (T = tmin or P = pmin), a forward difference scheme is used to account for the lack of data:
  #-----------------#
  H_up   = (lookup_data[iT+1   ,iP   , enth_column] + lookup_data[iT          ,iP   , enth_column]) / 2.
  H_down = (lookup_data[iT     ,iP   , enth_column] + lookup_data[max(0,iT-1) ,iP   , enth_column]) / 2.
  if((iT-1) >= 0):
    H11  = central_difference(H_up,H_down,tinc)
  else : # Forward difference
    H11  = (lookup_data[iT+1   ,iP   , enth_column] - lookup_data[iT ,iP   , enth_column]) / tinc
  #-----------------#
  H_up   = (lookup_data[iT+2   ,iP   , enth_column] + lookup_data[iT+1         ,iP   , enth_column]) / 2.
  H_down = (lookup_data[iT+1   ,iP   , enth_column] + lookup_data[iT           ,iP   , enth_column]) / 2.
  H21    = central_difference(H_up,H_down,tinc)
  #-----------------#
  H_up   = (lookup_data[iT+2   ,iP+1 , enth_column] + lookup_data[iT+1         ,iP+1 , enth_column]) / 2.
  H_down = (lookup_data[iT+1   ,iP+1 , enth_column] + lookup_data[iT           ,iP+1 , enth_column]) / 2.
  H22    = central_difference(H_up,H_down,tinc)
  #-----------------#
  H_up   = (lookup_data[iT+1   ,iP+1 , enth_column] + lookup_data[iT           ,iP+1 , enth_column]) / 2.
  H_down = (lookup_data[iT     ,iP+1 , enth_column] + lookup_data[max(0,iT-1)  ,iP+1 , enth_column]) / 2.
  if((iT-1) >= 0):
    H12  = central_difference(H_up,H_down,tinc)
  else : # Forward difference
    H12  = (lookup_data[iT+1   ,iP+1 , enth_column] - lookup_data[iT, iP+1 , enth_column]) / tinc
  #-----------------#
  # Calculate contributions from each vertex to nodal point:
  v1n = norm_temp     * norm_pres      * H11
  v2n = (1.-norm_temp)* norm_pres      * H21
  v3n = (1.-norm_temp)* (1.-norm_pres) * H22
  v4n = norm_temp     * (1.-norm_pres) * H12
  # Sum all contributions to get final value at this node:
  effective_cp = v1n+v2n+v3n+v4n
  # Convert in correct units:
  effective_cp = effective_cp * 1.0e6
  return effective_cp

def calculate_effective_alpha(lookup_data, T, P):
  # Use field name to determine which column of lookup table is the enthalpy:
  enth_column = get_attribute_column("Enthalpy")
  rho_column  = get_attribute_column("FullLookupDensity")
  # Find nearest T & P values in lookup tables (lower bounds):
  iP = int(numpy.floor((P-Lookup_Ps[0]) / pinc)) 
  iT = int(numpy.floor((T-Lookup_Ts[0]) / tinc))
  # Now interpolate data from surrounding (T,P) space to node:
  # First normalize values - begin with pressure:
  norm_pres = 1. - ((P - Lookup_Ps[iP]) / pinc)
  assert(norm_pres >= 0. and norm_pres <= 1.)
  # Next normalize temperature:
  norm_temp = 1. - ((T - Lookup_Ts[iT]) / tinc)
  assert(norm_temp >= 0. and norm_temp <= 1.)
  # Calculate Alpha eff on each surrounding point using a simple central difference scheme:
  #-----------------#
  H_up   = (lookup_data[iT     ,iP+1 , enth_column] + lookup_data[iT          ,iP            , enth_column]) / 2. 
  H_down = (lookup_data[iT     ,iP   , enth_column] + lookup_data[iT          ,max(0,iP-1)   , enth_column]) / 2.
  if((iP-1) >= 0):
    H11  = central_difference(H_up, H_down, pinc)
  else : # Forward difference
    H11  = (( lookup_data[iT   ,iP+1 , enth_column] - lookup_data[iT ,iP , enth_column] ) / pinc )
  rho    = lookup_data[iT     ,iP   , rho_column]
  T      = (1. / Lookup_Ts[iT]) # 1 over T
  H11    = T * (1. - (rho*H11))
  #-----------------#
  H_up   = (lookup_data[iT+1   ,iP+1 , enth_column] + lookup_data[iT+1        ,iP            , enth_column]) / 2.
  H_down = (lookup_data[iT+1   ,iP   , enth_column] + lookup_data[iT+1        ,max(0,iP-1)   , enth_column]) / 2.
  if((iP-1) >= 0):
    H21  = central_difference(H_up, H_down, pinc)
  else : # Forward difference
    H21  = (( lookup_data[iT+1 ,iP+1 , enth_column] - lookup_data[iT+1, iP , enth_column] )/ pinc )
  rho    = lookup_data[iT+1     ,iP   , rho_column]
  T      = (1. / Lookup_Ts[iT+1]) # 1 over T
  H21    = T * (1. - (rho*H21))
  #-----------------#
  H_up   = (lookup_data[iT+1   ,iP+2 , enth_column] + lookup_data[iT+1        ,iP+1          , enth_column]) / 2.
  H_down = (lookup_data[iT+1   ,iP+1 , enth_column] + lookup_data[iT+1        ,iP            , enth_column]) / 2.
  H22    = central_difference(H_up, H_down, pinc) 
  rho    = lookup_data[iT+1     ,iP+1   , rho_column]
  T      = (1. / Lookup_Ts[iT+1]) # 1 over T
  H22    = T * (1. - (rho*H22))
  #-----------------#
  H_up   = (lookup_data[iT     ,iP+2 , enth_column] + lookup_data[iT          ,iP+1          , enth_column]) / 2.
  H_down = (lookup_data[iT     ,iP+1 , enth_column] + lookup_data[iT          ,iP            , enth_column]) / 2.
  H12    = central_difference(H_up, H_down, pinc)
  rho    = lookup_data[iT     ,iP+1   , rho_column]
  T      = (1. / Lookup_Ts[iT]) # 1 over T
  H12    = T * (1. - (rho*H12))
  #-----------------#
  # Calculate contributions from each vertex to nodal point:
  v1n = norm_temp     * norm_pres      * H11
  v2n = (1.-norm_temp)* norm_pres      * H21
  v3n = (1.-norm_temp)* (1.-norm_pres) * H22
  v4n = norm_temp     * (1.-norm_pres) * H12
  # Sum all contributions to get final value at this node:
  effective_alpha = v1n+v2n+v3n+v4n
  return effective_alpha


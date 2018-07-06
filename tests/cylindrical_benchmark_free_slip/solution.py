from __future__ import division
import numpy
from math import sqrt, atan2, cos, sin, pi

nu = 1.
g = 1.
# outer and inner radius R_+ and R_-
Rp, Rm = 2.22, 1.22
# radial height of anomaly: r'
rp = Rm + 0.5
n = 2

if n<=1:
    raise NotImplemented()

alpha_pm = numpy.array([Rp/rp, Rm/rp])
alpha_mp = numpy.array([Rm/rp, Rp/rp])
alpha, beta = alpha_pm
pm = numpy.array([1, -1])
mp = -pm

# velocity solution: coefficients for n, -n, n+2, and -n+2 power of r
# since the alpha's are length-2 arrays (each entry corresponding to one halve of the domain)
# so will the coefficients and thus the solution
A_pm = -0.125*(alpha_mp**(2*n - 2) - 1)*g*pm*rp**(-n + 2)/((alpha_mp**(2*n - 2) - alpha_pm**(2*n - 2))*(n - 1)*nu)
B_pm = -0.125*(alpha_mp**(2*n + 2) - 1)*alpha_pm**(2*n + 2)*g*pm*rp**(n + 2)/((alpha_mp**(2*n + 2) - alpha_pm**(2*n + 2))*(n + 1)*nu)
C_pm = 0.125*(alpha_mp**(2*n + 2) - 1)*g*pm*rp**(-n)/((alpha_mp**(2*n + 2) - alpha_pm**(2*n + 2))*(n + 1)*nu)
D_pm = 0.125*(alpha_mp**(2*n - 2) - 1)*alpha_pm**(2*n - 2)*g*pm*rp**n/((alpha_mp**(2*n - 2) - alpha_pm**(2*n - 2))*(n - 1)*nu)


# pressure solution: coefficients for n and -n
G_pm = -4*nu*C_pm*(n+1)
H_pm = -4*nu*D_pm*(n-1)


def u_r(r, phi):
    dpsi_dphi = n*cos(n*phi)*(A_pm*r**n+B_pm*r**(-n)+C_pm*r**(n+2)+D_pm*r**(-n+2))
    return -dpsi_dphi/r

# some sanity checks:
# no-normal flow
numpy.testing.assert_almost_equal(u_r(Rm, 0.)[1], 0.)
numpy.testing.assert_almost_equal(u_r(Rp, 0.)[0], 0.)
# coninuity of u_r
numpy.testing.assert_almost_equal(u_r(rp, 0.)[0], u_r(rp, 0.)[1])

def u_phi(r, phi):
    dpsi_dr = sin(n*phi)*(A_pm*n*r**(n-1) + B_pm*-n*r**(-n-1) + C_pm*(n+2)*r**(n+1) + D_pm*(-n+2)*r**(-n+1))
    return dpsi_dr

# coninuity of u_phi
numpy.testing.assert_almost_equal(u_phi(rp, pi/(2*n))[0], u_phi(rp, pi/(2*n))[1])



def p(r, phi):
    return (G_pm*r**n + H_pm*r**(-n))*cos(n*phi)

def velocity_cartesian(X, i):
  # i==0: upper mantle, i==1: lower mantle
  r = sqrt(X[0]**2+X[1]**2)
  phi = atan2(X[1], X[0])
  ur = u_r(r, phi)[i]
  ut = u_phi(r, phi)[i]
  return [ur*X[0]/r - ut*X[1]/r, ur*X[1]/r + ut*X[0]/r]

def pressure_cartesian(X, i):
  # i==0: upper mantle, i==1: lower mantle
  r = sqrt(X[0]**2+X[1]**2)
  phi = atan2(X[1], X[0])
  return p(r, phi)[i]

from math import sin, cos, tanh, pi, e, sqrt

def u(X):
    return 0.150*sin(0.750*X[1]) + cos(0.900*X[0]) + 0.150*cos(0.750*X[1]) + 0.500

def v(X):
    return 0.9*X[1]*sin(0.9*X[0])

def p(X):
    return sin(X[0]*X[1]/pi) + sin(X[0]) + cos(X[1]) - 1.00

def s1(X):
    return 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.70*X[0]) - 0.0200*cos(2.10*X[1]) + 0.100

def s2(X):
    return -0.100*sin(0.600*X[0]*X[1]/pi) + 0.0100*sin(1.40*X[0]) - 0.0100*cos(3.00*X[1]) + 0.100

def s1_d(X):
    return -0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) + 10.0

def s2_d(X):
    return 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) - 0.200*cos(1.30*X[1]) + 5.00

def rho(X):
    return -0.0660*sin(0.600*X[0]*X[1]/pi) + 0.0231*sin(1.30*X[0]*X[1]/pi) + 0.00660*sin(1.40*X[0]) - 0.00330*sin(1.70*X[0]) - 0.00660*cos(2.10*X[1]) - 0.00660*cos(3.00*X[1]) + 0.0990

def nu(X):
    return (0.153846153846154*sin(0.600*X[0]*X[1]/pi) - 0.107692307692308*sin(1.30*X[0]*X[1]/pi) - 0.0153846153846154*sin(1.40*X[0]) + 0.0153846153846154*sin(1.70*X[0]) + 0.0307692307692308*cos(2.10*X[1]) + 0.0153846153846154*cos(3.00*X[1]) + 0.692307692307692)**(-1.62500000000000)

def mu(X):
    return (-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(2.30*X[0]*X[1]/pi) - 0.200*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.400*cos(1.30*X[1]) + 20.0)/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0)

def sigma(X):
    return (((-(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(2.30*X[0]*X[1]/pi) - 0.200*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.400*cos(1.30*X[1]) + 20.0)/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0) + 1.00)**2.00*(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) + 10.0) + (-(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(2.30*X[0]*X[1]/pi) - 0.200*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.400*cos(1.30*X[1]) + 20.0)/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0) + 2.00)**2.00*(0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) - 0.200*cos(1.30*X[1]) + 5.00))/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0))**0.500

def tau_b(X):
    return (0.150*sin(0.750*X[1]) + cos(0.900*X[0]) + 0.150*cos(0.750*X[1]) + 0.500)**2.00 + 0.81*(X[1]*sin(0.9*X[0]))**2.00

def E_1(X):
    if X[1] < 1e-6:
        return (2.07806200865388e-6)*(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) + 10.0)*(0.100*sin(0.600*X[0]*X[1]/pi) - 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.40*X[0]) + 0.0100*sin(1.70*X[0]) + 0.0200*cos(2.10*X[1]) + 0.0100*cos(3.00*X[1]) + 0.800)**2.39*((-0.288*(((-(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(2.30*X[0]*X[1]/pi) - 0.200*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.400*cos(1.30*X[1]) + 20.0)/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0) + 1.00)**2.00*(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) + 10.0) + (-(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(2.30*X[0]*X[1]/pi) - 0.200*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.400*cos(1.30*X[1]) + 20.0)/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0) + 2.00)**2.00*(0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) - 0.200*cos(1.30*X[1]) + 5.00))/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0))**0.500 + 1.00)*((0.150*sin(0.750*X[1]) + cos(0.900*X[0]) + 0.150*cos(0.750*X[1]) + 0.500)**2.00 + 0.81*(X[1]*sin(0.9*X[0]))**2.00)**0.500*((0.153846153846154*sin(0.600*X[0]*X[1]/pi) - 0.107692307692308*sin(1.30*X[0]*X[1]/pi) - 0.0153846153846154*sin(1.40*X[0]) + 0.0153846153846154*sin(1.70*X[0]) + 0.0307692307692308*cos(2.10*X[1]) + 0.0153846153846154*cos(3.00*X[1]) + 0.692307692307692)**1.62500000000000)**0.600/(0.100*sin(0.600*X[0]*X[1]/pi) - 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.40*X[0]) + 0.0100*sin(1.70*X[0]) + 0.0200*cos(2.10*X[1]) + 0.0100*cos(3.00*X[1]) + 0.800)**2.39)**5.00/((0.0000209905253399382*((-0.288*(((-(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(2.30*X[0]*X[1]/pi) - 0.200*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.400*cos(1.30*X[1]) + 20.0)/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0) + 1.00)**2.00*(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) + 10.0) + (-(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(2.30*X[0]*X[1]/pi) - 0.200*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.400*cos(1.30*X[1]) + 20.0)/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0) + 2.00)**2.00*(0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) - 0.200*cos(1.30*X[1]) + 5.00))/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0))**0.500 + 1.00)*((0.150*sin(0.750*X[1]) + cos(0.900*X[0]) + 0.150*cos(0.750*X[1]) + 0.500)**2.00 + 0.81*(X[1]*sin(0.9*X[0]))**2.00)**0.500*((0.153846153846154*sin(0.600*X[0]*X[1]/pi) - 0.107692307692308*sin(1.30*X[0]*X[1]/pi) - 0.0153846153846154*sin(1.40*X[0]) + 0.0153846153846154*sin(1.70*X[0]) + 0.0307692307692308*cos(2.10*X[1]) + 0.0153846153846154*cos(3.00*X[1]) + 0.692307692307692)**1.62500000000000)**0.600/(0.100*sin(0.600*X[0]*X[1]/pi) - 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.40*X[0]) + 0.0100*sin(1.70*X[0]) + 0.0200*cos(2.10*X[1]) + 0.0100*cos(3.00*X[1]) + 0.800)**2.39)**5.00 + 1.00)*(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0))
    else:
        return 0.0

def E_2(X):
    if X[1] < 1e-6:
        return 0.0000166244960692310*(0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) - 0.200*cos(1.30*X[1]) + 5.00)*(0.100*sin(0.600*X[0]*X[1]/pi) - 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.40*X[0]) + 0.0100*sin(1.70*X[0]) + 0.0200*cos(2.10*X[1]) + 0.0100*cos(3.00*X[1]) + 0.800)**2.39*((-0.288*(((-(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(2.30*X[0]*X[1]/pi) - 0.200*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.400*cos(1.30*X[1]) + 20.0)/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0) + 1.00)**2.00*(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) + 10.0) + (-(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(2.30*X[0]*X[1]/pi) - 0.200*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.400*cos(1.30*X[1]) + 20.0)/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0) + 2.00)**2.00*(0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) - 0.200*cos(1.30*X[1]) + 5.00))/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0))**0.500 + 1.00)*((0.150*sin(0.750*X[1]) + cos(0.900*X[0]) + 0.150*cos(0.750*X[1]) + 0.500)**2.00 + 0.81*(X[1]*sin(0.9*X[0]))**2.00)**0.500*((0.153846153846154*sin(0.600*X[0]*X[1]/pi) - 0.107692307692308*sin(1.30*X[0]*X[1]/pi) - 0.0153846153846154*sin(1.40*X[0]) + 0.0153846153846154*sin(1.70*X[0]) + 0.0307692307692308*cos(2.10*X[1]) + 0.0153846153846154*cos(3.00*X[1]) + 0.692307692307692)**1.62500000000000)**0.600/(0.100*sin(0.600*X[0]*X[1]/pi) - 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.40*X[0]) + 0.0100*sin(1.70*X[0]) + 0.0200*cos(2.10*X[1]) + 0.0100*cos(3.00*X[1]) + 0.800)**2.39)**5.00/((0.0000839621013597527*((-0.288*(((-(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(2.30*X[0]*X[1]/pi) - 0.200*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.400*cos(1.30*X[1]) + 20.0)/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0) + 1.00)**2.00*(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) + 10.0) + (-(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.200*sin(2.30*X[0]*X[1]/pi) - 0.200*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.400*cos(1.30*X[1]) + 20.0)/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0) + 2.00)**2.00*(0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) - 0.200*cos(1.30*X[1]) + 5.00))/(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0))**0.500 + 1.00)*((0.150*sin(0.750*X[1]) + cos(0.900*X[0]) + 0.150*cos(0.750*X[1]) + 0.500)**2.00 + 0.81*(X[1]*sin(0.9*X[0]))**2.00)**0.500*((0.153846153846154*sin(0.600*X[0]*X[1]/pi) - 0.107692307692308*sin(1.30*X[0]*X[1]/pi) - 0.0153846153846154*sin(1.40*X[0]) + 0.0153846153846154*sin(1.70*X[0]) + 0.0307692307692308*cos(2.10*X[1]) + 0.0153846153846154*cos(3.00*X[1]) + 0.692307692307692)**1.62500000000000)**0.600/(0.100*sin(0.600*X[0]*X[1]/pi) - 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.40*X[0]) + 0.0100*sin(1.70*X[0]) + 0.0200*cos(2.10*X[1]) + 0.0100*cos(3.00*X[1]) + 0.800)**2.39)**5.00 + 1.00)*(-0.200*sin(0.300*X[0]*X[1]/pi) + 0.100*sin(2.30*X[0]*X[1]/pi) - 0.100*sin(1.40*X[0]) + 0.200*sin(1.70*X[0]) + 0.400*cos(1.10*X[1]) - 0.200*cos(1.30*X[1]) + 15.0))
    else:
        return 0.0

def D_1(X):
    if X[1] < 1e-6:
        return -(0.330*(0.100*sin(0.600*X[0]*X[1]/pi) - 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.40*X[0]) + 0.0100*sin(1.70*X[0]) + 0.0200*cos(2.10*X[1]) + 0.0100*cos(3.00*X[1]) + 0.800)**2.39 + 0.9*X[1]*sin(0.9*X[0]))*(0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.70*X[0]) - 0.0200*cos(2.10*X[1]) + 0.100)
    else:
        return 0.0

def D_2(X):
    if X[1] < 1e-6:
        return -(0.660*(0.100*sin(0.600*X[0]*X[1]/pi) - 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.40*X[0]) + 0.0100*sin(1.70*X[0]) + 0.0200*cos(2.10*X[1]) + 0.0100*cos(3.00*X[1]) + 0.800)**2.39 + 0.9*X[1]*sin(0.9*X[0]))*(-0.100*sin(0.600*X[0]*X[1]/pi) + 0.0100*sin(1.40*X[0]) - 0.0100*cos(3.00*X[1]) + 0.100)
    else:
        return 0.0

def forcing_u(X):
    return 0.9*(-0.112500000000000*sin(0.750*X[1]) + 0.112500000000000*cos(0.750*X[1]))*X[1]*sin(0.9*X[0]) - 0.900*(0.150*sin(0.750*X[1]) + cos(0.900*X[0]) + 0.150*cos(0.750*X[1]) + 0.500)*sin(0.900*X[0]) + X[1]*cos(X[0]*X[1]/pi)/pi - (-0.0843750000000000*sin(0.750*X[1]) + 0.81*cos(0.9*X[0]) - 0.0843750000000000*cos(0.750*X[1]))/(0.153846153846154*sin(0.600*X[0]*X[1]/pi) - 0.107692307692308*sin(1.30*X[0]*X[1]/pi) - 0.0153846153846154*sin(1.40*X[0]) + 0.0153846153846154*sin(1.70*X[0]) + 0.0307692307692308*cos(2.10*X[1]) + 0.0153846153846154*cos(3.00*X[1]) + 0.692307692307692)**1.62500000000000 + 1.62500000000000*(0.81*X[1]*cos(0.9*X[0]) - 0.112500000000000*sin(0.750*X[1]) + 0.112500000000000*cos(0.750*X[1]))*(0.0923076923076923*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 0.140*X[0]*cos(1.30*X[0]*X[1]/pi)/pi - 0.0646153846153846*sin(2.10*X[1]) - 0.0461538461538461*sin(3.00*X[1]))/(0.153846153846154*sin(0.600*X[0]*X[1]/pi) - 0.107692307692308*sin(1.30*X[0]*X[1]/pi) - 0.0153846153846154*sin(1.40*X[0]) + 0.0153846153846154*sin(1.70*X[0]) + 0.0307692307692308*cos(2.10*X[1]) + 0.0153846153846154*cos(3.00*X[1]) + 0.692307692307692)**2.62500000000000 - 2.92500000000000*(0.0923076923076923*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 0.140*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 0.0215384615384615*cos(1.40*X[0]) + 0.0261538461538462*cos(1.70*X[0]))*sin(0.900*X[0])/(0.153846153846154*sin(0.600*X[0]*X[1]/pi) - 0.107692307692308*sin(1.30*X[0]*X[1]/pi) - 0.0153846153846154*sin(1.40*X[0]) + 0.0153846153846154*sin(1.70*X[0]) + 0.0307692307692308*cos(2.10*X[1]) + 0.0153846153846154*cos(3.00*X[1]) + 0.692307692307692)**2.62500000000000 + 1.62*cos(0.900*X[0])/(0.153846153846154*sin(0.600*X[0]*X[1]/pi) - 0.107692307692308*sin(1.30*X[0]*X[1]/pi) - 0.0153846153846154*sin(1.40*X[0]) + 0.0153846153846154*sin(1.70*X[0]) + 0.0307692307692308*cos(2.10*X[1]) + 0.0153846153846154*cos(3.00*X[1]) + 0.692307692307692)**1.62500000000000 + 0.0660*sin(0.600*X[0]*X[1]/pi) - 0.0231*sin(1.30*X[0]*X[1]/pi) - 0.00660*sin(1.40*X[0]) + 0.00330*sin(1.70*X[0]) + cos(X[0]) + 0.00660*cos(2.10*X[1]) + 0.00660*cos(3.00*X[1]) - 0.0990

def forcing_v(X):
    return 0.81*(0.150*sin(0.750*X[1]) + cos(0.900*X[0]) + 0.150*cos(0.750*X[1]) + 0.500)*X[1]*cos(0.9*X[0]) + 0.81*X[1]*sin(0.9*X[0])**2 + X[0]*cos(X[0]*X[1]/pi)/pi + 0.729*X[1]*sin(0.9*X[0])/(0.153846153846154*sin(0.600*X[0]*X[1]/pi) - 0.107692307692308*sin(1.30*X[0]*X[1]/pi) - 0.0153846153846154*sin(1.40*X[0]) + 0.0153846153846154*sin(1.70*X[0]) + 0.0307692307692308*cos(2.10*X[1]) + 0.0153846153846154*cos(3.00*X[1]) + 0.692307692307692)**1.62500000000000 + 1.62500000000000*(0.81*X[1]*cos(0.9*X[0]) - 0.112500000000000*sin(0.750*X[1]) + 0.112500000000000*cos(0.750*X[1]))*(0.0923076923076923*X[1]*cos(0.600*X[0]*X[1]/pi)/pi - 0.140*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 0.0215384615384615*cos(1.40*X[0]) + 0.0261538461538462*cos(1.70*X[0]))/(0.153846153846154*sin(0.600*X[0]*X[1]/pi) - 0.107692307692308*sin(1.30*X[0]*X[1]/pi) - 0.0153846153846154*sin(1.40*X[0]) + 0.0153846153846154*sin(1.70*X[0]) + 0.0307692307692308*cos(2.10*X[1]) + 0.0153846153846154*cos(3.00*X[1]) + 0.692307692307692)**2.62500000000000 + 2.925*(0.0923076923076923*X[0]*cos(0.600*X[0]*X[1]/pi)/pi - 0.140*X[0]*cos(1.30*X[0]*X[1]/pi)/pi - 0.0646153846153846*sin(2.10*X[1]) - 0.0461538461538461*sin(3.00*X[1]))*sin(0.9*X[0])/(0.153846153846154*sin(0.600*X[0]*X[1]/pi) - 0.107692307692308*sin(1.30*X[0]*X[1]/pi) - 0.0153846153846154*sin(1.40*X[0]) + 0.0153846153846154*sin(1.70*X[0]) + 0.0307692307692308*cos(2.10*X[1]) + 0.0153846153846154*cos(3.00*X[1]) + 0.692307692307692)**2.62500000000000 + 0.0660*sin(0.600*X[0]*X[1]/pi) - 0.0231*sin(1.30*X[0]*X[1]/pi) - 0.00660*sin(1.40*X[0]) + 0.00330*sin(1.70*X[0]) - sin(X[1]) + 0.00660*cos(2.10*X[1]) + 0.00660*cos(3.00*X[1]) - 0.0990

def forcing_s1(X):
    return (0.0910*X[1]*cos(1.30*X[0]*X[1]/pi)/pi - 0.0170*cos(1.70*X[0]))*(0.330*(0.100*sin(0.600*X[0]*X[1]/pi) - 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.40*X[0]) + 0.0100*sin(1.70*X[0]) + 0.0200*cos(2.10*X[1]) + 0.0100*cos(3.00*X[1]) + 0.800)**2.39 + 0.150*sin(0.750*X[1]) + cos(0.900*X[0]) + 0.150*cos(0.750*X[1]) + 0.500) + (0.0910*X[0]*cos(1.30*X[0]*X[1]/pi)/pi + 0.0420*sin(2.10*X[1]))*(0.330*(0.100*sin(0.600*X[0]*X[1]/pi) - 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.40*X[0]) + 0.0100*sin(1.70*X[0]) + 0.0200*cos(2.10*X[1]) + 0.0100*cos(3.00*X[1]) + 0.800)**2.39 + 0.9*X[1]*sin(0.9*X[0])) + 0.118300000000000*X[0]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 + 0.118300000000000*X[1]**2*sin(1.30*X[0]*X[1]/pi)/pi**2 - 0.0289*sin(1.70*X[0]) - 0.0882*cos(2.10*X[1])

def forcing_s2(X):
    return (-0.0600*X[1]*cos(0.600*X[0]*X[1]/pi)/pi + 0.0140*cos(1.40*X[0]))*(0.660*(0.100*sin(0.600*X[0]*X[1]/pi) - 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.40*X[0]) + 0.0100*sin(1.70*X[0]) + 0.0200*cos(2.10*X[1]) + 0.0100*cos(3.00*X[1]) + 0.800)**2.39 + 0.150*sin(0.750*X[1]) + cos(0.900*X[0]) + 0.150*cos(0.750*X[1]) + 0.500) + (-0.0600*X[0]*cos(0.600*X[0]*X[1]/pi)/pi + 0.0300*sin(3.00*X[1]))*(0.660*(0.100*sin(0.600*X[0]*X[1]/pi) - 0.0700*sin(1.30*X[0]*X[1]/pi) - 0.0100*sin(1.40*X[0]) + 0.0100*sin(1.70*X[0]) + 0.0200*cos(2.10*X[1]) + 0.0100*cos(3.00*X[1]) + 0.800)**2.39 + 0.9*X[1]*sin(0.9*X[0])) - 0.0360*X[0]**2*sin(0.600*X[0]*X[1]/pi)/pi**2 - 0.0360*X[1]**2*sin(0.600*X[0]*X[1]/pi)/pi**2 + 0.0196*sin(1.40*X[0]) - 0.0900*cos(3.00*X[1])

def velocity(X):
    return [u(X), v(X)]

def forcing_velocity(X):
    return [forcing_u(X), forcing_v(X)]

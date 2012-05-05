# -*- coding:utf-8 -*-

import numpy as np;
import re; # regexp
import matplotlib.pyplot as mp;

################################################################
# Airfoil : load profile of a wing
#
# Reads a file whose lines contain coordinates of points,
# separated by an empty line.
# Every line not containing a couple of floats is discarded. 
# Returns a couple constitued of the list of points of the
# extrados and the intrados. 
def load_foil(file):
    f = open(file, 'r')
    matchline = lambda line: re.match(r"\s*([\d\.-]+)\s*([\d\.-]+)", line)
    extra  = [];    intra = []
    rextra = False; rintra = False
    for line in f:
        m = matchline(line)
        if (m != None) and not(rextra):
            rextra = True
        if (m != None) and rextra and not(rintra):
            extra.append(m.groups())
        if (m != None) and rextra and rintra:
            intra.append(m.groups())
        if (m == None) and rextra:
            rintra = True
    ex = np.array(map(lambda t: float(t[0]),extra))
    ey = np.array(map(lambda t: float(t[1]),extra))
    ix = np.array(map(lambda t: float(t[0]),intra))
    iy = np.array(map(lambda t: float(t[1]),intra))
    return(ex,ey,ix,iy)


def plot_airfoil(ex, ey, ix, iy, file_name):
    """ Plots an airfoil (or any pair of curves) """
    mp.clf()
    mp.plot(ex, ey, 'ro')
    mp.plot(ix, iy, 'go')
    mp.savefig(file_name)
    mp.clf()


# -------------------------------------------------- #
# Cubic Splines :                                    #
# -------------------------------------------------- #

def cubic_spline_interpolation_precalc(X,Y):
    """ Pre-calculates the second derivative used for cubic_spline_interpolation for the X, Y cloud of points and returns it """

    N = X.size

    # --- The second derivative in each point ---
    Ypp = np.zeros(N)

    # --- Building the system to solve to find Ypp ---
    Mat_coeff = np.zeros([N,N])
    Images = np.zeros(N)
    Mat_coeff[0,0] = 1 # Limit conditions
    Mat_coeff[N-1,N-1] = 1
    Images[0] = 0
    Images[1] = 0
    
    for i in np.arange(1,N-1):
        Mat_coeff[i, i-1] = (X[i] - X[i-1])/6
        Mat_coeff[i, i] = (X[i+1] - X[i-1])/3
        Mat_coeff[i, i+1] = (X[i+1] - X[i])/6
        Images[i] = (Y[i+1] - Y[i])/(X[i+1] - X[i]) - (Y[i] - Y[i-1])/(X[i] - X[i-1])
    
    # --- Solving the system ---
    Ypp = np.linalg.solve(Mat_coeff, Images)

    return Ypp
    
def cubic_spline_interpolation_calc(X,Y,Ypp,x):
    """ Returns the cubic spline interpolation of the cloud of point X, Y, knowing its second derivative Ypp """

    N = X.size

    for i in np.arange(0,N-1):
        if (X[i] <= x) and (x < X[i+1]):
            # --- The 4 factors  ---
            A = (X[i+1] - x) / (X[i+1] - X[i])
            B = 1 - A
            opt_var = (1/6) * (X[i+1] - X[i])**2 # Optimisation légère du calcul
            C = (A**3 - A) * opt_var
            D = (B**3 - B) * opt_var
            return A*Y[i] + B*Y[i+1] + C*Ypp[i] + D*Ypp[i+1]

    if x >= X[N-1]:
        return X[N-1]
    else:
        return X[0]
    


# A terme sera fait avec un lambda
def cubic_spline_interpolation(X,Y):
    """ Returns a f(x) function that interpolates the cloud of points X, Y using cubic Splines """

    N = X.size
    
    # --- The second derivative in each point ---
    Ypp = cubic_spline_interpolation_precalc(X,Y)

    # --- Building the lambda function to return ---
    return lambda x: cubic_spline_interpolation_calc(X,Y,Ypp,x)


    
# ------------------- Tests ------------------------ #

(ex,ey,ix,iy) = load_foil("fx63145.dat")
plot_airfoil(ex, ey, ix, iy, "airfoil.png")

ext_airfoil_interpolation = cubic_spline_interpolation(ex,ey)
int_airfoil_interpolation = cubic_spline_interpolation(ex,ey)



# -------------------------------------------------- #

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
    mp.plot(ex, ey, 'r')
    mp.plot(ix, iy, 'g')
    mp.savefig(file_name)
    mp.clf()

# ------------------- Tests ------------------------ #

(ex,ey,ix,iy) = load_foil("fx63145.dat")
#plot_airfoil(ex, ey, ix, iy, "airfoil.png")

# -------------------------------------------------- #


# -------------------------------------------------- #
# Cubic Splines :                                    #
# -------------------------------------------------- #

# A terme sera fait avec un lambda
def cubic_spline_interpolation(X,Y):
    """ Returns a f(x) function that interpolates the cloud of points X, Y using cubic Splines """

    # The 4 factors for each interval
    A = np.zeros(X.size)
    B = np.zeros(X.size)
    C = np.zeros(X.size)
    D = np.zeros(X.size)
    
    # The second derivative in each point
    Ypp = np.zeros(X.size)

    # Building the system to solve to find Ypp


    # Solving the system


    # Calculating the 4 factors for each interval


    # Building the lambda function to return
        

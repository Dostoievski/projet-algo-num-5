# -*- coding:utf-8 -*-

import numpy as np;
import re; # regexp
import matplotlib.pyplot as mp;
from splines import *


# ------------------- Tests ------------------------ #

def plot_airfoil(ex, ey, ix, iy, file_name):
    """ Plots an airfoil (or any pair of curves) """
    mp.clf()
    mp.plot(ex, ey, 'ro')
    mp.plot(ix, iy, 'go')
    mp.savefig(file_name)
    mp.clf()


def plot_points_and_interpolation(ex, ey, ix, iy, file_name):
    """ Plots an airfoil (or any pair of curves) and its cubic splines interpolation """
    
    ext_interp_fun = cubic_spline_interpolation(ex,ey)
    int_interp_fun = cubic_spline_interpolation(ix,iy)
    
    nb_of_points = 1000.
    min_interval = min(ex)
    max_interval = max(ex)
    pas = (max_interval - min_interval) / nb_of_points

    x = np.arange(min_interval, max_interval, pas)

    ext_airfoil_interpolation = np.zeros(nb_of_points)
    int_airfoil_interpolation = np.zeros(nb_of_points)

    for i in np.arange(0, nb_of_points):
        ext_airfoil_interpolation[i] = ext_interp_fun(x[i])
        int_airfoil_interpolation[i] = int_interp_fun(x[i])

    mp.clf()
    mp.plot(x, ext_airfoil_interpolation, 'r')
    mp.plot(x, int_airfoil_interpolation, 'g')
    mp.plot(ex, ey, 'r+')
    mp.plot(ix, iy, 'g+')
    mp.savefig(file_name)
    mp.clf()


def plot_interpolation_and_derivative(ex, ey, ix, iy, file_name):
    """ Plots an airfoil (or any pair of curves) and its cubic splines interpolation """
    
    ext_interp_fun = cubic_spline_interpolation(ex,ey)
    int_interp_fun = cubic_spline_interpolation(ix,iy)
    ext_interp_fun_derivative = cubic_spline_interpolation_derivative(ex,ey)
    int_interp_fun_derivative = cubic_spline_interpolation_derivative(ix,iy)

    
    nb_of_points = 1000.
    min_interval = min(ex)
    max_interval = max(ex)
    pas = (max_interval - min_interval) / nb_of_points

    x = np.arange(min_interval, max_interval, pas)

    ext_airfoil_interpolation = np.zeros(nb_of_points)
    int_airfoil_interpolation = np.zeros(nb_of_points)
    ext_airfoil_interpolation_derivative = np.zeros(nb_of_points)
    int_airfoil_interpolation_derivative = np.zeros(nb_of_points)


    for i in np.arange(0, nb_of_points):
        ext_airfoil_interpolation[i] = ext_interp_fun(x[i])
        int_airfoil_interpolation[i] = int_interp_fun(x[i])
        ext_airfoil_interpolation_derivative[i] = ext_interp_fun_derivative(x[i])
        int_airfoil_interpolation_derivative[i] = int_interp_fun_derivative(x[i])

    mp.clf()
    mp.plot(x, ext_airfoil_interpolation, 'r')
    mp.plot(x, int_airfoil_interpolation, 'g')
    mp.plot(x, ext_airfoil_interpolation_derivative, 'orange')
    mp.plot(x, int_airfoil_interpolation_derivative, 'b')
    mp.savefig(file_name)
    mp.clf()

    

(ex,ey,ix,iy) = load_foil("fx63145.dat")
plot_airfoil(ex, ey, ix, iy, "airfoil.png")

output_image_name = "interpolation_of_the_airfoil.png"

print "Saving the results of the interpolation in", output_image_name, "..."
plot_points_and_interpolation(ex, ey, ix, iy, output_image_name)
print "Result written.\n"

output_image_name = "interpolation_derivative.png"

print "Saving the results of the derivative in", output_image_name, "..."
plot_interpolation_and_derivative(ex, ey, ix, iy, output_image_name)
print "Result written.\n"



# -------------------------------------------------- #

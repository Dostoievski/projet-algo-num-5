import numpy as np
import scipy.integrate as sp
import splines as spl
import pylab as mp
import integration as it

(ex,ey,ix,iy) = spl.load_foil("fx63145.dat")

#-------------------- Function ------------------#

#------------------------------------------------#
# SPLINE FUNCTION :                              #
#------------------------------------------------#

ext_interp_fun = spl.cubic_spline_interpolation(ex, ey)
int_interp_fun = spl.cubic_spline_interpolation(ix, iy)

ext_interp_fun_d = spl.cubic_spline_interpolation_derivative(ex, ey)
int_interp_fun_d = spl.cubic_spline_interpolation_derivative(ix, iy)

#------------------------------------------------#
# LAMBDA FUNCTIONS :                             #
#------------------------------------------------#
  
def f_lambda(f, l, hmax):
    """ Computes the curves followed by the air particules. f : the function defining one half of the airfoil, l : lambda (greater lambda is, closer to the airfoil is the curve)"""
    return lambda(x):3*l*hmax + (1-l)*f(x)

def f_lambda_d(f_d, l):
    """ Derivative of f_lambda. f_d : derivative of f, l : lambda """
    return lambda(x):(1-l)*f_d(x)

def length_of_curve(f_d, integration_method, min_x, max_x, n):
    """ Returns the length of a curve. f_d : derivative of the curve, integration method : integration method to use, min_x, max_x : min and max bound of the interval to find the length on, n : number of points to take for the integration. """
    g = lambda(x): np.sqrt(1+f_d(x)**2)
    t = integration_method(g, min_x, max_x, n)
    return t

def pressure_on_curve(Ps, length):
    """ Returns the pressure on a curve. Ps : static pressure, length : longueur de la courbe """
    r = 1
    v = length/1
    return Ps + 0.5 * r * v**2

#------------------------------------------------#
# DRAW :                                         #
#------------------------------------------------#
     
def curve(f,n,t):
    curve_values = range(n)
    for i in range(n):
        curve_values[i] = f(t[i])
    return curve_values


def airflow_model(f_int, f_ext):
    """ Generates the paths followed by the air particules around the airfoil defined by f_int and f_ext """
    t = np.arange(0.0, 1.0, 0.01)
    lambdaTab = np.arange(0.0,1.0,0.1)  
    nbLambda = len(lambdaTab)
    
    for i in range(nbLambda):
        f_lambda_int = f_lambda(f_int, lambdaTab[i], min(iy))
        f_lambda_values_int = it.curve(f_lambda_int,len(t),t)
        mp.plot(t,f_lambda_values_int,linewidth=1.0, color='#0000ff')
        
    for j in range(nbLambda):
        f_lambda_ext = f_lambda(f_ext, lambdaTab[j], max(ey))
        f_lambda_values_ext = it.curve(f_lambda_ext,len(t),t)
        mp.plot(t,f_lambda_values_ext,linewidth=1.0, color='#0000ff')

    mp.axis([min(ex), max(max(ex),max(ix)), 3*min(iy), 3*max(ey)])
    mp.xlabel('x')
    mp.ylabel('y')
    mp.title('Laminar airflow of the fx63145')
    mp.savefig("airflows.png")
    mp.clf()

def pressure_to_colour(P, Pbase, Pmax):
    return ((P-Pbase)*10/Pmax)**3
    

def curves_pressure_representation(file_src, file_dest, integration_method, precision, n):
    """ Plots colored curves according to pressure """
    
    (ex, ey, ix, iy) = spl.load_foil(file_src)
    f_int = spl.cubic_spline_interpolation(ix, iy)
    f_ext = spl.cubic_spline_interpolation(ex, ey)
    f_int_d = spl.cubic_spline_interpolation_derivative(ix, iy)
    f_ext_d = spl.cubic_spline_interpolation_derivative(ex, ey)
    
    t = np.arange(0.0, 1.0, 0.01)
    lambdaTab = np.arange(0.0,1.0,0.01)
    nbLambda = len(lambdaTab)
    Pressure_Tab = np.zeros(nbLambda)
    Pmax = 1
    Ps = 0

    f_lambda_int_d = f_lambda_d(f_int_d, 0)
    length_int = length_of_curve(f_lambda_int_d, integration_method, min(ix), max(ix), n)
    Pbase = 0.9*pressure_on_curve(Ps, length_int)
    
    mp.clf()
    
    for l in lambdaTab:
        f_lambda_int = f_lambda(f_int, l, min(iy))
        f_lambda_int_d = f_lambda_d(f_int_d, l)
        length_int = length_of_curve(f_lambda_int_d, integration_method, min(ix), max(ix), n)
        pressure = pressure_on_curve(Ps, length_int)
        color = pressure_to_colour(pressure, Pbase, Pmax)
        f_lambda_values_int = it.curve(f_lambda_int,len(t),t)
        mp.plot(t,f_lambda_values_int,linewidth=5.0, color=(color, color**3,0))
        
    for l in reversed(lambdaTab):
        f_lambda_ext = f_lambda(f_ext, l, max(ey))
        f_lambda_ext_d = f_lambda_d(f_ext_d, l)
        length_ext = length_of_curve(f_lambda_ext_d, integration_method, min(ex), max(ex), n)
        pressure = pressure_on_curve(Ps, length_ext)
        color = pressure_to_colour(pressure, Pbase, Pmax)
        f_lambda_values_ext = it.curve(f_lambda_ext,len(t),t)
        mp.plot(t,f_lambda_values_ext,linewidth=5.0, color=(color,color**3,0))

    mp.axis([min(ex), max(max(ex),max(ix)), 3*min(iy), 3*max(ey)])
    mp.xlabel('x')
    mp.ylabel('y')
    mp.title('Map pressure of the fx63145 (The lighter corresponds to high pressure)')
    mp.savefig(file_dest)
    mp.clf()


#------------------------------------------------#
# TESTS :                                        #
#------------------------------------------------#

print " \n ## With simpsons integrals ### \n"
print " \n ########### EXTRADOS ######### \n "

f_lambda_ext_d_a = f_lambda_d(ext_interp_fun_d, 0.1)
f_lambda_ext_d_b = f_lambda_d(ext_interp_fun_d, 1.0)

l = 0.1
print "Calculation of the length of the exterior airfoil with lambda =", l
f_lambda_ext_d = f_lambda_d(ext_interp_fun_d, l)
length_ext_a = length_of_curve(f_lambda_ext_d, it.simpson, min(ex), max(ex), 100)
print length_ext_a

l = 1
print "Calculation of the length of the exterior airfoil with lambda =", l
f_lambda_ext_d = f_lambda_d(ext_interp_fun_d, l)
length_ext_b = length_of_curve(f_lambda_ext_d, it.simpson, min(ex), max(ex), 100)
print length_ext_b

Ps = 0
print "\nCalculation of pressure_int for f_lambda_a with lambda = 0.1"
print pressure_on_curve(Ps, length_ext_a)

print "Calculation of pressure_int for f_lambda_b with lambda = 1.0"
print pressure_on_curve(Ps, length_ext_b)

    
#M = matrix_map_pressure(ext_interp_fun_d, 1, 1, 0.1)
#print matrix_to_png(M, "Map_pressure.png")

airflow_model(int_interp_fun, ext_interp_fun)

curves_pressure_representation("fx63145.dat", "pressure_curves_simpson.png", it.simpson, 0.1, 100)
#curves_pressure_representation("fx63145.dat", "pressure_curves_rect_right.png", it.rectangle_method_right, 0.1, 100)

Mat = plot_pressures_around_airfoil_mat("fx63145.dat", 10.)
matrix_to_png(Mat, "test_pressure.png")



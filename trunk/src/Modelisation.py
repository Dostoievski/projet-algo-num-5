import numpy as np
import scipy.integrate as sp
import splines as spl
import pylab as mp
import integration as it

#--------Five ways to integrate a function--------#

#### 1) quadrature(func, a, b, args=(), tol=1.5e-8, maxiter=50)
#Integrate func from a to b using Gaussian quadrature with absolute tolerance tol. args is extra arguments to pass to func, and maxiter is the maximum number of iterations. Returns (val, err), where val is the Gaussian quadrature approximation, and err is the difference between the last two estimates of the integral.

#### 2) romberg(function, a, b, tol=1.48e-8, show=0, divmax=10)
#Romberg integration of a callable function or method.

#### 3) dblquad(func, a, b, gfun, hfun, args=(), epsabs=1.5e-8, epsrel=1.5e-8)
#Compute the double (definite) integral of func(y,x,*args) from x=a..b and y=gfun(x)..hfun(x). Returns (val, err).

#### 4) tplquad(func, a, b, gfun, hfun, qfun, rfun, args=(), epsabs=1.5e-8, epsrel=1.5e-8)
#Compute the triple (definite) integral of func(z,y,x,*args) from x=a..b, y=gfun(x)..hfun(x), and z=qfun(x,y)..rfun(x,y). Returns (val, err).

#### 5) quad(func, a, b, args=(), full_output=0, epsabs=1.49e-8, epsrel=1.49e-8, limit=50, points=None, weight=None, wvar=None, wopts=None, maxp1=50, limlst=50)
#Computes an integral using a technique from the Fortan library QUADPACK. The function is integrated from a to b. Run scipy.integrate.quad_explain() for more information on the more esoteric inputs and outputs.

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
        mp.plot(t,f_lambda_values_int,linewidth=1.0)
        
    for j in range(nbLambda):
        f_lambda_ext = f_lambda(f_ext, lambdaTab[j], max(ey))
        f_lambda_values_ext = it.curve(f_lambda_ext,len(t),t)
        mp.plot(t,f_lambda_values_ext,linewidth=1.0)

    mp.axis([min(ex), max(max(ex),max(ix)), 3*min(iy), 3*max(ey)])
    mp.title('Laminar airflow of the fx63145')
    mp.savefig("airflows.png")
    mp.clf()

def pressure_to_colour(P, Pbase, Pmax):
    return (P-Pbase)*10/Pmax
    

def curves_pressure_representation(f_int, f_int_d, f_ext,  f_ext_d, integration_method, Ps, precision, n):
    t = np.arange(0.0, 1.0, 0.01)
    lambdaTab = np.arange(0.0,1.0,0.1)
    nbLambda = len(lambdaTab)
    Pressure_Tab = np.zeros(nbLambda)
    Pmax = 1

    f_lambda_int_d = f_lambda_d(f_int_d, 0)
    length_int = length_of_curve(f_lambda_int_d, integration_method, min(ix), max(ix), n)
    Pbase = 0.9*pressure_on_curve(Ps, length_int)
    
    mp.clf()
    
    for i in range(nbLambda):
        f_lambda_int = f_lambda(f_int, lambdaTab[i], min(iy))
        f_lambda_int_d = f_lambda_d(f_int_d, lambdaTab[i])
        length_int = length_of_curve(f_lambda_int_d, integration_method, min(ix), max(ix), n)
        Pressure_Tab[i] = pressure_on_curve(Ps, length_int)
        color = pressure_to_colour(Pressure_Tab[i], Pbase, Pmax)
        f_lambda_values_int = it.curve(f_lambda_int,len(t),t)
        mp.plot(t,f_lambda_values_int,linewidth=5.0, color=(color, 0, 0))
        
    for j in range(nbLambda):
        f_lambda_ext = f_lambda(f_ext, lambdaTab[j], max(ey))
        f_lambda_ext_d = f_lambda_d(f_ext_d, lambdaTab[j])
        length_ext = length_of_curve(f_lambda_ext_d, integration_method, min(ex), max(ex), n)
        Pressure_Tab[j] = pressure_on_curve(Ps, length_ext)
        color = pressure_to_colour(Pressure_Tab[j], Pbase, Pmax)
        f_lambda_values_ext = it.curve(f_lambda_ext,len(t),t)
        mp.plot(t,f_lambda_values_ext,linewidth=5.0, color=(color,0,0))

    mp.axis([min(ex), max(max(ex),max(ix)), 3*min(iy), 3*max(ey)])
    mp.title('map pressure of the fx63145')
    mp.savefig("map_pressure_test.png")
    mp.clf()



#def matrix_map_pressure(f_d, Ps, precision):
#    lambdaTab = np.arange(0.0,1.0,precision)
#    nbLambda = len(lambdaTab)
#    M = np.zeros([nbLambda,nbLambda])
#    for i in np.arange(0,nbLambda,1):
#        for j in np.arange(0,nbLambda,1):
#            f_lambda = f_lambda_ext(f_d, lambdaTab[i])
#            f_length = f_length_ext(f_lambda)
#            M[i,j] = pressure_ext(Ps, f_length)
#    return M

#def matrix_map_pressure(f_lambda_ext, f_lambda_int)

def plot_pressures_around_airfoil_mat(file_name, size_matrix):
    """ Plots the pressures around the airfoil defined in the file file_name """
    (ex, ey, ix, iy) = spl.load_foil(file_name)
    airfoil_interp_int = spl.cubic_spline_interpolation(ix, iy)
    airfoil_interp_ext = spl.cubic_spline_interpolation(ex, ey)
    airfoil_interp_int_d = spl.cubic_spline_interpolation_derivative(ix, iy)
    airfoil_interp_ext_d = spl.cubic_spline_interpolation_derivative(ex, ey)

    hmax = max(ey)
    hmin = min(iy)

    nb_lambda = 20.
    lambda_max = 1.
    f_lambda_ext = range(0)
    f_lambda_int = range(0)
    f_lambda_ext_d = range(0)
    f_lambda_int_d = range(0)

    t = np.arange(0.0, 1.0, 0.01)

    print "Pas :", lambda_max/nb_lambda

    for l in np.arange(0., lambda_max, (lambda_max / nb_lambda)):
        f_lambda_ext.append(f_lambda(airfoil_interp_ext, l, hmax))
        f_lambda_int.append(f_lambda(airfoil_interp_int, l, hmin))
        f_lambda_ext_d.append(f_lambda_d(airfoil_interp_ext_d, l))
        f_lambda_int_d.append(f_lambda_d(airfoil_interp_int_d, l))
        

    N = size_matrix
    #Matrice...
    M = np.zeros([N,N])

    for x in np.arange(0.,N):
        x_real = min(ex) + x * (max(ex) - min(ex))/N
        for y in np.arange(0.,N):
            y_real = min(ey) + y * (max(ey) - min(ey))/N
            # y_real is inside the airfoil
            if (y_real < airfoil_interp_ext(x_real)) and (y_real > airfoil_interp_int(x_real)):
                M[x,y] = 0;
            elif y_real >= airfoil_interp_int(x_real) :
                i = 0
                while (y_real >= f_lambda_ext[i](x_real)):
                    i += 1
                M[x, y] = pressure_on_curve(0,length_of_curve(f_lambda_ext_d[i], it.simpson, min(ex), max(ex), 100.))
            else : # Else y_real <= airfoil_interp_int(x_real)
                i = 0
                while (y_real <= f_lambda_int[i](x_real)):
                    i += 1
                M[x, y] = pressure_on_curve(0,length_of_curve(f_lambda_int_d[i], it.simpson, min(ex), max(ex), 100.))

    return M

            
def matrix_to_png(matrix, name):
    """ Convertit une matrice matrix en une image png """
    mp.clf()
    mp.imshow(matrix, cmap='hot', interpolation='bilinear')
    mp.savefig(name)
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


curves_pressure_representation(int_interp_fun, int_interp_fun_d, ext_interp_fun, ext_interp_fun_d, it.simpson, 0, 0.1, 50)

#airflow_model(int_interp_fun, ext_interp_fun)

Mat = plot_pressures_around_airfoil_mat("fx63145.dat", 10)
matrix_to_png(Mat, "test_pressure.png")



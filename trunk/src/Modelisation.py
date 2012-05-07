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
# AUXILIARY FUNCTION :                           #
#------------------------------------------------#

def deriv(f):
    "Computes the numerical derivative of a function."
    def df(x, h=0.1e-5):
        return lambda(x): ( f(x+h/2) - f(x-h/2) )/h
    return df

z = lambda(x):x**3
print deriv(z)(1)

#------------------------------------------------#
# EXTRADOS :                                     #
#------------------------------------------------#
  
def f_lambda_e(f, l):
    hmax = max(ey)
    return 3*l*hmax + (1-l)*f(x)

# def f_length_e(f_d): #returns a list (result, error)
#     g = lambda(x): np.sqrt(1+f_d(x)**2)
#     t = sp.quadrature(g, 0, max(ex), args=(), tol=1.5e-8, maxiter=50)
#     return t[0]

def P_e(Ps,r,f, l):
    "calculates Pressure"

    "Derivation calculation"
    f_d = deriv(f)

    "Length calculation"
    #g = np.sqrt(1+f_d**2)#Probleme de dérivation...
    calc = sp.quadrature(np.sqrt(1+f_d**2), 0, max(ex), args=(), tol=1.5e-8, maxiter=50) 

    "Speed calculation"
    v = calc/1
    
    return Ps + r*v**2

# def P_e(Ps,r,f_length_e):
#    v = f_length_e/1
#    return Ps + r*v**2

#------------------------------------------------#
# INTRADOS :                                     #
#------------------------------------------------#
  
def f_lambda_i(f, l):
    hmin = min(iy)
    return lambda(x):3*l*hmin + (1-l)*f(x)

def f_length_i(f_d): #returns a list (result, error)
    g = lambda(x): np.sqrt(1+f_d(x)**2)
    t = sp.quadrature(g, 0, max(ix), args=(), tol=1.5e-8, maxiter=50)
    return t[0]

  
# def P_i(Ps,r,f_length_i):
#     v = f_length_i/1
#     return Ps + r*v**2

 # def P_i(Ps, r, f_d):
 #     g = lambda(x): np.sqrt(1+f_d(x)**2)
 #     calc = sp.quadrature(g, 0, max(ix), args=(), tol=1.5e-8, maxiter=50)
 #     v = f_length_i[0]/1
 #     return Ps + r*v**2

# def P_i(Ps, r, f):
#     f_l = lambda(x): 3*l*hmax + (1-l)*f(x)
#     f_l_d = lambda(x): derivative(f_l, x, 0.0001)

#     g = lambda(x): np.sqrt(1+f_d_l(x)**2)
#     f_length = sp.quadrature(g, 0, max(ix), args=(), tol=1.5e-8, maxiter=50)
    



#------------------------------------------------#
# DRAW :                                         #
#------------------------------------------------#
     
def curve(f,n,t):
    curve_values = range(n)
    for i in range(n):
        curve_values[i] = f(t[i])
    return curve_values
    

def airflow_model(f_int, f_ext):

    t = np.arange(0.0, 1.0, 0.01)
    lambdaTab = np.arange(0.0,1.0,0.1)  
    nbLambda = len(lambdaTab)
    
    for i in range(nbLambda):
        g = f_lambda_i(f_int, lambdaTab[i])
        f_lambda_values_int = curve(g,len(t),t)
        mp.plot(t,f_lambda_values_int,linewidth=1.0)
    for j in range(nbLambda):
        h = f_lambda_e(f_ext, lambdaTab[j])
        f_lambda_values_ext = curve(h,len(t),t)
        mp.plot(t,f_lambda_values_ext,linewidth=1.0)

    mp.axis([min(ex), max(max(ex),max(ix)), 3*min(iy), 3*max(ey)])
    mp.title('Laminar airflow of the fx63145')
    mp.savefig("test.png")

    return 1

def matrix_to_png(matrix, name):
    """ Convertit une matrice matrix en une image png """
    mp.imshow(matrix, cmap='hot', interpolation='bilinear')
    mp.savefig(name)

#------------------------------------------------#
# TESTS :                                        #
#------------------------------------------------#
 

#fe_a = f_lambda_e(ext_interp_fun, 0.1)
#fe_b = f_lambda_e(ext_interp_fun, 1.0)

# print f_length_e(derivative(fe_a))
# print f_length_e(derivative(fe_b))
print P_e(1, 1, ext_interp_fun, 0.1)

# print f_lambda_i(int_interp_fun, 0.1)
# print f_length_i(int_interp_fun_d, 0.1)
#print P_i(1, 1, f_length_i)

airflow_model(int_interp_fun, ext_interp_fun)


import numpy as np
# import scipy.integrate as sp
# import splines as spl
import pylab as mp

# h = lambda(x): x
# g = sp.quadrature(h, 0, 3, args=(), tol=1.5e-8, maxiter=50)
# print g[0]

def matrix_to_png(matrix, name):
    """ Convertit une matrice matrix en une image png """
    mp.imshow(matrix, cmap='hot', interpolation='bilinear')
    mp.savefig(name)
    
def Matrix_square(size):
    M = np.zeros([size,size])
    for i in range(size):
        for j in range(size):
            if i == 1:
                M[i,j]=1
            elif i == 2:
                M[i,j]=5
            elif i == 3 :
                M[i,j]=100
            elif i == 4:
                M[i,j]=5000
            elif i == 5:
                M[i,j]=4
            else :
                M[i,j]=200
    return M


print matrix_to_png(Matrix_square(20), "toto.png")

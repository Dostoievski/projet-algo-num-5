import numpy as np
import scipy.integrate as sp
import splines as spl
import pylab as mp


#--------------Five ways to integrate a function--------------------#

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

#-------------------- Function --------------------#
##################################################
############ SPLINE  FUNCTION ####################
##################################################

f = lambda(x):x**3


##################################################
############ AUXILIAIRY FUNCTION #################
##################################################

def derivative(func, x, dx):
    return (f(x+dx)-f(x))/dx


##################################################
#################### EXTRADOS ####################
##################################################
  

def f_lambda_e(f, l):
    hmax = max(spl.ey)
    return lambda(x):3*l*hmax + (1-l)*f(x)

def f_length_e(f, l): #returns a list (result, error)
    g = lambda(x):sqrt(1+derivate(f_lambda_e(f, l), x, 0.0001)**2)
    return sp.quadrature(g, 0, max(spl.ex), args=(), tol=1.5e-8, maxiter=50)

  
def P_e(Ps,r,f_length_e):
    v = f_length_e/1
    return Ps + r*v**2

##################################################
#################### INTRADOS ####################
##################################################

def f_lambda_i(f, l):
    hmax = min(spl.iy)
    return lambda(x):3*l*hmax + (1-l)*f(x)

def f_length_i(f, l): #returns a list (result, error)
    g = lambda(x):sqrt(1+derivate(f_lambda_i(f, l), x, 0.0001)**2)
    return sp.quadrature(g, 0, max(spl.ex), args=(), tol=1.5e-8, maxiter=50)


  
def P(Ps, r, f_length_i):
    v = f_length_i/1
    return Ps + r*v**2


##################################################
#################### DRAW ########################
##################################################

   
def curve(f,n,t):
    curve_values = range(n)
    for i in range(n):
        curve_values[i] = f(t[i])
    return curve_values
    

def airflow_model(f):

    t = np.arange(-1.0, 1.0, 0.01)
    lambdaTab = np.arange(0.0,1.0,0.1)  
    nbLambda = len(lambdaTab)
    for i in range(nbLambda):
        g = f_lambda_i(f,lambdaTab[i])
        f_lambda_values = curve(g,len(t),t)
        mp.plot(t,f_lambda_values,linewidth=1.0)
    for j in range(nbLambda):
        h = f_lambda_e(f,lambdaTab[i])
        f_lambda_values2 = curve(h,len(t),t)
        mp.plot(t,f_lambda_values2,linewidth=1.0)


    mp.title('laminar airflow ')
    mp.savefig("test.png")

    return 1


##################################################
#################### TESTS #######################
##################################################

airflow_model(f)


h = lambda(x):x**3
print sp.quadrature(h, 0, 3, args=(), tol=1.5e-8, maxiter=50)

import numpy as np
import scipy as sp
import pylab as mp



#------------------------------------------------#
# FUNCTIONS :                                    #
#------------------------------------------------#


def rectangle_method_left(f, a, b, n):
    h = (b - a)/n
    add = 0
    for k in np.arange(0,n):
        add += f(a + k*h)
    add *= h
    return add


    
def rectangle_method_right(f, a, b, n):
    h = (b - a)/n
    add = 0
    for k in np.arange(1,n+1):
        add += f(a + k*h)
    add *= h
    return add


def rectangle_method_middle(f, a, b, n):
    h = (b - a)/n
    add = 0
    for k in np.arange(0,n):
        add += f(a + k*h + h/2)
    add *= h
    return add


def trapez_method(f, a, b, n):
    h = (b - a)/n
    add = 0
    for k in np.arange(1,n):
        add += f(a + k*h)
    add += (f(a)+f(b))/2
    add *= h
    return add

def simpson(f, a, b, n):
    h = (b - a)/n
    add = 0.
    for k in np.arange(1,n):
        add += (1./6.)*f(a+(k-1.)*h) + (2./3.)*f(a+(k-0.5)*h) + (1./6.)*(f(a+k*h))
    add *= h
    return add

#------------------------------------------------#
# TESTS :                                        #
#------------------------------------------------#

f = lambda(x): np.sqrt(1-np.sin(x)**2)

print rectangle_method_left(f, 0., 1., 10.)
print rectangle_method_right(f, 0., 1., 10.)
print rectangle_method_middle(f, 0., 1., 10.)
print trapez_method(f, 0., 1., 10.)
print simpson(f, 0., 1., 10.)


#------------------------------------------------#
# DRAW :                                         #
#------------------------------------------------#

mp.clf()

def curve(f,n,t):
    curve_values = np.arange(n)
    for i in np.arange(n):
        curve_values[i] = f(t[i])
    return curve_values

n=10.
x = np.arange(0.0, 1. , (1.-0.)/n)
y = curve(f, n, x)
mp.plot(x, y)
mp.title('Integration comparison')
mp.savefig("integration.png")


def trace_integrals(f, a, b, nmax, l_fcts):
    n = np.arange(0, nmax, 1)
    tab =
    
    for k in range( len(l_fcts)):

# -*- coding:utf-8 -*-

import numpy as np
import scipy as sp
import pylab as mp



#------------------------------------------------#
# FUNCTIONS WITH NUMBER OF POINTS :              #
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
    for k in np.arange(1,n+1):
        add += (1./6.)*f(a+(k-1.)*h) + (2./3.)*f(a+(k-0.5)*h) + (1./6.)*(f(a+k*h))
    add *= h
    return add

#------------------------------------------------#
# FUNCTIONS WITH PRECISION :                     #
#------------------------------------------------#

def rectangle_method_left_eps(f, a, b, e):
    
    n = 1
    h = (b - a)
    I = f(a) * h
    
    I2n = I+1.
    
    while (abs(I2n - I) > e):
        mem = I2n
        I2n = 0.
        for k in np.arange(0.,n):
            I2n += f(a + k*h + h/2.)
        h = h/2.
        I2n *= h
        I2n += I/2.
        I = mem
        n *= 2

    print "N =", n
    return I2n

def rectangle_method_right_eps(f, a, b, e):
    
    n = 1
    h = (b - a)
    I = f(a) * h
    
    I2n = I+1.
    
    while (abs(I2n - I) > e):
        mem = I2n
        I2n = 0.
        for k in np.arange(1,n+1):
            I2n += f(a + k*h + h/2.)
        h = h/2.
        I2n *= h
        I2n += I/2.
        I = mem
        n *= 2

    print "N =", n
    return I2n

def rectangle_method_middle_eps(f, a, b, e):
    
    n = 1
    h = (b - a)
    I = f((b+a)/2) * h
    
    I2n = I+1.
    
    while (abs(I2n - I) > e):
        mem = I2n
        I2n = 0.
        for k in np.arange(0.,n):
            I2n += f(a + k*h + h/6.) + f(a + k*h + (5.*h)/6.)
        h = h/3.
        I2n *= h
        I2n += I/3.
        I = mem
        n *= 3

    print "N =", n
    return I2n



#------------------------------------------------#
# TESTS :                                        #
#------------------------------------------------#

#f = lambda(x): np.sqrt(1-x**2)
f = lambda(x): x**5
a = -0.5
b = 1.

#real_value = np.pi/2
real_value = (b**6/6. - a**6/6.)

n = 1000.

print "Valeur théorique :", real_value
print rectangle_method_left(f, a, b, n)
print rectangle_method_right(f, a, b, n)
print rectangle_method_middle(f, a, b, n)
print trapez_method(f, a, b, n)
print simpson(f, a, b, n)

eps = 1e-2
#print "Méthode des rectangles (gauche) avec l'erreur :", rectangle_method_left_eps(f, a, b, eps)
#print "Méthode des rectangles (droite) avec l'erreur :", rectangle_method_right_eps(f, a, b, eps)
#print "Méthode des rectangles (milieu) avec l'erreur :", rectangle_method_middle_eps(f, a, b, eps)


#------------------------------------------------#
# DRAW :                                         #
#------------------------------------------------#

mp.clf()

def curve(f,n,t):
    curve_values = np.arange(0.,n)
    for i in np.arange(n):
        curve_values[i] = f(t[i])
        #print "curve :", curve_values[i], "f(t) :", f(t[i])
    return curve_values

n=10.
x = np.arange(0.0, 1. , (1.-0.)/n)
y = curve(f, n, x)
mp.plot(x, y)
mp.title('function trace')
mp.savefig("function.png")


def trace_integrals(f, a, b, nmax, l_fcts, l_labels):
    n = np.arange(0., nmax, 1.)
    m = len(l_fcts)
    integrals = range(0)
        
    for k in np.arange(m):
        integrals.append(np.arange(0, nmax, 1))
        for j in np.arange(nmax):
                integrals[k][j] = l_fcts[k](f, a, b, n[j])
                
    mp.clf()
    for l in range(m):
        mp.plot(n, integrals[l], label=l_labels[l])
    mp.legend(loc='best')
    mp.title('Integration method trace')
    mp.savefig("integration.png")

    mp.clf()


def trace_errors(f, a, b, nmax, l_fcts, l_labels, real_value):
    n = np.arange(2., nmax, 1.)
    m = len(l_fcts)
    errors = range(0)
    
    for k in np.arange(m):
        errors.append(np.arange(nmax-2))
        for j in np.arange(nmax-2):
            errors[k][j] = abs(l_fcts[k](f, a, b, n[j]) - real_value)
    mp.clf()
    for l in range(m):
        mp.plot(n, errors[l], label=l_labels[l])
    
    mp.legend(loc='best')
    mp.title('Error depending on integration methods')
    mp.savefig("errors_integration.png")

    mp.legend(loc='lower left')
    mp.xscale('log')
    mp.yscale('log')
    mp.savefig("log_errors_integration.png")

    mp.clf()



    
l_fcts = (rectangle_method_left,
          rectangle_method_right,
          rectangle_method_middle,
          trapez_method,
          simpson)

l_labels = ("rectangle_method_left",
            "rectangle_method_right",
            "rectangle_method_middle",
            "trapez_method",
            "simpson")


n = 300.
trace_integrals(f, 0., 1., 100., l_fcts, l_labels)
trace_errors(f, a, b, n, l_fcts, l_labels, real_value)

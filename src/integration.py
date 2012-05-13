# -*- coding:utf-8 -*-

import numpy as np
import scipy as sp
import pylab as mp



#------------------------------------------------#
# FUNCTIONS WITH NUMBER OF POINTS :              #
#------------------------------------------------#


def rectangle_method_left(f, a, b, n):
    """returns the integral values of f between a and b with n iteration using the left rectangle method"""
    h = (b - a)/n
    add = 0
    for k in np.arange(0,n):
        add += f(a + k*h)
    add *= h
    return add


    
def rectangle_method_right(f, a, b, n):
    """returns the integral values of f between a and b with n iteration using the right rectangle method"""
    h = (b - a)/n
    add = 0
    for k in np.arange(1,n+1):
        add += f(a + k*h)
    add *= h
    return add


def rectangle_method_middle(f, a, b, n):
    """returns the integral values of f between a and b with n iteration using the middle rectangle method"""
    h = (b - a)/n
    add = 0
    for k in np.arange(0,n):
        add += f(a + k*h + h/2)
    add *= h
    return add


def trapez_method(f, a, b, n):
    """returns the integral values of f between a and b with n iteration using the trapez method"""
    h = (b - a)/n
    add = 0
    for k in np.arange(1,n):
        add += f(a + k*h)
    add += (f(a)+f(b))/2
    add *= h
    return add

def simpson(f, a, b, n):
    """returns the integral values of f between a and b with n iteration using the simpson method"""
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
    """returns the integral values of f between a and b with e precision using the left rectangle method"""
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
    """returns the integral values of f between a and b with e precision using the right rectangle method"""
    
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
    """returns the integral values of f between a and b with e precision using the middle rectangle method"""
    
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

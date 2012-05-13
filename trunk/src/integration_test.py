# -*- coding:utf-8 -*-

from integration import *

#------------------------------------------------#
# TESTS :                                        #
#------------------------------------------------#

f = lambda(x): x**5
a = -0.5
b = 1.

real_value = (b**6/6. - a**6/6.)

# Special integral, all the methods converges at the same speed with it
# f = lambda(x): np.sqrt(1-x**2)
# real_value = np.pi/2
# a = -1
# b = 1

n = 1000.


print "Theoretical Value :", real_value
print "Value for n =", n, "left rectangle method :", rectangle_method_left(f, a, b, n)
print "Value for n =", n, "right rectangle method :", rectangle_method_right(f, a, b, n)
print "Value for n =", n, "middle rectangle method :", rectangle_method_middle(f, a, b, n)
print "Value for n =", n, "trapezium method :", trapez_method(f, a, b, n)
print "Value for n =", n, "simpson :", simpson(f, a, b, n)


# These methods are very long to compute and not very precise, it may be an error in the algorithm
eps = 1e-2
#print "Left rectangle method with error (epsilon =", eps,") :", rectangle_method_left_eps(f, a, b, eps)
#print "Right rectangle method with error (epsilon =", eps,") :", rectangle_method_right_eps(f, a, b, eps)
#print "Middle rectangle method with error (epsilon =", eps,") :", rectangle_method_middle_eps(f, a, b, eps)


#------------------------------------------------#
# DRAW :                                         #
#------------------------------------------------#

mp.clf()

def curve(f,n,t):
    curve_values = np.arange(0.,n)
    for i in np.arange(n):
        curve_values[i] = f(t[i])
    return curve_values

n=10.
x = np.arange(0.0, 1. , (1.-0.)/n)
y = curve(f, n, x)
mp.plot(x, y)
mp.title('function trace')
mp.savefig("function.png")


def trace_integrals(f, a, b, nmax, l_fcts, l_labels):
    """Plots a graph comparing the convergence speed of each integral method"""
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
    mp.xlabel('n')
    mp.ylabel('integral value')
    mp.title('Integration methods comparison')
    mp.savefig("integration.png")

    mp.clf()


def trace_errors(f, a, b, nmax, l_fcts, l_labels, real_value):
    """Plots a graph comparing the evolution of the error of each integral method"""
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
    mp.xlabel('n')
    mp.ylabel('error')
    mp.title('Error depending on integration methods')
    mp.savefig("errors_integration.png")
    
    mp.xlabel('n')
    mp.ylabel('error')
    mp.legend(loc='lower left')
    mp.xscale('log')
    mp.yscale('log')
    mp.title('Error depending on integration methods (logarithmic scale)')
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
            "trapezium_method",
            "Simpson")


n = 300.
trace_integrals(f, 0., 1., 100., l_fcts, l_labels)
trace_errors(f, a, b, n, l_fcts, l_labels, real_value)

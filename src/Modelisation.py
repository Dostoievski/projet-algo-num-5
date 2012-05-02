import numpy as np
import scipy as sp
import "splines.py" as spl


#-------------------- Function --------------------#

##################################################
#################### EXTRADOS ####################
##################################################


def f_lambda_e(f, lambda, spl.ey):
    hmax = max(spl.ey)
    return 3*lambda*hmax + (1-lambda)*f

def f_length_e(f):
    return integrale(sqrt(1+derivate(f)**2))




##################################################
#################### INTRADOS ####################
##################################################

def f_lambda_i(f, lambda, spl.iy):
    hmax = min(spl.iy)
    return 3*lambda*hmax + (1-lambda)*f

def f_length_i(f):
    return integrale(sqrt(1+derivate(f)**2))



##################################################
#################### PRESSURE ####################
##################################################

    
def P(Ps,r,f_length):
    v = f_length/1
    return Ps + r*v**2


##################################################
#################### DRAW ########################
##################################################


def f_lambda(f,l,hmax):
    return lambda(x):3*l*hmax + (1-l)*f(x)

   
def curve(f,n,t):
    curve_values = range(n)
    for i in range(n):
        curve_values[i] = f(t[i])
    return curve_values
    

def airflow_model(f):

    hmax = 1
    t = np.arange(0.0, 1.0, 0.01)
    lambdaTab = np.arange(0.0,1.0,0.1)  
    nbLambda = len(lambdaTab)
    for i in range(nbLambda):
        g = f_lambda(f,lambdaTab[i],hmax)
        f_lambda_values = curve(g,len(t),t)
        mp.plot(t,f_lambda_values,linewidth=1.0)


    mp.title('laminar airflow ')
    mp.savefig("laminar_airflow.png")

    return 1



h = lambda(x):x**3

airflow_model(h)





##################################################
#################### TESTS #######################
##################################################

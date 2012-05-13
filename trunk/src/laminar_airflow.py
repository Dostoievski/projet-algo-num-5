import numpy as np;
import scipy as sp;
import matplotlib.pyplot as mp;
#import splines.py as sp;


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

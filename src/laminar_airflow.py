import numpy as np;
import scipy as sp;
import matplotlib.pyplot as mp;
import splines as spl;


def f_lambda(f,l,hmax):
    return lambda(x):3*l*hmax + (1-l)*f(x)

   
def curve(f,n,t):
    curve_values = range(n)
    for i in range(n):
        curve_values[i] = f(t[i])
    return curve_values
    

def airflow_model(file_src):
    """ Generates the paths followed by the air particles around the airfoil defined by f_int and f_ext """
    
    (ex,ey,ix,iy) = spl.load_foil(file_src)

    f_ext = spl.cubic_spline_interpolation(ex, ey)
    f_int = spl.cubic_spline_interpolation(ix, iy)
    
    t = np.arange(0.0, 1.0, 0.01)
    lambdaTab = np.arange(0.0,1.0,0.1)  
    nbLambda = len(lambdaTab)
    
    for i in range(nbLambda):
        f_lambda_int = f_lambda(f_int, lambdaTab[i], min(iy))
        f_lambda_values_int = curve(f_lambda_int,len(t),t)
        mp.plot(t,f_lambda_values_int,linewidth=1.0, color='#0000ff')
        
    for j in range(nbLambda):
        f_lambda_ext = f_lambda(f_ext, lambdaTab[j], max(ey))
        f_lambda_values_ext = curve(f_lambda_ext,len(t),t)
        mp.plot(t,f_lambda_values_ext,linewidth=1.0, color='#0000ff')

    mp.axis([min(ex), max(max(ex),max(ix)), 3*min(iy), 3*max(ey)])
    mp.xlabel('x')
    mp.ylabel('y')
    mp.title('Laminar airflow of the fx63145')
    mp.savefig("airflow.png")
    mp.clf()

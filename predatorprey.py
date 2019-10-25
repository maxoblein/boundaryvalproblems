import sys
from scipy.integrate import ode
import scipy.signal as signal
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal


def odefuncPP(X,t,parameters):
    a = parameters[0]
    b = parameters[1]
    d = parameters[2]
    dx = X[0]*(1-X[0]) - ((a*X[0]*X[1])/(d+X[0]))
    dy = b * X[1] * (1 - (X[1]/X[0]))
    dXdt = [dx,dy]
    return np.array(dXdt)

# Function to implement one step of the runge kutter fourth order method
# X1 is the x at start of step and X2 is the end with respective t1 and t2
def rk4(f,X1,t1,t2,parameters):

    h = t2-t1
    k1 = h*f(X1,t1,parameters)
    k2 = h*f(X1+(k1/2),t1+(h/2),parameters)
    k3 = h*f(X1+(k2/2),t1+(h/2),parameters)
    k4 = h*f(X1+k3,t1+h,parameters)
    X2 = X1 + ((1/6) * (k1+(2*k2)+(2*k3)+k4))
    return X2

def rk4solver(f,X0,t,parameters):
    #decide on a step size
    h =  t[-1]/(np.size(t)-1)
    t1 = t[0]
    rk4_plot = np.array(X0)
    X_rk4 = X0
    #check that rk4 does not go pat the tend specified
    while t1 < t[-1]:
        tnext = min(t[-1],t1+h)
        X_rk4 = rk4(f,X_rk4,t1,tnext,parameters)
        rk4_plot = np.vstack((rk4_plot,X_rk4))
        t1 = t1+h
        t1 =round(t1,1)

    return rk4_plot

if __name__ == '__main__':
    X0 = [0.4,0.4]
    t = np.linspace(0,500,5001)
    equaltime = []

    fig = plt.figure()
    ax = fig.add_axes([0.20, 0.20, 0.70, 0.70])

    parameters = [1,0.26,0.1]

    plot_array = rk4solver(odefuncPP,X0,t,parameters)

    peak_array = signal.find_peaks(plot_array[:,0])
    
    print(peak_array)
    print(period)
    ax.plot(t,plot_array[:,0])
    ax.plot(t,plot_array[:,1])
    plt.show()

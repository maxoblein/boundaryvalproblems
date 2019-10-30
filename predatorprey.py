import sys
from scipy.integrate import odeint
import scipy.signal as signal
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal


def odefuncPP(X,t,a,b,d):
    '''
        function to implement ode for predator prey system

    '''
    X.tolist()
    '''a = parameters[0]
    b = parameters[1]
    d = parameters[2]'''
    dx = X[0]*(1-X[0]) - ((a*X[0]*X[1])/(d+X[0]))
    dy = b * X[1] * (1 - (X[1]/X[0]))
    dXdt = [dx,dy]
    return np.array(dXdt)

# Function to implement one step of the runge kutter fourth order method
# X1 is the x at start of step and X2 is the end with respective t1 and t2
def rk4(f,X1,t1,t2,parameters):
    '''
    function to implement the runge kutta fourth order integrator

    inputs: -f ode Function
            -X1 approx solution at start of timestep
            -t1,t2 start and end of timestep
            -parameters for ode

    output: -X2 approx solution after timestep

    '''

    h = t2-t1
    k1 = h*f(X1,t1,parameters)
    k2 = h*f(X1+(k1/2),t1+(h/2),parameters)
    k3 = h*f(X1+(k2/2),t1+(h/2),parameters)
    k4 = h*f(X1+k3,t1+h,parameters)
    X2 = X1 + ((1/6) * (k1+(2*k2)+(2*k3)+k4))
    return X2

def rk4solver(f,X0,t,parameters):
    '''
    function to implement rk4 over a time period

    inputs: -f ode function
            -X0 initial conditions
            -t timeperiod
            -parameters for ode

    outputs: -the solutions after every timestep as an array

    '''
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

def find_period(sol_array):
    '''
    function to find the period of a given function

    inputs: sol_array array of the solutions to the ode as found in integrator

    outputs: the period of the given function
    '''
    peak_array, properties = signal.find_peaks(sol_array[:,0])

    peak_times = []

    for i in peak_array:
        peak_times.append(t[i])


    periods = []

    for j in range(1,len(peak_times)):
        periods.append(peak_times[j] - peak_times[j-1])

    period_array = np.array(periods)
    return np.mean(period_array)

def sol_after_period(X0,f,t,parameters):
    '''
    function to find value after a calculated time period

    inputs: -X0 initial conditions
            -f ode function
            -t timespan
            -parameters for ode

    outputs: -difference between inital conditions and the value after time period
    '''
    sol_array = rk4solver(f,X0,t,parameters)
    period = find_period(sol_array)
    index_of_period = np.argwhere(abs(t-period)< 0.05)
    if np.size(index_of_period) != 1:
        return [0,0]
    index_of_period = index_of_period[0,0]

    return abs(X0-sol_array[index_of_period,:])

def sol_after_given_period(X0_T,f,t,parameters):
    '''given a set of initial conditions and a Time guessed to be the period
    calculate the solution a this time

    inputs: -X0_T array of initial guess conditions for x,y, and period
            -f the function of the ode being analysed
            -t timeperiod to be integrated through
            -parameters of the ode

    outputs: -return solution after time period

    '''

    X0 = X0_T[0:2]
    period = X0_T[2]
    sol_array = odeint(f,X0,t,args = parameters)
    index_of_period = np.argwhere(abs(t-period)< 0.05)
    if np.size(index_of_period) != 1:
        return [0,0]
    index_of_period = index_of_period[0,0]

    return sol_array[index_of_period,:]



def optimal_start(X0,t,parameters):
    '''
    function to find optimal initial conditions after calculating a time period

    inputs: -X0 initial conditions
            - t timespan
            -parameters for ode

    outputs: -plots solution for optimal starting conditions
    '''
    fig = plt.figure()
    ax = fig.add_axes([0.20, 0.20, 0.70, 0.70])
    solution = fsolve(sol_after_period,X0,(odefuncPP,t,parameters))
    print(solution)

    plot_array = rk4solver(odefuncPP,solution,t,parameters)
    ax.plot(t,plot_array[:,0])
    ax.plot(t,plot_array[:,1])
    ax.hlines(solution[0],0,50)
    ax.hlines(solution[1],0,50)
    ax.vlines(18,0.25,0.35)
    plt.show()

def phase_condition(X0_T,f,t,parameters):
    '''
    function that implements the constraints on the ode

    inputs: -X0_T array of initial conditions and guess at period
            -f ode Function
            -t timespan
            -parameters for ode

    outputs: -phi constraints to be made to zero
    '''
    phi = np.zeros([3,1])
    phi[2,0] = X0_T[0] - 0.32
    phi[0,0] = (X0_T[0:2] - sol_after_given_period(X0_T,f,t,parameters))[0]
    phi[1,0] = (X0_T[0:2] - sol_after_given_period(X0_T,f,t,parameters))[1]
    phi = phi.flatten()
    return phi

if __name__ == '__main__':
    #X0 = [0.3,0.3]
    t = np.linspace(0,500,5001)

    fig = plt.figure()
    ax = fig.add_axes([0.20, 0.20, 0.70, 0.70])

    parameters = (1,0.26,0.1)
    X0_T = np.array([0.32,0.28,20])
    solution = fsolve(phase_condition,X0_T,(odefuncPP,t,parameters))
    print(solution)
    plot_array = odeint(odefuncPP,solution[0:2],t,args = parameters)
    ax.plot(t,plot_array[:,0])
    ax.plot(t,plot_array[:,1])
    ax.hlines(solution[0],0,50)
    ax.hlines(solution[1],0,50)
    ax.vlines(solution[2],0.2,0.4)
    plt.show()

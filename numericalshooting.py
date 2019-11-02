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
    dx = X[0]*(1-X[0]) - ((a*X[0]*X[1])/(d+X[0]))
    dy = b * X[1] * (1 - (X[1]/X[0]))
    dXdt = [dx,dy]
    return np.array(dXdt)

'''def odefuncHOPF(X,t,alpha,beta):
    ''''''
        function to implement the ode for the hopf bifurcation
    ''''''
    du1 = beta*X[0] - X[1] + alpha*X[0]*((X[0]**2) + (X[1]**2))
    du2 = X[0] + beta*X[1] + alpha*X[1]*((X[0]**2) + (X[1]**2))
    dXdt = [du1,du2]
    return np.array(dXdt)'''


def sol_after_given_period(X0_T,f,parameters):
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
    t = np.linspace(0,period)
    sol_array = odeint(f,X0,t,args = parameters)


    return sol_array[-1,:]

def phaseconditionPP(X0_T):
    return X0_T[0] - 0.32

def constraints(X0_T,f,phasecondition,parameters):
    '''
    function that implements the constraints on the ode

    inputs: -X0_T array of initial conditions and guess at period
            -f ode Function
            -t timespan
            -parameters for ode

    outputs: -phi constraints to be made to zero
    '''
    phi = np.zeros([3,1])
    phi[2,0] = phasecondition(X0_T,parameters)
    phi[0,0] = (X0_T[0:2] - sol_after_given_period(X0_T,f,parameters))[0]
    phi[1,0] = (X0_T[0:2] - sol_after_given_period(X0_T,f,parameters))[1]
    phi = phi.flatten()
    return phi

def shooting(odefunc,phasecond,parameters,X0_T):
    solution = fsolve(constraints,X0_T,(odefunc,phasecond,parameters))
    return(solution)

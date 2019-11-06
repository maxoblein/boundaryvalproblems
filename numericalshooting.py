import sys
from scipy.integrate import odeint
import scipy.signal as signal
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal

def sol_after_given_period(X0_T,f,parameters):
    '''given a set of initial conditions and a Time guessed to be the period
    calculate the solution a this time

    inputs: -X0_T array of initial guess conditions for x,y, and period
            -f the function of the ode being analysed
            -t timeperiod to be integrated through
            -parameters of the ode

    outputs: -return solution after time period

    '''

    X0 = X0_T[0:-1]
    period = X0_T[-1]
    t = np.linspace(0,period)
    sol_array = odeint(f,X0,t,args = parameters)


    return sol_array[-1,:]


def constraints(X0_T,f,phasecondition,parameters):
    '''
    function that implements the constraints on the ode

    inputs: -X0_T array of initial conditions and guess at period
            -f ode Function
            -t timespan
            -parameters for ode

    outputs: -phi constraints to be made to zero
    '''


    phi = (X0_T[0:-1] - sol_after_given_period(X0_T,f,parameters))
    phi = phi.flatten()
    phi = np.hstack((phi,phasecondition(X0_T,parameters)))

    return phi

def shooting(odefunc,phasecond,parameters,X0_T):
    solution = fsolve(constraints,X0_T,(odefunc,phasecond,parameters))
    return(solution)

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


def constraints(X0_T,f,parameters):
    '''
    function that implements the constraints on the ode

    inputs: -X0_T array of initial conditions and guess at period
            -f ode Function
            -t timespan
            -parameters for ode

    outputs: -phi constraints to be made to zero
    '''


    phi = (X0_T[0:-1] - sol_after_given_period(X0_T,f,parameters))
    t = 0

    phi = np.hstack((phi,f(X0_T[:-1],t,*parameters)[0]))

    return phi

def constraints_cont(v,f,params,dv,v_tilde):
    phi = v[:-2] - sol_after_given_period(v[:-1],f,params)
    t = 0
    phi = np.hstack((phi,f(v[:-2],t,*params)[0]))
    phi = np.hstack((phi,np.dot(v-v_tilde,dv)))
    return phi


def shooting(odefunc,parameters,X0_T):
    solution = fsolve(constraints,X0_T,(odefunc,parameters))
    return(solution)

def natural_continuation(u0,params,odefunc,vary_param = 0,steps = 100, discretisation = lambda odefunc,parameters,X0_T : X0_T,plot = False ):
    pspan = params[vary_param]
    delta = (pspan[1] - pspan[0])/steps
    p0 = pspan[0]
    params = np.delete(params,vary_param)
    params = np.insert(params,vary_param,p0)
    u0_tilde = u0
    param_list = []
    plot_list = []
    sol_list = []
    for i in range(steps):
        u0 = discretisation(odefunc,tuple(params),u0_tilde)
        tspan = np.linspace(0,u0_tilde[-1])
        sol_array = odeint(odefunc,u0_tilde[:-1],tspan,args = tuple(params))
        mag_list = []
        for i in sol_array:
            mag_list.append(np.linalg.norm(i))
        plot_list.append(max(mag_list))
        param_list.append(params[vary_param])
        sol_list.append(u0)
        u0_tilde = np.copy(u0)
        params[vary_param] += delta
    if plot == True:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(param_list,plot_list)
        plt.show()
    return np.array(param_list), np.array(sol_list)

def pseudo_continuation(u0,params,odefunc,vary_param = 0,steps = 100,discretisation = lambda odefunc,parameters,X0_T : X0_T,plot = False):
    pspan = params[vary_param]
    delta = (pspan[1] - pspan[0])/steps
    p0 = pspan[0]
    params = np.delete(params,vary_param)
    params = np.insert(params,vary_param,p0)
    u0 = discretisation(odefunc,tuple(params),u0)
    p1 = p0 + delta
    params = np.delete(params,vary_param)
    params = np.insert(params,vary_param,p1)
    u1 = discretisation(odefunc,tuple(params),u0)
    v0 = np.hstack((u0,p0))
    v1 = np.hstack((u1,p1))
    param_list = []
    plot_list = []
    sol_list = [v0,v1]
    for i in range(steps):
        dv = v1 - v0
        v2_tilde = v1 + dv
        params = np.delete(params,vary_param)
        params = np. insert(params,vary_param,v1[-1])
        tspan = np.linspace(0,v1[-2])
        sol_array = odeint(odefunc,v1[:-2],tspan,args = tuple(params))
        mag_list = []
        for i in sol_array:
            mag_list.append(np.linalg.norm(i))
        plot_list.append(max(mag_list))
        param_list.append(params[vary_param])
        v2 = fsolve(constraints_cont,v2_tilde,(odefunc,tuple(params),dv,v2_tilde))
        sol_list.append(v2)
        v0 = np.copy(v1)
        v1 = np.copy(v2)
    fig =plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(param_list,plot_list)
    plt.show()

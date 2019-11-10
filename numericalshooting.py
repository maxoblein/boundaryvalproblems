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

    phi = np.hstack((phi,f(X0_T,t,*parameters)[0]))

    return phi

def constraints_cont(v,f,phasecondition,params,dv,v_tilde,vary_param=0):
    params = np.insert(params,vary_param,v[-1])
    params = tuple(params)
    phi = v[:-2] - sol_after_given_period(v[:-1],f,params)
    phi = np.hstack((phi,phasecondition(v[:-1],params)))
    phi = np.hstack((phi,np.dot(v-v_tilde,dv)))
    return phi


def shooting(odefunc,parameters,X0_T):
    solution = fsolve(constraints,X0_T,(odefunc,parameters))
    return(solution)

def natural_continuation(u0,params,odefunc,vary_param = 0,delta = 0.01,discretisation = lambda odefunc,phasecond,parameters,X0_T : X0_T ):
    pspan = params[vary_param]
    delta = (pspan[1] - pspan[0])/100
    p0 = pspan[0]
    params[vary_param] = p0
    u0_tilde = u0
    param_list = []
    plot_list = []
    for i in range(100):
        u0 = discretisation(odefunc,tuple(params),u0_tilde)
        tspan = np.linspace(0,u0_tilde[-1])
        sol_array = odeint(odefunc,u0_tilde[:-1],tspan,args = tuple(params))
        mag_list = []
        for i in sol_array:
            mag_list.append(np.linalg.norm(i))
        plot_list.append(max(mag_list))
        param_list.append(params[vary_param])

        u0_tilde = np.copy(u0)
        params[vary_param] += delta
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(param_list,plot_list)
    plt.show()

def pseudo_continuation(u0,params,odefunc,phasecond,vary_param = 0,delta = 0.01):
    param_span = params[vary_param]
    delta = (param_span[1] - param_span[0])/100
    p0 = param_span[0]
    p1 = p0
    params = np.delete(params,vary_param)
    param_list = []
    plot_list = []
    while p1 <=2:
        params = np.insert(params,vary_param,p0)
        param_list.append(p0)
        tspan = np.linspace(0,u0[-1])
        plot_list.append(np.max(odeint(odefunc,u0[:-1],tspan,args = tuple(params))))
        u1 = shooting(odefunc,phasecond,tuple(params),u0)
        p1 = p0 + delta
        params = np.delete(params,vary_param)
        v0 = np.hstack((u0,p0))
        v1 = np.hstack((u1,p1))
        dv = v1-v0
        v2_tilde = v1 + delta*dv
        print(v2_tilde)
        solution = fsolve(constraints_cont,v2_tilde,args = (odefunc,phasecond,params,dv,v2_tilde,1))

        u0 = u1
        p0 = p1
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(param_list,plot_list)
    plt.show()

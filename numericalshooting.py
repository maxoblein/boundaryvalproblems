import sys
from scipy.integrate import odeint
import scipy.signal as signal
from scipy.optimize import fsolve
import sys
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")

def sol_after_given_period(v,f,params):
    '''given a set of initial conditions and a Time guessed to be the period
    calculate the solution a this time

    inputs: -v ndarray of initial guess conditions for x,y, and period
            -f the function of the ode being analysed
            -t timeperiod to be integrated through
            -parameters of the ode

    outputs: -return solution after time period

    '''

    X0 = v[0:-1]
    period = v[-1]
    t = np.linspace(0,period)
    sol_array = odeint(f,X0,t,args = params)


    return sol_array[-1,:]

def constraints(v,f,params,dv = None,v_tilde = None,vary_param = None,pseudo = False):
    '''
    function that implements the constraints on the ode

    inputs: -v ndarray of initial conditions and guess at period
            -f ode Function
            -t timespan
            -parameters for ode

    outputs: -phi constraints to be made to zero
    '''

    if pseudo == False:
        phi = (v[0:-1] - sol_after_given_period(v,f,params))
        t = 0

        phi = np.hstack((phi,f(v[:-1],t,*params)[0]))
    if pseudo == True:
        params = np.delete(params,vary_param)
        params = np.insert(params,vary_param,v[-1])
        params = tuple(params)
        phi = v[:-2] - sol_after_given_period(v[:-1],f,params)
        t = 0
        phi = np.hstack((phi,f(v[:-2],t,*params)[0]))
        phi = np.hstack((phi,np.dot(v-v_tilde,dv)))
    return phi

def shooting(odefunc,params,v):
    '''
    function that implements the numerical shooting method

    inputs:
            -odefunc ode Function
            -parameters for ode
            --X0_T array of initial conditions and guess at period

    outputs: -array of correct initial conditions and timeperiod
    '''
    solution = fsolve(constraints,v,(odefunc,params))
    return(solution)

def check_input(u0,odefunc,params,vary_param):
    pspan = params[vary_param]
    if np.size(params[vary_param]) != 2:
        sys.stderr.write('Incorrect index vary_param\n')
        return 1
    else:
        p0 = pspan[0]

        params = np.delete(params,vary_param)
        params = np.insert(params,vary_param,p0)
        output = odefunc(u0,0,*tuple(params))
        if len(u0)-1 == len(output):
            return None
        else:
            sys.stderr.write('Incorrect u0 dimensions\n')
            return 1

def natural_continuation(u0,params,odefunc,vary_param = 0,steps = 100, discretisation = lambda odefunc,parameters,X0_T : X0_T,plot = False ):
    '''
    function that implements natural parameter continuation

    inputs:
            -u0 initial state variables and period ndarray
            -params parameters of sysystem one is a list to be varied ndarray
            -odefunc ode to be analysed
            -vary_param index of parameter to vary default 0
            -steps number of steps default 100
            -discretisation type of method to use e.g. shooting default is none
            -plot option to show plot default no plot

    outputs: -ndarray of state variables, timeperiod and parameter values at each parameter step
    '''
    if check_input(u0,odefunc,params,vary_param) == 1:
        return False

    param_list = []
    plot_list = []
    sol_list = []
    if discretisation.__name__ == 'shooting':
        pspan = params[vary_param]
        delta = (pspan[1] - pspan[0])/steps
        p0 = pspan[0]

        params = np.delete(params,vary_param)
        params = np.insert(params,vary_param,p0)

        u0_tilde = u0

        for i in range(steps):
            u0 = discretisation(odefunc,tuple(params),u0_tilde)
            tspan = np.linspace(0,u0[-1])
            sol_array = odeint(odefunc,u0[:-1],tspan,args = tuple(params))
            mag_list = []
            for i in sol_array:
                mag_list.append(np.linalg.norm(i))
            plot_list.append(max(mag_list))
            param_list.append(params[vary_param])
            sol = np.append(u0,params[vary_param])
            sol_list.append(sol)
            u0_tilde = np.copy(u0)
            params[vary_param] += delta
        if plot == True:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(param_list,plot_list)
            ax.set_xlabel('Bifucation parameter',fontsize = 14)
            ax.set_ylabel(r'$|u|$',fontsize=14)
            ax.set_title('Branch of limit cycles',fontsize=16)
            plt.show()

    elif discretisation.__name__ == '<lambda>':
        pspan = params[vary_param]
        delta = (pspan[1] - pspan[0])/steps
        p0 = pspan[0]
        params = np.delete(params,vary_param,axis = 0)
        for i in range(steps):
            params = np.insert(params,vary_param,float(p0))
            sol = fsolve(odefunc,u0,args = (tuple(params)))
            param_list.append(p0)
            plot_list.append(sol)
            sol_list.append(sol)
            sol_list.append(params[vary_param])
            params = np.delete(params,vary_param)
            u0 = sol
            p0 += delta
        if plot == True:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(param_list,plot_list)
            ax.set_xlabel('Bifucation parameter',fontsize = 14)
            ax.set_ylabel(r'$|u|$',fontsize=14)
            ax.set_title('Branch of limit cycles',fontsize=16)
            plt.show()

    return np.array(sol_list)

def pseudo_continuation(u0,params,odefunc,vary_param = 0,steps = 100,discretisation = lambda odefunc,parameters,X0_T : X0_T,plot = False):
    '''
    function that implements pseudo arclength parameter continuation

    inputs:
            -u0 initial state variables and period ndarray
            -params parameters of sysystem one is a list to be varied ndarray
            -odefunc ode to be analysed
            -vary_param index of parameter to vary default 0
            -steps number of steps default 100
            -discretisation type of method to use e.g. shooting default is none
            -plot option to show plot default no plot

    outputs: -ndarray of state variables, timeperiod and parameter values at each parameter step
    '''
    if check_input(u0,odefunc,params,vary_param) == 1:
        return False
    param_list = []
    plot_list = []
    sol_list = []
    if discretisation.__name__ == 'shooting':
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
        sol_list.append(v0)
        sol_list.append(v1)
        for i in range(steps-1):
            dv = v1 - v0
            v2_tilde = v1 + dv
            params = np.delete(params,vary_param)
            params = np.insert(params,vary_param,v2_tilde[-1])
            v2 = fsolve(constraints,v2_tilde,(odefunc,params,dv,v2_tilde,vary_param, True))
            params = np.delete(params,vary_param)
            params = np.insert(params,vary_param,v2[-1])
            tspan = np.linspace(0,v2[-2])
            sol_array = odeint(odefunc,v2[:-2],tspan,args = tuple(params))
            mag_list = []
            for j in sol_array:
                mag_list.append(np.linalg.norm(j))
            plot_list.append(max(mag_list))
            param_list.append(params[vary_param])
            sol_list.append(v2)
            v0 = np.copy(v1)
            v1 = np.copy(v2)

        if plot == True:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(param_list,plot_list)
            ax.set_xlabel('Bifucation parameter',fontsize = 14)
            ax.set_ylabel(r'$|u|$',fontsize=14)
            ax.set_title('Branch of limit cycles',fontsize=16)
            plt.show()

    elif discretisation.__name__ == '<lambda>':
        pspan = params[vary_param]
        delta = (pspan[1] - pspan[0])/steps
        p0 = pspan[0]
        params = np.delete(params,vary_param,axis = 0)
        for i in range(steps):
            params = np.insert(params,vary_param,float(p0))
            sol = fsolve(odefunc,u0,args = (tuple(params)))
            param_list.append(p0)
            plot_list.append(sol)
            sol_list.append(sol)
            sol_list.append(params[vary_param])
            params = np.delete(params,vary_param)
            u0 = sol
            p0 += delta



        if plot == True:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(param_list,plot_list)
            ax.set_xlabel('Bifucation parameter',fontsize = 14)
            ax.set_ylabel(r'$|u|$',fontsize=14)
            ax.set_title('Branch of limit cycles',fontsize=16)
            plt.show()
    return np.array(sol_list)

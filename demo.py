import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from numericalshooting import *

def odefuncHOPF(X,t,alpha,beta):
    '''
        function to implement the ode for the hopf bifurcation
    '''
    du1 = beta*X[0] - X[1] + alpha*X[0]*((X[0]**2) + (X[1]**2))
    du2 = X[0] + beta*X[1] + alpha*X[1]*((X[0]**2) + (X[1]**2))
    dXdt = [du1,du2]
    return np.array(dXdt)

def odefuncHOPFMOD(X,t,alpha,beta):
    du1 = beta*X[0] - X[1] +X[0]*((X[0]**2) + (X[1]**2)) - X[0]*(((X[0]**2)+(X[1]**2))**2)
    du2 = X[0] + beta*X[1] + X[1]*((X[0]**2) + (X[1]**2)) - X[1]*(((X[0]**2)+(X[1]**2))**2)
    dXdt = [du1,du2]
    return np.array(dXdt)

def odefuncPP(X,t,a,b,d):
    '''
        function to implement ode for predator prey system

    '''
    X.tolist()
    dx = X[0]*(1-X[0]) - ((a*X[0]*X[1])/(d+X[0]))
    dy = b * X[1] * (1 - (X[1]/X[0]))
    dXdt = [dx,dy]
    return np.array(dXdt)

def cubic(X,t,c):
    return  X**3 - X + c

def HOPFanalytic(beta,t):
    u1 = np.sqrt(beta) * np.cos(t)
    u2 = np.sqrt(beta) * np.sin(t)
    return[u1,u2]
#shooting(odefunc,phasecond,parameters,X0_T)
if __name__ == '__main__':

    parameters = (-1,0.1)
    X0_T = [0.3,0,6.28]
    solution = shooting(odefuncHOPF,parameters,X0_T)
    analytic_sol = HOPFanalytic(parameters[1],solution[2])
    if np.isclose(solution[0:-1],analytic_sol,atol = 1e-5).all() == True:
        print('Test passed')
    else:
        print('Test failed')


    params = np.array([-1,[2,-1]])
    #sol = natural_continuation([0.4,0.4,20],params,odefuncPP,vary_param=1,discretisation = shooting,plot = True)
    sol = pseudo_continuation([0.3,0,6.3],params,cubic,vary_param=1,discretisation = shooting,plot = True)

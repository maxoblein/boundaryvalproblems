import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import sys
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

#shooting(odefunc,phasecond,parameters,X0_T)
if __name__ == '__main__':

    if sys.argv[1] == 'natural':
        params = np.array([-1,[2,0]])
        sol = natural_continuation([0.3,0,6.3],params,odefuncHOPF,vary_param=1,discretisation = shooting,plot = True)
    if sys.argv[1] == 'pseudo':
        pseudo_continuation([0.3,0,6.3],np.array([-1,[2,-1]]),odefuncHOPFMOD,vary_param=1,discretisation = shooting,plot =True)

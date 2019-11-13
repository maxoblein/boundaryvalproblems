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

def test_shooting_output():
    solution = shooting(odefuncHOPF,(-1,0.1),[0.3,0,6.28])
    analytic_sol = HOPFanalytic((-1,0.1)[1],solution[2])
    if np.isclose(solution[0:-1],analytic_sol,atol = 1e-5).all() == True:
        print('Shooting output test passed')
    else:
        print('Shooting output test failed')
    return(None)

def test_natural_continuation_otuput():
    sol = natural_continuation([0.3,0,6.3],np.array([-1,[2,0]]),odefuncHOPF,vary_param=1,discretisation = shooting)
    y_list= []
    x_list = []
    for i in sol:
        ode_sol = odeint(odefuncHOPF,i[:-2],np.linspace(0,i[-2]),args = (-1,i[-1]))
        mag_list = []
        for k in ode_sol:
            mag_list.append(np.linalg.norm(k))
        x_list.append(i[-1])
        y_list.append(max(mag_list))
    test_y_list = []
    for j in x_list:
        test_y_list.append(np.sqrt(j))
    if np.isclose(np.array(y_list),np.array(test_y_list),atol = 1e-5).all() == True:
        print('Natural continuation output test passed')
    else:
        print('Natural continuation output test failed')
    return None

def test_incorrect_dimensions():
    if natural_continuation([0.3,6.3],np.array([-1,[2,0]]),odefuncHOPF,vary_param=1,discretisation = shooting) == False:
        print('Incorrect input test passed')
    else:
        print('Incorrect input test failed')

#shooting(odefunc,phasecond,parameters,X0_T)
if __name__ == '__main__':

    test_shooting_output()
    test_natural_continuation_otuput()
    test_incorrect_dimensions()

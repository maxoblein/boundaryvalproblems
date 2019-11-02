import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from numericalshooting import sol_after_given_period, constraints, shooting

def odefuncHOPF(X,t,alpha,beta):
    '''
        function to implement the ode for the hopf bifurcation
    '''
    du1 = beta*X[0] - X[1] + alpha*X[0]*((X[0]**2) + (X[1]**2))
    du2 = X[0] + beta*X[1] + alpha*X[1]*((X[0]**2) + (X[1]**2))
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

def phaseconditionPP(X0_T,parameters):
    return X0_T[0] - 0.32

def phaseconditionHOPF(X0_T,parameters):
    alpha, beta = parameters
    X = X0_T[0:2]
    return beta*X[0] - X[1] + alpha*X[0]*((X[0]**2) + (X[1]**2))

#shooting(odefunc,phasecond,parameters,X0_T)
solution = shooting(odefuncHOPF,phaseconditionHOPF,(-1,0.8),[0.32,0.2,10])
print(solution)
t = np.linspace(0,solution[2])
parameters = (-1,0.8)
plot_array = odeint(odefuncHOPF,solution[0:2],t, args = parameters)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(t,plot_array[:,0])
ax.plot(t,plot_array[:,1])
ax.hlines(solution[0],0,solution[2])
ax.hlines(solution[1],0,solution[2])
plt.show()

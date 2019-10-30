import numpy as np

def odefuncHOPF(X,t,alpha,beta):
    '''
        function to implement the ode for the hopf bifurcation
    '''
    du1 = beta*X[0] - X[1] + alpha*X[0]*((X[0]**2) + (X[1]**2))
    du2 = X[0] + beta*X[1] + alpha*X[1]*((X[0]**2) + (X[1]**2))
    dXdt = [du1,du2]
    return np.array(dXdt)

def odefuncPP(X,t,parameters):
    a = parameters[0]
    b = parameters[1]
    d = parameters[2]
    dx = X[0]*(1-X[0]) - ((a*X[0]*X[1])/(d+X[0]))
    dy = b * X[1] * (1 - (X[1]/X[0]))
    dXdt = [dx dy]
    return dXdt

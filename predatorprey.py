def odefuncPP(X,t,parameters):
    a = parameters[0]
    b = parameters[1]
    d = parameters[2]
    dx = X[0]*(1-X[0]) - ((a*X[0]*X[1])/(d+X[0]))
    dy = b * X[1] * (1 - (X[1]/X[0]))
    dXdt = [dx dy]
    return dXdt

# Function to implement one step of the runge kutter fourth order method
# X1 is the x at start of step and X2 is the end with respective t1 and t2
def rk4(f,X1,t1,t2):

    h = t2-t1
    k1 = h*f(X1,t1)
    k2 = h*f(X1+(k1/2),t1+(h/2))
    k3 = h*f(X1+(k2/2),t1+(h/2))
    k4 = h*f(X1+k3,t1+h)
    X2 = X1 + ((1/6) * (k1+(2*k2)+(2*k3)+k4))

    return X2

def rk4solver(f,X0,t):
    #decide on a step size
    h = 0.1
    t1 = t[0]
    X_rk4 = X0
    #check that rk4 does not go pat the tend specified
    while t1 <= t[-1]:
        tnext = min(t[-1],t1+h)
        X_rk4 = rk4(odefunc,X_rk4,t,tnext)
        t1 = t1+h
    return X_rk4

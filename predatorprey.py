import sys
from scipy.integrate import ode
import scipy.signal as signal
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import numpy as np
from decimal import Decimal


def odefuncPP(X,t,parameters):
    X.tolist()
    a = parameters[0]
    b = parameters[1]
    d = parameters[2]
    dx = X[0]*(1-X[0]) - ((a*X[0]*X[1])/(d+X[0]))
    dy = b * X[1] * (1 - (X[1]/X[0]))
    dXdt = [dx,dy]
    return np.array(dXdt)

# Function to implement one step of the runge kutter fourth order method
# X1 is the x at start of step and X2 is the end with respective t1 and t2
def rk4(f,X1,t1,t2,parameters):

    h = t2-t1
    k1 = h*f(X1,t1,parameters)
    k2 = h*f(X1+(k1/2),t1+(h/2),parameters)
    k3 = h*f(X1+(k2/2),t1+(h/2),parameters)
    k4 = h*f(X1+k3,t1+h,parameters)
    X2 = X1 + ((1/6) * (k1+(2*k2)+(2*k3)+k4))
    return X2

def rk4solver(f,X0,t,parameters):
    #decide on a step size
    h =  t[-1]/(np.size(t)-1)
    t1 = t[0]
    rk4_plot = np.array(X0)
    X_rk4 = X0
    #check that rk4 does not go pat the tend specified
    while t1 < t[-1]:
        tnext = min(t[-1],t1+h)
        X_rk4 = rk4(f,X_rk4,t1,tnext,parameters)
        rk4_plot = np.vstack((rk4_plot,X_rk4))
        t1 = t1+h
        t1 =round(t1,1)

    return rk4_plot

def find_period(sol_array):
    peak_array, properties = signal.find_peaks(sol_array[:,0])

    peak_times = []

    for i in peak_array:
        peak_times.append(t[i])


    periods = []

    for j in range(1,len(peak_times)):
        periods.append(peak_times[j] - peak_times[j-1])

    period_array = np.array(periods)
    return np.mean(period_array)

def sol_after_period(X0,f,t,parameters):
    sol_array = rk4solver(odefuncPP,X0,t,parameters)
    period = find_period(sol_array)
    index_of_period = np.argwhere(abs(t-period)< 0.05)
    if np.size(index_of_period) != 1:
        return [0,0]
    index_of_period = index_of_period[0,0]


    return abs(X0-sol_array[index_of_period,:])



def plot_peaks(t,parameters):
    mean_x_diff = []
    mean_y_diff = []
    X0_guess = np.linspace(0.25,0.4)
    for i in X0_guess:
        X0 = [i,i]
        plot_array = rk4solver(odefuncPP,X0,t,parameters)
        peak_array_x, properties = signal.find_peaks(plot_array[:,0])
        peak_array_y, properties = signal.find_peaks(plot_array[:,0])
        x_peaks = []
        y_peaks = []
        for j in peak_array_x:
            x_peaks.append(plot_array[j,0])

        for k in peak_array_y:
            y_peaks.append(plot_array[k,1])

        period = find_period(plot_array)

        x_diff = []
        y_diff = []
        for m in range(1,len(x_peaks)):
            x_diff.append(x_peaks[m] - x_peaks[m-1])

        for n in range(1,len(y_peaks)):
            y_diff.append(y_peaks[n]-y_peaks[n-1])

        x_diff_array = np.array(x_diff)
        y_diff_array = np.array(y_diff)
        mean_x_diff.append(np.mean(x_diff_array))
        mean_y_diff.append(np.mean(y_diff_array))


    ax.plot(X0_guess,mean_x_diff,label = 'xdiff')
    ax.plot(X0_guess,mean_y_diff,label = 'ydiff')
    ax.legend()
    ax.hlines(0,0.25,0.4)
    plt.show()

if __name__ == '__main__':
    #X0 = [0.3,0.3]
    t = np.linspace(0,50,501)

    fig = plt.figure()
    ax = fig.add_axes([0.20, 0.20, 0.70, 0.70])

    parameters = [1,0.26,0.1]
    X0 = np.array([0.4,0.4])
    solution = fsolve(sol_after_period,X0,(odefuncPP,t,parameters))
    print(solution)

    plot_array = rk4solver(odefuncPP,solution,t,parameters)
    ax.plot(t,plot_array[:,0])
    ax.plot(t,plot_array[:,1])
    ax.hlines(solution[0],0,50)
    ax.hlines(solution[1],0,50)
    ax.vlines(18,0.25,0.35)
    plt.show()

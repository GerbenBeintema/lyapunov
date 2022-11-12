import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
import scipy
from scipy import linalg

def _lyapunov(f, dfdx, x0, tmax, tmin=0, dt=None, intercept=0 , lamb0 = 0):
    nx = len(x0)
    Ilamb0 = np.eye(nx)*lamb0
    def dxadt(t, xa, lamb0 = 0):
        x = xa[:nx]
        a = xa[nx:].reshape((nx, nx))
        dxdt = f(t, x)
        dadt = (np.array(dfdx(t, x)) - Ilamb0)@a
        return np.concatenate([dxdt, dadt.flat])
    
    a0 = np.eye(nx)*np.exp(-intercept)
    xa0flat = np.concatenate((np.array(x0), a0.flat))
    t_eval = None if dt==None else np.linspace(tmin,tmax, int((tmax-tmin)/dt)+1)
    
    #compute:
    sim = integrate.solve_ivp(lambda t, xa: dxadt(t, xa, lamb0=lamb0), [tmin,tmax], xa0flat, t_eval=t_eval)
    
    #post process:
    t = sim.t.T
    a = np.transpose(sim.y[nx:].reshape(nx, nx, -1),[2,0,1])
    x = sim.y[:nx].T
    logmaxeigs = np.log(np.max(abs(np.real(np.linalg.eigvals(a))),axis=1)) + lamb0*t + intercept
    lamb0, intercept = np.polyfit(t,logmaxeigs, 1)
    return t, x, a, logmaxeigs, sim, lamb0, intercept

def lyapunov(f, dfdx, x0, tmax, tmin=0, dt=None, intercept=0, lamb0=0, nsteps=15, full_return=False, plot=False, verbose=False):
    #presimulate until tmin
    if tmin>0:
        sim = integrate.solve_ivp(f, [0,tmin], x0)
        x0 = sim.y.T[-1]
    
    # simulate in itteratively such that it becomes never unstable
    tdiff = tmax - tmin
    assert nsteps>0
    for n in range(nsteps-1,-1,-1):
        tmax_now = tmin + tdiff/2**n
        t, x, a, logmaxeigs, sim, lamb0, intercept = \
               _lyapunov(f, dfdx, x0, tmax_now, tmin=tmin, dt=dt, intercept=intercept , lamb0 = lamb0)
        if verbose: print(n, lamb0, tdiff/2**n)
        
    if plot:
        plt.plot(t,logmaxeigs,label='Log largest eigen value')
        plt.plot(t,t*lamb0+intercept,label=f'fit {lamb0:.3f}*t+{intercept:.2f}')
        plt.legend()
        plt.title(f'$\\lambda_0$ = {lamb0:.5f}')
        plt.grid()
        plt.xlabel('time')
        plt.ylabel('log max eigen value of a')
        plt.show()
    
    if full_return:
        return t, x, a, logmaxeigs, sim, lamb0, intercept
    else:
        return lamb0

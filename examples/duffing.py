import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
import scipy
from scipy import linalg

from lyapunov import lyapunov


class Duffing:
    def __init__(self, delta, alpha, beta, gamma, omega):
        self.delta = delta
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.omega = omega
    
    def f(self, t, x):
        x,v = x
        dvdt = -self.delta*v - self.alpha*x - self.beta*x**3 + self.gamma*np.cos(self.omega*t)
        dxdt = v
        return [dxdt, dvdt]
    def dfdx(self, t, x):
        x,v = x
        dxdtdx = [0, 1]
        dvdtdx = [-self.alpha-3*self.beta*x**2, -self.delta]
        return [dxdtdx, dvdtdx]


sys = Duffing(delta=0.02, beta=0, alpha=1, gamma=8, omega=0.5)
lambs = []
from tqdm.auto import tqdm
for x0 in tqdm(np.random.uniform(1.3,1.5,10)):
    t, x, a, logmaxeigs, sim, lamb0, intercept = lyapunov(sys.f, sys.dfdx, [x0, 0], tmin=200, tmax=200+10000, lamb0=0.046, nsteps=10, plot=True,full_return=True)
    lambs.append(lamb0)

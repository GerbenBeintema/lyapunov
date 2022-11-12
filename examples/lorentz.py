import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
import scipy
from scipy import linalg

from lyapunov import lyapunov

class LA_system():
    def __init__(self, sigma=10, beta=8/3, rho=28):
        self.sigma = sigma
        self.beta = beta
        self.rho = rho

    def f(self, t, x):
        x,y,z = x
        dxdt = self.sigma*(y - x)
        dydt = x*(self.rho - z) - y
        dzdt = x*y - self.beta*z
        return np.array([dxdt, dydt, dzdt])

    def dfdx(self, t, x):
        x,y,z = x
        dxdtdx = [-self.sigma,  self.sigma, 0    ]
        dydtdx = [self.rho - z, -1,    -x   ]
        dzdtdx = [y,       x,     -self.beta]
        return np.array([dxdtdx,dydtdx,dzdtdx])


if False:
    sigma = 10
    beta = 8/3
    rho = 28
    lamb_real = 0.90566
elif True:
    sigma = 16
    rho = 45.92
    beta = 4
    lamb_real = 1.50255
else:
    sigma = 16
    rho = 40
    beta = 4
    lamb_real = 1.37446

la = LA_system(sigma=sigma, rho=rho, beta=beta)

for dt in [0.05, 0.01, 0.001]:
    np.random.seed(2)
    x0_3s = np.random.uniform(1,10,size=7)
    from tqdm.auto import tqdm

    lambs = [lyapunov(la.f, la.dfdx, x0=[ 1.83763056,  3.33011506, x0_3], tmin=10, tmax=200, dt=dt, nsteps=5) for x0_3 in tqdm(x0_3s)]

    print(f'##### dt={dt} ######')
    std_est = np.std(lambs)/(len(lambs)-1)**0.5
    mean_est = np.mean(lambs)
    print('est:',mean_est, 'pm', std_est)
    print('real:',lamb_real)
    diff = (mean_est - lamb_real)
    print('difference:', diff, 'scaled difference:',diff/std_est, '<-- should be smaller then 1')

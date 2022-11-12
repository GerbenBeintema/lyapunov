
## lyapunov: An easy-going maximum lyapunov exponent computer for continuous systems
 
Define systems where the state derivative and Jacobian of the state-derivative are defined. 

```python
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

la = LA_system(sigma=10, rho=8/3, beta=28)
```

Then we can compute lyapunov exponent as;

```python
from lyapunov import lyapunov

lamb = lyapunov(la.f, la.dfdx, x0=[ 1.83763056,  3.33011506, 10], tmax=200)
print('computed lyapunov exponent:',lamb)
print('Actual lyapunov exponent:',0.90566)
```

You can also get more details by doing the following options

```python
t, x, a, logmaxeigs, sim, lamb0, intercept = \
    lyapunov(la.f, la.dfdx, x0=[ 1.83763056,  3.33011506, 10], tmax=200, full_return=True, plot=True)
```

Warning: I have verified that this code can work decently, but I have noticed that the computations are a little biased (about 1% of error).

todo's:

 * Function Docs
 * Method Docs
 * More validation
 * More examples
 * Discrete time lyapunov exponent
 * Lyapunov spectrum
 * Lyapunov Dimension
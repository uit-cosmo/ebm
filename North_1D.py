import numpy as np 
from scipy import sparse
from scipy.special import eval_legendre
from scipy.integrate import odeint
from matplotlib import pyplot as plt
import sys
import matplotlib.animation as animation


# Constants
A_ = 211.2 - 18.
B = 1/0.32
D = 0.38
s2 = -0.482
Tc = -10
b0 = 0.38
a0 = 0.697
a2 = -0.0779
Q = 334.4
c = 13.2
tid = 50
tmax = 131


""" 
Laplacian Construction
"""
n = 500 
Delta_x = 1/n 
xx = np.arange(0, n ) * Delta_x 
X = (np.arange(1, n + 1) - 0.5) * Delta_x
lam = (1 - xx ** 2)/(Delta_x ** 2)

v = np.ones(n)
#lam = (v)/(Delta_x ** 2)

# Representation of sparse matrix and right-hand side
main  = np.zeros(n)
lower = np.zeros(n)
upper = np.zeros(n)


# Precompute sparse matrix
main[:] = -np.roll(lam, 1) - lam
lower[:] = lam  
upper[:] = lam  


A = sparse.diags(
    diagonals=[main, lower, upper],
    offsets=[0, -1, 1], shape=(n, n),
    format='csr')

# Insert boundary conditions
A[0, 0] = -lam[0]
A[0, 1] = lam[0]
A[-1, -2] = lam[-2]
A[-1, -1] = -lam[-2]
A = A.toarray()


US = lambda x: (1+ np.tanh(x))/2
LP = eval_legendre(2, X)


def equation(T, t, forcing): 
    """
    ODE to solve.
    """
    
    a = b0*US(Tc-T) + (a0 + a2 * LP) * US(T-Tc)
    S = 1 + s2*LP
    dTdt = (-A_*v - B*T + D*np.dot(A,T) + Q*S*a + forcing(t) * v)/c
    
    return dTdt


def forcing(t):
    return -0.5 + 0.03*t
# initial condition
T0 = 30 - ((np.arange(1, n + 1)-0.5)*Delta_x) * 50

# time vector
tinit = 200
t = np.linspace(0, tinit, 100)

sol = odeint(equation, T0, t, args= (forcing,))

fig, ax = plt.subplots()
def animate_plot(): 
    x = np.arange(500)
    line, = ax.plot(x, sol[0,:])
    return line,
    
def animate(i):
    if i < 100:
        line.set_ydata(sol[i, :])  # update the data.
    return line,

line, = animate_plot()

ani = animation.FuncAnimation(
    fig, animate, interval=20, blit=True, save_count=10, repeat = True)
plt.show()

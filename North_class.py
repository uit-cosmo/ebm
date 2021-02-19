import numpy as np
from scipy import sparse
from scipy.special import eval_legendre
from scipy.integrate import odeint
from matplotlib import pyplot as plt
import matplotlib.animation as animation

class North_1D(): 
    
    def  __init__(self, A, B, D, s1, s2, s22, Tc, b0, a0, a2, Q, c, n, T_0): 
        self.A = A 
        self.B = B
        self.D = D
        self.s1 = s1
        self.s2 = s2
        self.s22 = s22
        self.Tc = Tc
        self.b0 = b0
        self.a0 = a0
        self.a2 = a2
        self.Q = Q
        self.c = c
        self.n = n 
        self.T_0 = T_0
        self.Delta_theta = 1/self.n 
        self.X =  (np.arange(1, self.n + 1) - 0.5) * self.Delta_theta
        self.theta =  np.arange(0, self.n) * self.Delta_theta 
        self.r_E = 6371
        self.lam = (1 - self.theta ** 2) / (self.Delta_theta ** 2)
        self.v = np.ones(self.n)
        self.laplacian =self.__construct_laplacian().toarray()
        self.__set_boundary_conditions()
        self.sol = []
        self.fig, self.ax = plt.subplots()
        self.line = []
        self.t = t

    def __construct_laplacian(self): 
        """ 
        Laplacian Construction.
        """
        
        # Representation of sparse matrix and right-hand side
        main  = np.zeros(self.n)
        lower = np.zeros(self.n)
        upper = np.zeros(self.n)

        # Precompute sparse matrix
        main[:] = -np.roll(self.lam, 1) - self.lam
        lower[:] = self.lam  
        upper[:] = self.lam  

        laplacian = sparse.diags(
            diagonals=[main, lower, upper],
            offsets=[0, -1, 1], shape=(self.n, self.n),
            format='csr')
        
        return laplacian
    
    def __incident_solar_radiation(self, t): 
        """ 
        Calculates the seasonal incident solar radiation at latitude X, taken from North & Coakley 1979.
        
        Parameters:
        
        t: Time vector in years.    
        """
        
        return 1 + self.s1 * np.cos(2 * np.pi * t) * eval_legendre(1, self.X) + \
            (self.s2 + self.s22 * np.cos(4 * np.pi * t)) * eval_legendre(2, self.X)
            
    US = lambda self, x: (1+ np.tanh(x))/2
        
    def __set_boundary_conditions(self): 
        """
        Inserts the boundary conditions into the Laplacian matrix. 
        """
        
        self.laplacian[0, 0] = -self.lam[0]
        self.laplacian[0, 1] = self.lam[0]
        self.laplacian[-1, -2] = self.lam[-2]
        self.laplacian[-1, -1] = -self.lam[-2]
        

    def __forcing(self, t):
        return -0.5 + 0.03*t
    

    def __define_equation_seasonal(self, T, t): 
        """
        ODE to solve with seasonality.
        """
        
        a = self.b0*self.US(self.Tc-T) + (self.a0 + self.a2 * eval_legendre(2, self.X)) * self.US(T-self.Tc)
        dTdt = (-self.A*self.v - self.B*T + self.D*np.dot(self.laplacian, T) + self.Q*self.__incident_solar_radiation(t)*a + self.__forcing(t) * self.v)/self.c
        return dTdt
    
    def __define_equation(self, T, t): 
        """
        ODE to solve without seasonality.
        """

        a = self.b0*self.US(self.Tc-T) + (self.a0 + self.a2 * eval_legendre(2, self.X)) * self.US(T-self.Tc)
        S = 1 + self.s2 * eval_legendre(2, self.X)
        dTdt = (-self.A*self.v - self.B*T + self.D*np.dot(self.laplacian, T) + self.Q*S*a + self.__forcing(t) * self.v)/self.c
        return dTdt
    
    
    def solve_model(self, t, seasonality): 
        """
        Solve the 1D (latitudinal) EBM on the northern hemisphere using scipy.odeint. 
        Returns the solution array (see odeint documentation for the structure of the solution).
        
        Parameters:
        
        t: time array
        
        seasonality: True or False for EBM with or without seasonality.
        """
        
        self.t = t
        if seasonality == True: 
            self.sol = odeint(self.__define_equation_seasonal, self.T_0, t)
        elif seasonality == False: 
            self.sol = odeint(self.__define_equation, self.T_0, t)
        else: 
            raise Exception("No valid value.")
        return self.sol
    
    
    
    def __animate_plot(self): 
        x = np.arange(self.n)
        self.line, = self.ax.plot(x, self.sol[0, :])
        self.ax.set_ylim(-40,40)
        self.ax.axhline(-10, ls = "--")
        return self.line,
    
    
    def __animate(self, i):
        if i < self.t.size:
            self.line.set_ydata(self.sol[i, :])   # update the data.
        return self.line,

    
    def animate_solution(self): 
        """
        Animates the solution of the EBM from start to end.
        """
        self.line, = self.__animate_plot()
        ani = animation.FuncAnimation(
            self.fig, self.__animate, interval=20, blit=True, save_count=10, repeat = True)
        plt.show()
         
    
    def calculate_SIA(self): 
        """
        Calculates the sea ice area on the northern hemisphere. 
        
        Returns the area of the sea ice and the index of the latitude.
        """
        
        ind = (np.argmax(self.sol<-10, axis = 1))
        ind[ind == 0] = self.n - 1 
        A = 2 * np.pi * self.r_E**2 *(1 - self.theta[ind])
        return A, ind
        
        


tinit = 200
num = 3000
A = 211.2 - 18.
B = 1/0.32
D = 0.38
s1 = -0.796
s2 = -0.482
s22 = 0.147
Tc = -10
b0 = 0.38
a0 = 0.697
a2 = -0.0779
Q = 334.4
c = 13.2
tid = 50
tmax = 131
n = 500 

T0 = 30 - ((np.arange(1, n + 1) - 0.5) * 1/n) * 50
t = np.linspace(0, tinit, num)

m = North_1D(A = A, B = B, D = D, s1 = s1, s2 = s2, s22 = s22, Tc = Tc, b0 = b0, a0 = a0, a2 = a2, Q = Q, c = c, n = n, T_0 = T0)
m.solve_model(t = t, seasonality = True)
m.animate_solution()
A, ind = m.calculate_SIA()

plt.plot(A)
plt.show()

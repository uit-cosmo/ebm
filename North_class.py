import numpy as np
from scipy import sparse
from scipy.special import eval_legendre
from scipy.integrate import odeint
from matplotlib import pyplot as plt
import matplotlib.animation as animation
import os
import sys
class North_1D(): 
    
    def  __init__(self, A, B, D, s1, s2, s22, Tc, b0, a0, a2, Q, c, n, T_0, CO2_parameter): 
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
        self.CO2_parameter = CO2_parameter
        self.CO2_Base = 280 #1750 280 ppm CO2
        self.global_temp = np.zeros(t.size)
        plt.rcParams['font.size'] = 10
        plt.rcParams['axes.linewidth'] = 2
        plt.rcParams['axes.labelsize'] = 10
        plt.rcParams['axes.labelweight'] = 'bold'
        plt.rcParams['font.weight'] = 'bold'

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
    
    def __forcing_CO2(self, t): 
        """
        Formula for CO2 forcing. 
        
        Reference:
        Myhre, G., E.J. Highwood, K.P. Shine, and F. Stordal, 1998: 
        New estimates of radiative forcing due to well mixed greenhouse gases. 
        Geophysical Research Letters, 25, 2715-2718.
        """
        
        CO2ppm = 280
        return self.CO2_parameter*np.log((CO2ppm*1.02**(t)/self.CO2_Base))

    def __define_equation_seasonal(self, T, t, forcing, time_dependent_albedo): 
        """
        ODE to solve with seasonality.
        """
        if time_dependent_albedo == False:
            a = self.b0*self.US(self.Tc-self.T_0) + (self.a0 + self.a2 * eval_legendre(2, self.X)) * self.US(self.T_0-self.Tc)
        elif time_dependent_albedo == True: 
            a = self.b0*self.US(self.Tc-T) + (self.a0 + self.a2 * eval_legendre(2, self.X)) * self.US(T-self.Tc)
        else: 
            raise Exception("No valid value.")
        if forcing == False:
            dTdt = (-self.A*self.v - self.B*T + self.D*np.dot(self.laplacian, T) + self.Q*self.__incident_solar_radiation(t)*a + 0 * self.v)/self.c
        else: 
            dTdt = (-self.A*self.v - self.B*T + self.D*np.dot(self.laplacian, T) + self.Q*self.__incident_solar_radiation(t)*a + self.__forcing_CO2(t) * self.v)/self.c
        return dTdt
    
    def __define_equation(self, T, t, forcing, time_dependent_albedo): 
        """
        ODE to solve without seasonality.
        """
        if time_dependent_albedo == False:
            a = self.b0*self.US(self.Tc-self.T_0) + (self.a0 + self.a2 * eval_legendre(2, self.X)) * self.US(self.T_0-self.Tc)
        elif time_dependent_albedo == True: 
            a = self.b0*self.US(self.Tc-T) + (self.a0 + self.a2 * eval_legendre(2, self.X)) * self.US(T-self.Tc)
        else: 
            raise Exception("No valid value.")
        S = 1 + self.s2 * eval_legendre(2, self.X)
        if forcing == False:
            dTdt = (-self.A*self.v - self.B*T + self.D*np.dot(self.laplacian, T) + self.Q*S*a + 0 * self.v)/self.c
        else:
            dTdt = (-self.A*self.v - self.B*T + self.D*np.dot(self.laplacian, T) + self.Q*S*a + self.__forcing_CO2(t) * self.v)/self.c
        return dTdt
    
    
    def solve_model(self, t, seasonality, time_dependent_albedo): 
        """
        Solve the 1D (latitudinal) EBM on the northern hemisphere using scipy.odeint. 
        Returns the solution array (see odeint documentation for the structure of the solution).
        
        Parameters:
        
        t: time array
        
        seasonality: True or False for EBM with or without seasonality.
        
        time_dependent_albedo: True for time dependent albedo, False for constant albedo. Only for seasonal model. 
        """
        
        self.seasonality = seasonality 
        self.time_dependent_albedo = time_dependent_albedo
        self.t = t
        if seasonality == True: 
            # initial values calculated by non-seasonal model without forcing
            self.sol = odeint(self.__define_equation, self.T_0, t, args= (False, time_dependent_albedo))
            self.T_0 = self.sol[-1, :]
            self.sol = odeint(self.__define_equation_seasonal, self.T_0, t, args = (True, time_dependent_albedo))
        elif seasonality == False: 
            self.sol = odeint(self.__define_equation, self.T_0, t, args =  (True, time_dependent_albedo))
        else: 
            raise Exception("No valid value.")
        return self.sol
    
    
    
    def __animate_plot(self): 
        x = np.arange(self.n)
        self.line, = self.ax.plot(x, self.sol[0, :])
        self.ax.set_ylim(-40,40)
        self.ax.axhline(-10, ls = "--", color = "black")
        xticklabels = [str(x) for x in range(0, 90, 10)]
        self.ax.set_xlabel(r"Latitude $\theta$")
        self.ax.set_xticklabels(xticklabels)
        self.ax.set_xticks(np.arange(0, self.n + 1, int(self.n/90*10)))
        self.ax.grid(True, ls = '--')
        self.ax.set_ylabel(r"Temperature [$^\circ$C]")
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
        
        ind = (np.argmax(self.sol<-10, axis = 1)) #threshold temperature for ice
        ind[ind == 0] = self.n - 1 
        self.area = 2 * np.pi * self.r_E**2 *(1 - self.theta[ind])
        return self.area, ind
    
    def plot_SIA(self): 
        """
        Function to plot the sea ice area in the northern hemisphere against time. 
        """
        
        year_steps = 25
        self.fig, self.ax = plt.subplots()
        xticklabels = [str(x) for x in range(0, int(self.t[-1]), year_steps)]
        self.ax.plot(self.area, lw = 2)
        self.ax.set_xlabel("Modeltime [year]")
        self.ax.set_xticklabels(xticklabels)
        self.ax.set_xticks(np.arange(0, self.t.size, year_steps * self.t.size/int(self.t[-1])))
        self.ax.grid(True, ls = '--')
        self.ax.set_ylabel("SIA [km²]")
        plt.show()
        
    def __calculate_global_average_temp(self): 
        """
        Calculating the average global temperature for every year, weighted by the areas of the latitude bands. 
        """
        weight_area = 2 * np.pi * self.r_E**2 *(1 - self.theta)
        self.global_temp = np.average(self.sol, weights = weight_area, axis = 1)  
    
    def plot_SIA_month(self, month): 
        """
        Function to plot the sea ice area of a specific month in the northern hemisphere against time.
        
        Parameters:
        
        month = 0-11 (Jan-Dec) 
        """
        
        year_steps = 25
        month_list = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
        self.fig, self.ax = plt.subplots()
        xticklabels = [str(x) for x in range(0, int(self.t[-1]), year_steps)]

        self.ax.plot(self.area[month:-1:int(self.t.size/int(self.t[-1]))], lw = 2)
        self.ax.set_xlabel("Modeltime [year]")
        self.ax.set_xticklabels(xticklabels)
        self.ax.set_xticks(np.arange(0, self.t.size, year_steps * self.t.size/int(self.t[-1])))
        self.ax.grid(True, ls = '--')
        self.ax.set_ylabel("SIA in " + month_list[month] + " [km²]")
        plt.show()
        
    def save_solution(self):
        """
        Saving the solution array.
        """
        
        if self.time_dependent_albedo == True:
            if self.seasonality == True:  
                np.save('../ebm/sol_albd1_seasonality1.npy', self.sol)
            else: 
                np.save('../ebm/sol_td_alb1_seasonality0.npy', self.sol)
        else: 
            if self.seasonality == True:  
                np.save('../ebm/sol_albd0_seasonality1.npy', self.sol)
            else: 
                np.save('../ebm/sol_td_alb0_seasonality0.npy', self.sol)
                
    def save_SIA_against_T(self):
        """
        Saving the global weighted average temperature against the SIA (northern hemisphere).
        """
        self.__calculate_global_average_temp()
        T_against_SIA = np.vstack((self.global_temp, self.area))
        if self.time_dependent_albedo == True:
            if self.seasonality == True:  
                np.savetxt(f'../ebm/T_SIA_albd1_seasonality1_a{int(self.CO2_parameter*1e5)}.txt', T_against_SIA, fmt= "%e", delimiter=',')
            else: 
                np.savetxt(f'../ebm/T_SIA_albd1_seasonality0_a{int(self.CO2_parameter*1e5)}.txt', T_against_SIA, fmt= "%e", delimiter=',')
        else: 
            if self.seasonality == True:  
                np.savetxt(f'../ebm/T_SIA_albd0_seasonality1_a{int(self.CO2_parameter*1e5)}.txt', T_against_SIA, fmt= "%e", delimiter=',')
            else: 
                np.savetxt(f'../ebm/T_SIA_albd0_seasonality0_a{int(self.CO2_parameter*1e5)}.txt', T_against_SIA, fmt= "%e", delimiter=',')
    
        

# All the constants for the calculation. 
tinit = 100 #number of years
num = tinit*12 #12 steps per year
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
# Parameters chosen as for ACCCESS and CanESM model
a_parameter = 5.73471 #5.48945  
time_dependent_param  = False

# Initialize model
m = North_1D(A = A, B = B, D = D, s1 = s1, s2 = s2, s22 = s22, Tc = Tc, b0 = b0, a0 = a0, a2 = a2, Q = Q, c = c, n = n, T_0 = T0, CO2_parameter = a_parameter)

# Solve model
m.solve_model(t = t, seasonality = False, time_dependent_albedo = time_dependent_param)
# Calculating, saving and plotting the SIA 
A, ind = m.calculate_SIA()
m.save_SIA_against_T()
m.plot_SIA()
m.save_solution()
#m.animate_solution()
#m.plot_SIA_month(month = 8)
sys.exit()

# Saving the solution
m.save_solution()
# Animating the solution
m.animate_solution()


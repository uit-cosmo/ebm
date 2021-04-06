A Python script to solve the 1D-Energy Balance Equation on the Northern hemisphere, adapted from North et al. 

Options include:
- seasonal solar insolation
- time dependent albedo
- CO2-forcing
- noise

To initialize the model run: 

```m = North_1D(A = A, B = B, D = D, s1 = s1, s2 = s2, s22 = s22, Tc = Tc, b0 = b0, a0 = a0, a2 = a2, Q = Q, c = c, n = n, T_0 = T0, CO2_parameter = a_parameter, stdv_noise = stdv_noise)```

with your parameters. Example parameters are given in the file. 

Solve the model via: 

```
m.solve_model(t = t, seasonality = True, time_dependent_albedo = True, noise =  True)
``` 

where ```t``` is the time vector, ```seasonality``` decides if you want seasonal solar insolation (winter/summer cycle), ```time_dependent_albedo = True``` for albedo-feedback turned on, ```noise = True``` for Gaussian noise with ```stdv_noise```.

A Julia script to solve the 1D-Energy Balance Equation on the Northern hemisphere, adapted from North et al. 

Options include:
- seasonal solar insolation
- time dependent albedo
- CO2-forcing
- noise
- historical forcing
- ensemble runs

To solve the model run the file ```main.jl``` with your chosen options. Example parameters are given in the file. The number of runs for ensemble simulations is given in ```North.jl```. Running ```ensemble_analysis.jl``` calculates the autocorrelation at lag-1 and variance of the sea ice area. 

Note that the Python version is deprecated! 

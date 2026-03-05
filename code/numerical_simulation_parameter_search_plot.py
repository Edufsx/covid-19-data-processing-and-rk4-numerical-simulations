#%%
# Numerical operations
import numpy as np
# Data manipulation 
import pandas as pd
# Data visualization  
import matplotlib.pyplot as plt
# Progress Tracking
from tqdm import tqdm

# Define the ODE's for the simplified SIR epidemiological Model

# Differential Equation for the susceptible compartment
def S_dot(t, S, I, R, alpha, beta):
    return -alpha * I

# Differential Equation for the Infectious compartment
def I_dot(t, S, I, R, alpha, beta):
    return (alpha - beta) * I

# Differential Equation for the Recovered compartment
def R_dot(t, S, I, R, alpha, beta):
    return beta * I
  
# Fourth-order Runge-Kutta Method (RK4)
def rk4_step_3d(t, S_dot, I_dot, R_dot, S, I, R, alpha, beta, h):
    """ 
    Initial Value Problem (IVP):
    ----------------------------------------
    S' = S_dot; I' = I_dot; and R' = R_dot
    S_dot, I_dot and R_dot define the SIR system of ODEs
    t: current time step (before advancing to t + h)
    S, I, R: State variables at time t (before RK4 update)
    h: Integration step size 
    """

    half_h = h / 2
    k1_S = S_dot(t,          S,                   I,                   R,                   alpha, beta)
    k1_I = I_dot(t,          S,                   I,                   R,                   alpha, beta)
    k1_R = R_dot(t,          S,                   I,                   R,                   alpha, beta)
    
    k2_S = S_dot(t + half_h, S + (half_h * k1_S), I + (half_h * k1_I), R + (half_h * k1_R), alpha, beta)
    k2_I = I_dot(t + half_h, S + (half_h * k1_S), I + (half_h * k1_I), R + (half_h * k1_R), alpha, beta)  
    k2_R = R_dot(t + half_h, S + (half_h * k1_S), I + (half_h * k1_I), R + (half_h * k1_R), alpha, beta)     
    
    k3_S = S_dot(t + half_h, S + (half_h * k2_S), I + (half_h * k2_I), R + (half_h * k2_R), alpha, beta)
    k3_I = I_dot(t + half_h, S + (half_h * k2_S), I + (half_h * k2_I), R + (half_h * k2_R), alpha, beta)  
    k3_R = R_dot(t + half_h, S + (half_h * k2_S), I + (half_h * k2_I), R + (half_h * k2_R), alpha, beta)     

    k4_S = S_dot(t + h,      S + (h * k3_S),      I + (h * k3_I),      R + (h * k3_R),      alpha, beta)      
    k4_I = I_dot(t + h,      S + (h * k3_S),      I + (h * k3_I),      R + (h * k3_R),      alpha, beta)   
    k4_R = R_dot(t + h,      S + (h * k3_S),      I + (h * k3_I),      R + (h * k3_R),      alpha, beta)     

    S_next = S + (h * (k1_S + (2 * k2_S) + (2 * k3_S) + k4_S)) / 6    
    I_next = I + (h * (k1_I + (2 * k2_I) + (2 * k3_I) + k4_I)) / 6 
    R_next = R + (h * (k1_R + (2 * k2_R) + (2 * k3_R) + k4_R)) / 6     

    return S_next, I_next, R_next

# Fourth-order Runge-Kutta (RK4) System
def rk4_sir_system(S_dot, I_dot, R_dot, initial_conditions, t_start, t_end, alpha, beta, h):
    """
    Integrates a 3-equation ODE system using the RK4

    Parameters
    ----------------------------------------
    S_dot, I_dot and R_dot: functions that define the SIR system
    (x0, y0, z0): Initial condition for t = a
    [a, b]: Integration interval 
    h: Integration step

    Returns
    ----------------------------------------
    t_vector, S_vector, I_vector, R_vector : Time vector and numerical solutions.
    """

    # Step 0 : Initialize the S, I, R and time vector
    t_vector = [t_start]
    S0, I0, R0 = initial_conditions
    S_vector = [S0]
    I_vector = [I0]
    R_vector = [R0]    

    # Step 1 : Calculate the number of iterations (N) 
    N = int((t_end - t_start) / h)

    # Step 2 : Perform N times the rk4_step_3d
    for i in range(N):
        t_next = t_vector[i] + h
        S_next, I_next, R_next = rk4_step_3d(t_vector[i], S_dot, I_dot, R_dot, \
                                                 S_vector[i], I_vector[i], R_vector[i],\
                                                 alpha, beta, h)
        t_vector.append(t_next)
        S_vector.append(S_next)
        I_vector.append(I_next)
        R_vector.append(R_next)
        
    return t_vector, S_vector, I_vector, R_vector

def get_parameters(municipality, year):
    """
    Load the COVID-19 7-day moving average observed data and model Parameters
    for the selected municipality and year.
    """
     
    real_data = pd.read_csv(
        f"../data/filtered_data/{municipality}_casos_mm7d_{year}"
        )["cases_7-dma"].to_list() 
   
    params = municipality_year_params[municipality][year]
    t_end = params["t_end"]
    initial_conditions = params["initial_conditions"]
    alpha = params["alpha"]
    beta = params["beta"]

    return real_data, t_end, initial_conditions, alpha, beta

def extract_daily(I_vec, t_end, steps_per_day = 999):
    # Extract discrete daily time points (1 day= 0, 2 days = 999, 3 days = 2*999,...)
    daily_I = []
    for i in range(0, t_end+1):
        t_discreto = i*steps_per_day
        daily_I.append(I_vec[t_discreto])

    return daily_I

def search_best_parameters(municipality, year):
    """
    Search the best alpha and beta parameters for the SIR Model 
    for a selected municipality and year 
    by minimizing the cumulative squared error between
    the simulated infected cases and real COVID-19 data.
    ----------------------------------------
    Note: This grid search is computationally heavy (~80k simulations).
    tqdm is used to track progress and avoid uncertainty during long runs. 
    """
    
    # Get the COVID data and model parameters
    real_data, t_end, initial_conditions, alpha, beta = get_parameters(municipality, year)
    
    # Parameter grid for search of optimal alpha and beta
    alpha_vec = np.linspace(0.3, 0.7, 401)
    beta_vec = np.linspace(0.3, 0.7, 401)

    # Initialize matrix to store results: columns = [alpha, beta, error]
    simulation_matrix = np.zeros([len(alpha_vec) * len(beta_vec), 3])
    row = 0

    # Grid search over alpha and beta
    for alpha_candidate in tqdm(alpha_vec):
        for beta_candidate in beta_vec:
            # Infection rate should be greater than recovery rate: alpha > beta
            if alpha_candidate > beta_candidate:

                # Simulate SIR Model for the current alpha and beta parameters 
                t_vec, S_vec, I_vec, R_vec = rk4_sir_system(S_dot, I_dot, R_dot, \
                                                                initial_conditions, t_start, t_end, \
                                                                alpha_candidate, beta_candidate, h)
                
                # Extracted daily infected values
                daily_I = extract_daily(I_vec, t_end)
                     
                # Compute cumulative squared error between simulated and observed data
                error = np.sum((np.array(daily_I) - np.array(real_data))**2)

                # Store results in simulation matrix           
                simulation_matrix[row, :] = np.array([alpha_candidate, beta_candidate, error])
                row += 1    

    # Trim unused rows
    trimmed_matrix= simulation_matrix[0:row,:]
    
    # Find parameters with minimal error
    min_index=np.argmin(trimmed_matrix[:,2])
    optimal_parameters = simulation_matrix[min_index,:]

    print("Optimal [alpha, beta, error]:", optimal_parameters)
    
    return optimal_parameters, trimmed_matrix

def get_best_I0(I0_end_interval, municipality, year):
    """
    Search the best I0 (Initial number of infected people) for the SIR Model 
    for a selected municipality and year 
    (considering that the real data could be missing the I0) 
    by minimizing the cumulative squared error between
    the simulated infected cases and real COVID-19 data. 
    """

    # Get the COVID data and model parameters
    real_data, t_end, initial_conditions, alpha, beta = get_parameters(municipality, year)

    # Initialize matrix to store results of I0 and error
    error_vec=np.zeros([I0_end_interval,2])
    I0 = 0

    for k in tqdm(range(0,I0_end_interval)):

        # Change the I0 value
        I0 += 1
        initial_conditions[1] = I0

        # Simulate SIR Model for the current I0  
        t_vector, S_vector, I_vector, R_vector = rk4_sir_system(S_dot, I_dot, R_dot, \
                                                    initial_conditions, t_start, t_end, \
                                                    alpha, beta, h)
        
       # Extracted daily infected values
        daily_I = extract_daily(I_vector, t_end)

        # Compute cumulative squared error between simulated and observed data# Compute
        error = np.sum((np.array(daily_I) - np.array(real_data))**2)
        
        # Saves the current I0 and cumulative squared error on the array result
        result = np.array([I0, error])               
        error_vec[k, :] = result    

    # Find parameters with minimal error
    min_index=np.argmin(error_vec[:,1])
    optimal_I0 = error_vec[min_index,:]
    print(error_vec[min_index,:])

    return optimal_I0, error_vec

def plot_cases(municipality, year):

    # Get the COVID data and model parameters
    real_data, t_end, initial_conditions, alpha, beta  = get_parameters(municipality, year)

    # Simulate SIR Model  
    t_vector, S_vector, I_vector, R_vector = rk4_sir_system(S_dot, I_dot, R_dot, \
                                               initial_conditions, t_start, t_end, \
                                               alpha, beta, h)
    
    # Plot graphs comparing the simulated and actual COVID-19 cumulative cases over time 
    plt.figure()
    plt.plot(t_vector, I_vector)
    plt.plot(real_data)
    plt.xlabel(f'Dias - {year}')
    plt.ylabel(f'Número de casos Confirmados em {municipality}')
    plt.title(f"Casos de Covid em {municipality} (Oficial e Simulado) ao longo de {year}")
    plt.legend(["Simulado", "Oficial"])
    plt.show()

# Plot graphs using "plot_cases()" for São Paulo and Bauru (2020-2022)  
def plot_6_graphs(municipalities, years):
    for municipality in municipalities:
        for year in years:
            plot_cases(municipality, year)

# Dictionary containing the parameters obtained from research 
# and simulations for each municipality and year
municipality_year_params = {
    "Bauru" : {
        2020 : {"alpha" : 0.576, 
                "beta" : 0.541, 
                "t_end" : 272,
                "initial_conditions" :[364223, 2, 0]
                },
        2021 : {"alpha" :  0.541, 
                "beta" : 0.516, 
                "t_end" : 365,
                "initial_conditions" :[345928, 8, 18289]
                },
        2022 : {"alpha" : 0.360, 
                "beta" : 0.330,
                "t_end": 365,
                "initial_conditions" : [310282, 1, 53942]
                }
},   
    "São Paulo" : {
        2020 : {"alpha" : 0.584, 
                "beta" : 0.541, 
                "t_end" : 310,
                "initial_conditions" : [11869659, 1, 0]
                },
        2021 : {"alpha" :  0.541, 
                "beta" : 0.522, 
                "t_end" : 365,
                "initial_conditions" : [11466825, 1117, 401718]
                },
        2022 : {"alpha" : 0.541, 
                "beta" : 0.518,
                "t_end": 365,
                "initial_conditions" : [10891674, 68, 977918]
                }}
}

# Numeric Integration parameters that remain constant across municipalities or years once defined
h = 10**(-3)
t_start = 0

# Define municipalities to filter COVID-19 data for São Paulo and Bauru (2020 - 2022)
municipalities = ["São Paulo", "Bauru"]
years = [2020, 2021, 2022]

# Plot graphs for São Paulo and Bauru (2020-2022)
plot_6_graphs(municipalities, years) 

""" 
# This function takes a long time to run:
search_best_parameters(municipality=municipalities[0], year=years[0])
"""


""" # This function just makes sense once we have the best alpha and beta parameters
get_best_I0(I0_end_interval=13, municipality=municipalities[0], year=years[0])
 """
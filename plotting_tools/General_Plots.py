import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd

def pendulum_ode(t, y):
    """
    Pendulum ODE system
    y[0] = θ
    y[1] = dθ/dt
    """
    g = 9.81
    L = 1.0
    
    dydt = np.zeros(2)
    dydt[0] = y[1]
    dydt[1] = -(g/L) * np.sin(y[0])
    return dydt

# Initial conditions
y0 = [np.pi/4, 0]  # Initial angle = π/4, initial angular velocity = 0
t_span = [0, 5]   # Time span from 0 to 10 seconds

# Solve using RK45 method
sol = solve_ivp(pendulum_ode, t_span, y0, method='RK45', 
                rtol=1e-10, atol=1e-10,  # Tight tolerances for accuracy
                dense_output=True)        # Enable dense output for smooth plotting

# Plot solution
plt.figure(figsize=(12, 8))
t = np.linspace(t_span[0], t_span[1], 1000)
y = sol.sol(t)  # Dense output evaluation
plt.plot(t, y[0], label='RK45')

plt.xlabel('Time (s)')
plt.ylabel('Angle (θ)')
plt.title('Pendulum Solutions with Different ODE Solvers')
plt.legend()
plt.grid(True)
plt.show()

# Read the CSV file. Make sure ODE_Pendulum_Solution.csv is in the same directory or adjust the path accordingly.
data = pd.read_csv('plotting_tools/ODE_Pendulum_Solution.csv')

# Extract the first column (time) and second column (angle theta)
time  = data.iloc[:, 0]
theta = data.iloc[:, 1]

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(time, theta, label=r'$\theta(t)$')
plt.xlabel('Time')
plt.ylabel('Angle (θ)')
plt.title('Pendulum Angle vs. Time')
plt.legend()
plt.grid(True)
plt.show()
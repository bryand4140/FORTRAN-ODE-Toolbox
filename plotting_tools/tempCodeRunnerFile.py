import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file. Make sure ODE_Pendulum_Solution.csv is in the same directory or adjust the path accordingly.
data = pd.read_csv('ODE_Pendulum_Solution.csv')

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
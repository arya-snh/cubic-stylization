import matplotlib.pyplot as plt
import pandas as pd

# Replace with your CSV file path
file_path = 'objective_values.csv'

# Read the data from the CSV file
# Assuming no header in the file and columns are iteration, objective value
data = pd.read_csv(file_path, header=None, names=["Iteration", "Objective"])

plt.style.use('seaborn-v0_8-pastel')
# Number of rows for each parameter
rows_per_param = 5 # Adjust this based on your data

# Assigning parameter values manually
data['Parameter'] = 0.9
data.loc[rows_per_param:rows_per_param*2-1, 'Parameter'] = 0.7
data.loc[rows_per_param*2:, 'Parameter'] = 0.5

# Plotting
plt.figure(figsize=(10, 6))

# Plot each parameter's data
for param in data['Parameter'].unique():
    subset = data[data['Parameter'] == param]
    plt.plot(subset["Iteration"], subset["Objective"], marker='o', label=f'Parameter {param}')

plt.title("Objective Value vs Iteration for Different Parameters")
plt.xlabel("Iteration")
plt.ylabel("Objective Value")
plt.legend()
plt.grid(True)
plt.show()

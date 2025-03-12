import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import os

# Load data from the CSV file
csv_path = "./TimeStepTest/timestep_results.csv"

# Check if file exists
if not os.path.exists(csv_path):
    raise FileNotFoundError(f"CSV file not found at {csv_path}")

# Read the CSV into a DataFrame
df = pd.read_csv(csv_path)

# Clean and prepare data
df = df.dropna()

# Extract data for logarithmic linear regression
dt_col = "TimeStep" if "TimeStep" in df.columns else "dt"
final_energy_col = "FinalEnergy" if "FinalEnergy" in df.columns else "final_energy"

# Original time step values, log-transformed final energy
dt_values = df[dt_col].values
log_final_energies = np.log10(np.abs(df[final_energy_col].values))

# Perform linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(dt_values, log_final_energies)

# Print regression results
print("Partial Logarithmic Linear Regression Analysis Results")
print("---------------------------------------------")
print(f"Slope: {slope:.6f}")
print(f"Intercept: {intercept:.6f}")
print(f"R-squared: {r_value**2:.6f}")
print(f"P-value: {p_value:.6f}")
print(f"Standard Error: {std_err:.6f}")

# Create plot
plt.figure(figsize=(12, 8))

# Scatter plot of actual data points
plt.scatter(dt_values, log_final_energies, color='blue', label='Log(Final Energy) Data Points')

# Linear regression line
line = slope * dt_values + intercept
plt.plot(dt_values, line, color='red', 
         label=f'Regression Line (log(y) = {slope:.4f}*x + {intercept:.4f})')

# Add labels and title
plt.xlabel('Time Step Size')
plt.ylabel('Log(Final Energy Absolute Value)')
plt.title(f'Linear Regression with Logarithmic Energy: Slope = {slope:.2f}')
plt.grid(True)
plt.legend()

# Add text annotation with regression details
plt.text(0.05, 0.95, 
         f'Slope: {slope:.4f}\n'
         f'Intercept: {intercept:.4f}\n'
         f'R²: {r_value**2:.4f}',
         transform=plt.gca().transAxes,
         verticalalignment='top',
         bbox=dict(facecolor='white', alpha=0.8))

plt.tight_layout()

# Save the plot
output_dir = os.path.dirname(csv_path)
plt.savefig(os.path.join(output_dir, "partial_log_regression_final_energy.png"), dpi=300)

# Create a second plot to show original energy scale
plt.figure(figsize=(12, 8))
plt.scatter(dt_values, np.abs(df[final_energy_col].values), color='blue', label='Final Energy Data Points')
y_pred = 10**(intercept) * dt_values**slope
plt.plot(dt_values, y_pred, color='red', label=f'Predicted Curve')
plt.xlabel('Time Step Size')
plt.ylabel('Final Energy (Absolute Value)')
plt.title(f'Predicted Energy Trend: Power Law Slope = {slope:.2f}')
plt.yscale('log')  # Log scale for energy
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "predicted_energy_trend.png"), dpi=300)

print(f"Partial log regression analysis complete. Plots saved to {output_dir}")
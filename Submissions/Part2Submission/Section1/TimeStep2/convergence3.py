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

# Extract data for linear regression
dt_col = "TimeStep" if "TimeStep" in df.columns else "dt"
final_energy_col = "FinalEnergy" if "FinalEnergy" in df.columns else "final_energy"

dt_values = df[dt_col].values
final_energies = df[final_energy_col].values

# Filter out any invalid values
valid_indices = ~np.isnan(final_energies) & ~np.isinf(final_energies) & (dt_values > 0)
dt_values = dt_values[valid_indices]
final_energies = final_energies[valid_indices]

# Perform linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(dt_values, final_energies)

# Print regression results
print("Linear Regression Analysis Results")
print("--------------------------------")
print(f"Slope: {slope:.6f}")
print(f"Intercept: {intercept:.6f}")
print(f"R-squared: {r_value**2:.6f}")
print(f"P-value: {p_value:.6f}")
print(f"Standard Error: {std_err:.6f}")

# Create linear regression plot
plt.figure(figsize=(12, 8))

# Scatter plot of actual data points
plt.scatter(dt_values, final_energies, color='blue', label='Data Points')

# Linear regression line
line = slope * dt_values + intercept
plt.plot(dt_values, line, color='red', label=f'Regression Line (y = {slope:.4f}x + {intercept:.4f})')

# Add labels and title
plt.xlabel('Time Step Size')
plt.ylabel('Final Energy')
plt.title('Linear Regression: Final Energy vs Time Step')
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
plt.savefig(os.path.join(output_dir, "linear_regression_final_energy.png"), dpi=300)

print(f"Linear regression analysis complete. Plot saved to {output_dir}")
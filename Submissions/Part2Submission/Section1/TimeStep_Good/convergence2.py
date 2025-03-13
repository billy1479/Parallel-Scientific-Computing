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

# Print some information about the loaded data
print(f"Loaded {len(df)} data points from {csv_path}")
print(f"Columns: {df.columns.tolist()}")
print(f"Sample data:\n{df.head()}")

# Clean and prepare data
# Remove rows with NaN values
df = df.dropna()
print(f"After removing NaN values: {len(df)} data points remain")

# Extract data for plotting - column names may be different in the CSV
dt_col = "TimeStep" if "TimeStep" in df.columns else "dt"
final_energy_col = "FinalEnergy" if "FinalEnergy" in df.columns else "final_energy"

dt_values = df[dt_col].values
final_energies = df[final_energy_col].values

# Filter out any invalid values
valid_indices = ~np.isnan(final_energies) & ~np.isinf(final_energies) & (dt_values > 0)
dt_values = dt_values[valid_indices]
final_energies = final_energies[valid_indices]

print(f"Valid data points for convergence analysis: {len(dt_values)}")

# Perform log-log linear regression to find convergence order
log_dt = np.log10(dt_values)
log_final_energy = np.log10(np.abs(final_energies))  # Take absolute value for log

# Use scipy.stats for linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(log_dt, log_final_energy)

print(f"Convergence Order Analysis Results (using Final Energy)")
print(f"--------------------------------")
print(f"Slope (Order): {slope:.6f}")
print(f"R-squared: {r_value**2:.6f}")
print(f"Standard Error: {std_err:.6f}")
print(f"P-value: {p_value:.6f}")
print(f"Intercept: {intercept:.6f}")

# Create log-log plot
plt.figure(figsize=(12, 8))

# Plot data points
plt.loglog(dt_values, np.abs(final_energies), 'bo', label='Final Energy (absolute value)')

# Plot regression line
x_regression = np.logspace(np.log10(dt_values.min()), np.log10(dt_values.max()), 100)
y_regression = 10**(slope * np.log10(x_regression) + intercept)
plt.loglog(x_regression, y_regression, 'r-', 
         label=f'Regression Line (Order = {slope:.2f})')

# Reference lines for 1st and 2nd order slopes
x_ref = np.logspace(np.log10(dt_values.min()), np.log10(dt_values.max()), 100)

# Find suitable scaling factors for reference lines by matching at the midpoint
mid_idx = len(dt_values) // 2
mid_dt = dt_values[mid_idx]
mid_energy = np.abs(final_energies[mid_idx])

# First order reference (slope = 1)
scale_first = mid_energy / mid_dt
y_first = x_ref * scale_first
plt.loglog(x_ref, y_first, 'g--', label='1st Order Reference (slope=1)')

# Second order reference (slope = 2)
scale_second = mid_energy / (mid_dt**2)
y_second = x_ref**2 * scale_second
plt.loglog(x_ref, y_second, 'm--', label='2nd Order Reference (slope=2)')

# Add labels and grid
plt.xlabel('Time Step Size (dt)')
plt.ylabel('Final Energy (absolute value)')
plt.grid(True, which="both", ls="-")
plt.title(f'Final Energy vs Time Step: Order = {slope:.2f} (R² = {r_value**2:.3f})')

# Add text annotation with detailed info
plt.text(0.05, 0.05, 
       f'Order: {slope:.4f}\n'
       f'R²: {r_value**2:.4f}\n'
       f'Standard Error: {std_err:.4f}',
       transform=plt.gca().transAxes,
       bbox=dict(facecolor='white', alpha=0.8))

plt.legend()
plt.tight_layout()

# Save the output directory
output_dir = os.path.dirname(csv_path)
plt.savefig(os.path.join(output_dir, "final_energy_convergence.png"), dpi=300)

# Also plot on linear scale
plt.figure(figsize=(12, 6))
plt.plot(dt_values, np.abs(final_energies), 'bo-')
plt.xlabel('Time Step Size (dt)')
plt.ylabel('Final Energy (absolute value)')
plt.grid(True)
plt.title('Final Energy vs Time Step Size (Linear Scale)')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "final_energy_linear.png"), dpi=300)

print(f"Analysis complete. Results saved to {output_dir}")
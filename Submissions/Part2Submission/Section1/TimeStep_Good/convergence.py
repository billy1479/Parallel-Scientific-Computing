import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import os

# Load data from the CSV file
csv_path = "./TimeStepTest_EnergyError/timestep_results.csv"

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
init_energy_col = "InitialEnergy" if "InitialEnergy" in df.columns else "initial_energy"
final_energy_col = "FinalEnergy" if "FinalEnergy" in df.columns else "final_energy"
rel_error_col = "RelativeError" if "RelativeError" in df.columns else "relative_error"

dt_values = df[dt_col].values
initial_energies = df[init_energy_col].values
final_energies = df[final_energy_col].values

# Check if we have RelativeError column, otherwise calculate it
if rel_error_col in df.columns:
    rel_errors = df[rel_error_col].values
else:
    # Calculate relative error
    rel_errors = np.abs((final_energies - initial_energies) / initial_energies)
    
# Filter out any invalid values
valid_indices = ~np.isnan(rel_errors) & ~np.isinf(rel_errors) & (rel_errors > 0) & (dt_values > 0)
dt_values = dt_values[valid_indices]
rel_errors = rel_errors[valid_indices]

print(f"Valid data points for convergence analysis: {len(dt_values)}")

# Perform log-log linear regression to find convergence order
log_dt = np.log10(dt_values)
log_errors = np.log10(rel_errors)

# Use scipy.stats for linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(log_dt, log_errors)

print(f"Convergence Order Analysis Results")
print(f"--------------------------------")
print(f"Convergence Order: {slope:.6f}")
print(f"R-squared: {r_value**2:.6f}")
print(f"Standard Error: {std_err:.6f}")
print(f"P-value: {p_value:.6f}")
print(f"Intercept: {intercept:.6f}")
print()

print(f"Interpretation:")
print(f"For an N-body simulation with Euler integration, the expected order is 1.0.")
print(f"For a second-order method (e.g., leapfrog), the expected order is 2.0.")
print(f"For a fourth-order method (e.g., RK4), the expected order is 4.0.")
print()

if abs(slope - 1.0) < 0.2:
    print(f"This simulation appears to use a FIRST-ORDER method (like Euler integration).")
elif abs(slope - 2.0) < 0.2:
    print(f"This simulation appears to use a SECOND-ORDER method (like leapfrog/Verlet).")
elif abs(slope - 4.0) < 0.2:
    print(f"This simulation appears to use a FOURTH-ORDER method (like RK4).")
else:
    print(f"The convergence order of {slope:.2f} doesn't match standard methods.")
    print(f"This might indicate mixed-order terms or non-standard integration.")

# Create log-log plot
plt.figure(figsize=(12, 8))

# Plot data points
plt.loglog(dt_values, rel_errors, 'bo', label='Simulation Data')

# Plot regression line
x_regression = np.logspace(np.log10(dt_values.min()), np.log10(dt_values.max()), 100)
y_regression = 10**(slope * np.log10(x_regression) + intercept)
plt.loglog(x_regression, y_regression, 'r-', 
         label=f'Regression Line (Order = {slope:.2f})')

# Reference lines for 1st and 2nd order convergence
x_ref = np.logspace(np.log10(dt_values.min()), np.log10(dt_values.max()), 100)

# Find suitable scaling factors for reference lines by matching at the midpoint
mid_idx = len(dt_values) // 2
mid_dt = dt_values[mid_idx]
mid_err = rel_errors[mid_idx]

# First order reference (slope = 1)
scale_first = mid_err / mid_dt
y_first = x_ref * scale_first
plt.loglog(x_ref, y_first, 'g--', label='1st Order Reference (slope=1)')

# Second order reference (slope = 2)
scale_second = mid_err / (mid_dt**2)
y_second = x_ref**2 * scale_second
plt.loglog(x_ref, y_second, 'm--', label='2nd Order Reference (slope=2)')

# Add labels and grid
plt.xlabel('Time Step Size (dt)')
plt.ylabel('Relative Energy Error')
plt.grid(True, which="both", ls="-")
plt.title(f'Convergence Order Analysis: Order = {slope:.2f} (R² = {r_value**2:.3f})')

# Add text annotation with detailed info
plt.text(0.05, 0.05, 
       f'Convergence Order: {slope:.4f}\n'
       f'R²: {r_value**2:.4f}\n'
       f'Standard Error: {std_err:.4f}',
       transform=plt.gca().transAxes,
       bbox=dict(facecolor='white', alpha=0.8))

plt.legend()
plt.tight_layout()
plt.savefig("convergence_order.png", dpi=300)
# Save figures to the same directory as the input CSV
output_dir = os.path.dirname(csv_path)

# Also plot the absolute errors vs dt (not in log scale)
plt.figure(figsize=(12, 6))
plt.plot(dt_values, rel_errors, 'bo-')
plt.xlabel('Time Step Size (dt)')
plt.ylabel('Relative Energy Error')
plt.grid(True)
plt.title('Relative Error vs Time Step Size (Linear Scale)')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "error_vs_dt_linear.png"), dpi=300)

# Save the convergence order plot to the same directory
plt.figure(figsize=(12, 8))
plt.title(f'Convergence Order = {slope:.2f}')
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "convergence_order.png"), dpi=300)

# Save the convergence analysis results to a text file
results_path = os.path.join(output_dir, "convergence_results.txt")
with open(results_path, 'w') as f:
    f.write(f"Convergence Order Analysis Results\n")
    f.write(f"--------------------------------\n")
    f.write(f"Convergence Order: {slope:.6f}\n")
    f.write(f"R-squared: {r_value**2:.6f}\n")
    f.write(f"Standard Error: {std_err:.6f}\n")
    f.write(f"P-value: {p_value:.6f}\n")
    f.write(f"Intercept: {intercept:.6f}\n\n")
    
    f.write(f"Interpretation:\n")
    f.write(f"For an N-body simulation with Euler integration, the expected order is 1.0.\n")
    f.write(f"For a second-order method (e.g., leapfrog), the expected order is 2.0.\n")
    f.write(f"For a fourth-order method (e.g., RK4), the expected order is 4.0.\n\n")
    
    if abs(slope - 1.0) < 0.2:
        f.write(f"This simulation appears to use a FIRST-ORDER method (like Euler integration).\n")
    elif abs(slope - 2.0) < 0.2:
        f.write(f"This simulation appears to use a SECOND-ORDER method (like leapfrog/Verlet).\n")
    elif abs(slope - 4.0) < 0.2:
        f.write(f"This simulation appears to use a FOURTH-ORDER method (like RK4).\n")
    else:
        f.write(f"The convergence order of {slope:.2f} doesn't match standard methods.\n")
        f.write(f"This might indicate mixed-order terms or non-standard integration.\n")

print(f"Analysis complete. Results saved to {output_dir}")
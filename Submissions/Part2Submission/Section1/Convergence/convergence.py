import subprocess
import os
import shutil
import argparse
from pathlib import Path
import logging
import csv
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

class TimeStepEnergyTester:
    """
    A class to test N-Body simulations with different time steps and monitor energy conservation.
    """
    
    def __init__(self, log_level=logging.INFO):
        """Initialize the tester with logging configuration."""
        # Configure logging
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger('TimeStepTester')
    
    def run_test_series(self, 
                        python_script="create_initial_conditions.py",
                        N=10,
                        min_mass=1000.0,
                        max_mass=1000000.0,
                        start_dt=0.001,
                        end_dt=1.5,
                        dt_increment=0.0001,
                        final_time=60.0,
                        executable_name="./step-0",
                        output_dir="TimeStepTest",
                        result_file="timestep_results.csv"):
        """
        Run a series of simulations with increasing time steps and record energy values.
        
        Args:
            python_script: Path to the script that generates initial conditions
            N: Number of bodies in the simulation
            min_mass: Minimum mass for bodies
            max_mass: Maximum mass for bodies
            start_dt: Starting time step size
            end_dt: Ending time step size
            dt_increment: Increment for time step between tests
            final_time: Final simulation time
            executable_name: Name of the executable to use
            output_dir: Directory to store all outputs
            result_file: CSV file to store the results
        """
        # Create output directory
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Create results file with header
        results_path = output_path / result_file
        with open(results_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["TimeStep", "InitialEnergy", "FinalEnergy", "RelativeError"])
        
        # Calculate time steps to test
        time_steps = np.arange(start_dt, end_dt + dt_increment/2, dt_increment)
        total_tests = len(time_steps)
        
        self.logger.info(f"Starting test series with {total_tests} different time steps")
        self.logger.info(f"Time steps from {start_dt} to {end_dt} with increment {dt_increment}")
        
        # Store data for plotting
        dt_values = []
        initial_energy_values = []
        final_energy_values = []
        relative_errors = []
        
        # Run tests for each time step
        for i, dt in enumerate(time_steps):
            self.logger.info(f"Test {i+1}/{total_tests}: Time step = {dt:.6f}")
            
            # Run simulation with this time step
            initial_energy, final_energy = self._run_single_test(
                python_script=python_script,
                N=N,
                min_mass=min_mass,
                max_mass=max_mass,
                dt=dt,
                final_time=final_time,
                executable_name=executable_name,
                output_dir=output_dir
            )
            
            # Calculate relative error
            if initial_energy is not None and final_energy is not None and initial_energy != 0:
                rel_error = abs((final_energy - initial_energy) / initial_energy)
            else:
                rel_error = float('nan')
            
            # Store for plotting
            dt_values.append(dt)
            initial_energy_values.append(initial_energy)
            final_energy_values.append(final_energy)
            relative_errors.append(rel_error)
            
            # Write results to CSV
            with open(results_path, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([dt, initial_energy, final_energy, rel_error])
            
            # Print simple output as requested
            print(f"Time Step: {dt:.6f}, Initial Energy: {initial_energy}, Final Energy: {final_energy}")
        
        self.logger.info(f"Test series completed. Results saved to {results_path}")
        
        # Create plot of the results
        self._create_energy_plot(dt_values, initial_energy_values, final_energy_values, 
                               relative_errors, output_dir)
        
        # Create convergence order plot
        self._create_convergence_plot(dt_values, relative_errors, output_dir)
        
        # Analyze and recommend stable time step
        self._analyze_stability(dt_values, initial_energy_values, final_energy_values, 
                              relative_errors)
        
        # Analyze convergence order
        self._analyze_convergence_order(dt_values, relative_errors, output_dir)
    
    def _run_single_test(self, 
                         python_script,
                         N,
                         min_mass,
                         max_mass,
                         dt,
                         final_time,
                         executable_name,
                         output_dir):
        """
        Run a single simulation test and extract energy values.
        
        Returns:
            tuple: (initial_energy, final_energy)
        """
        try:
            # Create a unique subdirectory for this test
            test_dir = f"{output_dir}/dt_{dt:.6f}"
            os.makedirs(test_dir, exist_ok=True)
            
            # Generate the shell script with create_initial_conditions.py
            python_cmd = [
                "python3", python_script,
                "--final-time", str(final_time),
                "--snapshots", "0",
                "--executable-name", executable_name,
                "--min-mass", str(min_mass),
                "--max-mass", str(max_mass),
                "--dt", str(dt),
                "--N", str(N)
            ]
            
            subprocess.run(
                python_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            
            # Get the shell script name
            shell_script = f"{os.path.basename(executable_name)}.sh"
            
            # Ensure it's executable
            os.chmod(shell_script, 0o755)
            
            # Output file path
            output_file_path = os.path.join(test_dir, "simulation_output.txt")
            
            # Run the simulation with explicit bash call
            with open(output_file_path, "w") as output_file:
                subprocess.run(
                    f"bash ./{shell_script}",
                    shell=True,
                    stdout=output_file,
                    stderr=subprocess.STDOUT,
                    env=os.environ.copy()
                )
            
            # Extract energy values from output
            initial_energy, final_energy = self._extract_energy_values(output_file_path)
            return initial_energy, final_energy
            
        except Exception as e:
            self.logger.error(f"Error in test with dt={dt}: {str(e)}")
            return None, None
    
    def _extract_energy_values(self, output_file):
        """
        Extract initial and final energy values from simulation output.
        
        Args:
            output_file: Path to the simulation output file
            
        Returns:
            tuple: (initial_energy, final_energy)
        """
        initial_energy = None
        final_energy = None
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
                
            # Debug help: If we didn't find energy values, log the file content
            if "Initial energy:" not in content and "Total energy:" not in content:
                self.logger.debug(f"Energy values not found in output file. First 500 chars: {content[:500]}")
                
            # Try to extract energy values using various patterns
            with open(output_file, 'r') as f:
                for line in f:
                    if "Initial energy:" in line:
                        parts = line.strip().split(":")
                        if len(parts) > 1:
                            try:
                                initial_energy = float(parts[1].strip())
                            except ValueError:
                                self.logger.warning(f"Could not parse initial energy value: {parts[1].strip()}")
                    
                    # Try to match either "Final energy:" or "Total energy:" for the final value
                    elif "Final energy:" in line or "Total energy:" in line:
                        parts = line.strip().split(":")
                        if len(parts) > 1:
                            try:
                                final_energy = float(parts[1].strip())
                            except ValueError:
                                self.logger.warning(f"Could not parse final energy value: {parts[1].strip()}")
        except Exception as e:
            self.logger.warning(f"Error extracting energy values: {str(e)}")
        
        return initial_energy, final_energy
        
    def _create_energy_plot(self, dt_values, initial_energies, final_energies, rel_errors, output_dir):
        """
        Create a plot of energy values vs time step size.
        
        Args:
            dt_values: List of time step values
            initial_energies: List of initial energy values
            final_energies: List of final energy values
            rel_errors: List of relative errors
            output_dir: Directory to save the plot
        """
        try:
            # Convert to numpy arrays
            dt_values = np.array(dt_values)
            initial_energies = np.array(initial_energies)
            final_energies = np.array(final_energies)
            rel_errors = np.array(rel_errors)
            
            # Create the plot
            plt.figure(figsize=(12, 8))
            
            # Plot initial and final energy values
            plt.plot(dt_values, initial_energies, 'b-', label='Initial Energy')
            plt.plot(dt_values, final_energies, 'r-', label='Final Energy')
            
            # Use scientific notation for y-axis if values are very large
            if np.nanmax(np.abs(initial_energies)) > 1000:
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            
            # Plot relative error on secondary axis if we have valid values
            valid_indices = ~np.isnan(rel_errors)
            if np.any(valid_indices):
                ax2 = plt.gca().twinx()
                ax2.plot(dt_values[valid_indices], rel_errors[valid_indices], 'g--', label='Relative Error')
                ax2.set_ylabel('Relative Error')
                ax2.set_yscale('log')
                
                # Add reference lines for 0.1% and 1% error thresholds
                ax2.axhline(y=0.001, color='g', linestyle=':', alpha=0.5, label='0.1% Error')
                ax2.axhline(y=0.01, color='r', linestyle=':', alpha=0.5, label='1% Error')
                
                # Get handles and labels for second axis legend
                lines2, labels2 = ax2.get_legend_handles_labels()
            else:
                lines2, labels2 = [], []
            
            # Find where NaN values start (instability point)
            nan_indices = np.isnan(final_energies)
            if np.any(nan_indices) and not np.all(nan_indices):
                first_nan_idx = np.where(nan_indices)[0][0]
                if first_nan_idx > 0:
                    instability_dt = dt_values[first_nan_idx]
                    plt.axvline(x=instability_dt, color='r', linestyle='-', alpha=0.5)
                    plt.text(instability_dt, plt.ylim()[0], f'  Instability at dt={instability_dt:.4f}', 
                            verticalalignment='bottom', horizontalalignment='left')
            
            # Add labels and grid
            plt.xlabel('Time Step Size (dt)')
            plt.ylabel('Energy')
            plt.grid(True)
            plt.title('Energy vs Time Step Size')
            
            # Combine legends from both axes
            lines1, labels1 = plt.gca().get_legend_handles_labels()
            plt.legend(lines1 + lines2, labels1 + labels2, loc='best')
            
            # Save the plot
            plt.tight_layout()
            plt.savefig(f"{output_dir}/energy_vs_timestep.png", dpi=300)
            plt.close()
            
            self.logger.info(f"Energy plot saved to {output_dir}/energy_vs_timestep.png")
            
        except Exception as e:
            self.logger.error(f"Error creating energy plot: {str(e)}")
    
    def _create_convergence_plot(self, dt_values, rel_errors, output_dir):
        """
        Create a log-log plot of relative error vs time step size to show convergence order.
        
        Args:
            dt_values: List of time step values
            rel_errors: List of relative errors
            output_dir: Directory to save the plot
        """
        try:
            # Convert to numpy arrays
            dt_values = np.array(dt_values)
            rel_errors = np.array(rel_errors)
            
            # Filter out NaN or invalid values
            valid_indices = ~np.isnan(rel_errors) & (rel_errors > 0) & (dt_values > 0)
            
            if not np.any(valid_indices):
                self.logger.warning("No valid data points found for convergence analysis.")
                return
            
            # Extract valid data
            valid_dt = dt_values[valid_indices]
            valid_errors = rel_errors[valid_indices]
            
            # Create log-log plot
            plt.figure(figsize=(12, 8))
            
            # Plot data points
            plt.loglog(valid_dt, valid_errors, 'bo', label='Simulation Data')
            
            # Perform linear regression on log-log data to find convergence order
            log_dt = np.log10(valid_dt)
            log_errors = np.log10(valid_errors)
            
            # Filter out any remaining non-finite values
            finite_indices = np.isfinite(log_dt) & np.isfinite(log_errors)
            if not np.any(finite_indices):
                self.logger.warning("No finite data points found for linear regression.")
                return
            
            # Use scipy.stats for linear regression
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                log_dt[finite_indices], 
                log_errors[finite_indices]
            )
            
            # Plot regression line
            x_regression = np.logspace(np.log10(valid_dt.min()), np.log10(valid_dt.max()), 100)
            y_regression = 10**(slope * np.log10(x_regression) + intercept)
            plt.loglog(x_regression, y_regression, 'r-', 
                     label=f'Regression Line (Order = {slope:.2f})')
            
            # Reference lines for 1st and 2nd order convergence
            x_ref = np.logspace(np.log10(valid_dt.min()), np.log10(valid_dt.max()), 100)
            
            # Position reference lines relative to the data
            ref_offset = intercept - np.log10(valid_dt.min())
            
            # First order reference (slope = 1)
            y_first = x_ref * (valid_errors.min() / valid_dt.min())
            plt.loglog(x_ref, y_first, 'g--', label='1st Order Reference (slope=1)')
            
            # Second order reference (slope = 2)
            y_second = x_ref**2 * (valid_errors.min() / valid_dt.min()**2)
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
            plt.savefig(f"{output_dir}/convergence_order.png", dpi=300)
            plt.close()
            
            self.logger.info(f"Convergence plot saved to {output_dir}/convergence_order.png")
            
        except Exception as e:
            self.logger.error(f"Error creating convergence plot: {str(e)}")
    
    def _analyze_stability(self, dt_values, initial_energies, final_energies, rel_errors):
        """
        Analyze the results to determine the largest stable time step.
        
        Args:
            dt_values: List of time step values
            initial_energies: List of initial energy values
            final_energies: List of final energy values
            rel_errors: List of relative errors
        """
        try:
            # Convert to numpy arrays
            dt_values = np.array(dt_values)
            initial_energies = np.array(initial_energies)
            final_energies = np.array(final_energies)
            rel_errors = np.array(rel_errors)
            
            # Find valid data points (where we have both energies and can compute error)
            valid_indices = ~np.isnan(rel_errors)
            
            if not np.any(valid_indices):
                self.logger.warning("No valid data points found for stability analysis.")
                return
            
            # Find largest time step with relative error < 0.1%
            stable_indices = (rel_errors < 0.001) & valid_indices
            if np.any(stable_indices):
                largest_stable_dt = np.max(dt_values[stable_indices])
                self.logger.info(f"Largest stable time step (error < 0.1%): {largest_stable_dt:.6f}")
            else:
                self.logger.info("No time steps meet the 0.1% error threshold.")
            
            # Find largest time step with relative error < 1%
            moderate_indices = (rel_errors < 0.01) & valid_indices
            if np.any(moderate_indices):
                largest_moderate_dt = np.max(dt_values[moderate_indices])
                self.logger.info(f"Largest moderately stable time step (error < 1%): {largest_moderate_dt:.6f}")
            else:
                self.logger.info("No time steps meet the 1% error threshold.")
            
            # Find where instability (NaN values) begins
            nan_indices = np.isnan(final_energies)
            if np.any(nan_indices) and not np.all(nan_indices):
                first_nan_idx = np.where(nan_indices)[0][0]
                if first_nan_idx > 0:
                    instability_dt = dt_values[first_nan_idx]
                    # Recommend a time step well below the instability threshold
                    recommended_dt = min(instability_dt * 0.1, 
                                        largest_moderate_dt if 'largest_moderate_dt' in locals() else 0.01)
                    self.logger.info(f"Instability occurs at time step: {instability_dt:.6f}")
                    self.logger.info(f"Recommended safe time step: {recommended_dt:.6f}")
            
        except Exception as e:
            self.logger.error(f"Error analyzing stability: {str(e)}")
    
    def _analyze_convergence_order(self, dt_values, rel_errors, output_dir):
        """
        Analyze the convergence order of the simulation.
        
        Args:
            dt_values: List of time step values
            rel_errors: List of relative errors
            output_dir: Directory to save the results
        """
        try:
            # Convert to numpy arrays
            dt_values = np.array(dt_values)
            rel_errors = np.array(rel_errors)
            
            # Filter out NaN or invalid values
            valid_indices = ~np.isnan(rel_errors) & (rel_errors > 0) & (dt_values > 0)
            
            if not np.any(valid_indices):
                self.logger.warning("No valid data points found for convergence analysis.")
                return
            
            # Extract valid data
            valid_dt = dt_values[valid_indices]
            valid_errors = rel_errors[valid_indices]
            
            # Take logarithms for regression
            log_dt = np.log10(valid_dt)
            log_errors = np.log10(valid_errors)
            
            # Filter out any remaining non-finite values
            finite_indices = np.isfinite(log_dt) & np.isfinite(log_errors)
            if not np.any(finite_indices):
                self.logger.warning("No finite data points found for linear regression.")
                return
            
            # Use scipy.stats for linear regression
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                log_dt[finite_indices], 
                log_errors[finite_indices]
            )
            
            # Write convergence results to file
            with open(f"{output_dir}/convergence_results.txt", 'w') as f:
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
            
            self.logger.info(f"Experimental convergence order: {slope:.4f} (R² = {r_value**2:.4f})")
            self.logger.info(f"Detailed results saved to {output_dir}/convergence_results.txt")
            
        except Exception as e:
            self.logger.error(f"Error analyzing convergence order: {str(e)}")

# Example usage as a script
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Test N-Body simulations with different time steps")
    parser.add_argument("--script", default="create_initial_conditions.py", 
                       help="Python script for generating initial conditions")
    parser.add_argument("--N", type=int, default=5, help="Number of bodies")
    parser.add_argument("--min-mass", type=float, default=1.0, help="Minimum body mass")
    parser.add_argument("--max-mass", type=float, default=1000000.0, help="Maximum body mass")
    parser.add_argument("--start-dt", type=float, default=0.001, help="Starting time step")
    parser.add_argument("--end-dt", type=float, default=2, help="Ending time step")
    parser.add_argument("--dt-increment", type=float, default=0.001, help="Time step increment")
    parser.add_argument("--final-time", type=float, default=10.0, help="Final simulation time")
    parser.add_argument("--executable", default="./step-0", help="Executable name")
    parser.add_argument("--output", default="TimeStepTest", help="Output directory")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    
    args = parser.parse_args()
    
    # Create tester with appropriate log level
    log_level = logging.DEBUG if args.debug else logging.INFO
    tester = TimeStepEnergyTester(log_level=log_level)
    
    # Run test series
    tester.run_test_series(
        python_script=args.script,
        N=args.N,
        min_mass=args.min_mass,
        max_mass=args.max_mass,
        start_dt=args.start_dt,
        end_dt=args.end_dt,
        dt_increment=args.dt_increment,
        final_time=args.final_time,
        executable_name=args.executable,
        output_dir=args.output
    )
import subprocess
import os
import shutil
import argparse
from pathlib import Path
import logging
import csv
from matplotlib import pyplot as plt
import numpy as np

class TimeStepEnergyTester:
    """
    A class to test N-Body simulations with different time steps and monitor energy conservation
    compared to a truth simulation with a very small time step.
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
                        N=5,
                        min_mass=1000.0,
                        max_mass=1000000.0,
                        truth_dt=0.0001,  # Very small time step for truth simulation
                        start_dt=0.001,
                        end_dt=1.5,
                        dt_increment=0.0001,
                        final_time=60.0,
                        executable_name="./step-4",
                        output_dir="TimeStepTest",
                        result_file="timestep_results.csv"):
        """
        Run a series of simulations with increasing time steps and compare to truth simulation.
        
        Args:
            python_script: Path to the script that generates initial conditions
            N: Number of bodies in the simulation
            min_mass: Minimum mass for bodies
            max_mass: Maximum mass for bodies
            truth_dt: Time step for truth simulation (very small)
            start_dt: Starting time step size for test simulations
            end_dt: Ending time step size for test simulations
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
            writer.writerow(["TimeStep", "InitialEnergy", "FinalEnergy", "TruthDeviation", "RelativeError"])
        
        # Calculate time steps to test
        time_steps = np.arange(start_dt, end_dt + dt_increment/2, dt_increment)
        total_tests = len(time_steps)
        
        self.logger.info(f"Starting test series with {total_tests} different time steps")
        self.logger.info(f"Time steps from {start_dt} to {end_dt} with increment {dt_increment}")
        
        # First, run the truth simulation
        self.logger.info(f"Running truth simulation with dt={truth_dt}")
        truth_initial_energy, truth_final_energy = self._run_single_test(
            python_script=python_script,
            N=N,
            min_mass=min_mass,
            max_mass=max_mass,
            dt=truth_dt,
            final_time=final_time,
            executable_name=executable_name,
            output_dir=os.path.join(output_dir, "truth_simulation")
        )
        
        if truth_initial_energy is None or truth_final_energy is None:
            self.logger.error("Truth simulation failed. Cannot proceed with comparisons.")
            return
            
        self.logger.info(f"Truth simulation completed. Initial energy: {truth_initial_energy}, Final energy: {truth_final_energy}")
        
        # Store data for plotting
        dt_values = []
        initial_energy_values = []
        final_energy_values = []
        truth_deviations = []  # Deviation from truth simulation
        relative_errors = []   # Relative to initial energy (original metric)
        
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
                output_dir=os.path.join(output_dir, f"dt_{dt:.6f}")
            )
            
            # Calculate relative error (compared to own initial energy)
            if initial_energy is not None and final_energy is not None and initial_energy != 0:
                rel_error = abs((final_energy - initial_energy) / initial_energy)
            else:
                rel_error = float('nan')
                
            # Calculate deviation from truth simulation
            if final_energy is not None and truth_final_energy != 0:
                truth_deviation = abs((final_energy - truth_final_energy) / truth_final_energy)
            else:
                truth_deviation = float('nan')
            
            # Store for plotting
            dt_values.append(dt)
            initial_energy_values.append(initial_energy)
            final_energy_values.append(final_energy)
            truth_deviations.append(truth_deviation)
            relative_errors.append(rel_error)
            
            # Write results to CSV
            with open(results_path, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([dt, initial_energy, final_energy, truth_deviation, rel_error])
            
            # Print simple output
            print(f"Time Step: {dt:.6f}, Initial Energy: {initial_energy}, Final Energy: {final_energy}, Truth Deviation: {truth_deviation:.6e}")
        
        self.logger.info(f"Test series completed. Results saved to {results_path}")
        
        # Create plots comparing to truth simulation
        self._create_truth_comparison_plot(dt_values, final_energy_values, truth_final_energy, 
                                         truth_deviations, output_dir, truth_dt)
        
        # Create logarithmic plot
        self._create_log_truth_comparison_plot(dt_values, final_energy_values, truth_final_energy, 
                                             truth_deviations, output_dir, truth_dt)
                                   
        # Create original energy conservation plots (for comparison)
        self._create_energy_plot(dt_values, initial_energy_values, final_energy_values, 
                               relative_errors, output_dir)
        self._create_log_energy_plot(dt_values, initial_energy_values, final_energy_values, 
                                   relative_errors, output_dir)
        
        # Analyze and recommend stable time step
        self._analyze_truth_stability(dt_values, truth_deviations)
    
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
            test_dir = output_dir
            os.makedirs(test_dir, exist_ok=True)
            
            # Generate the shell script with create_initial_conditions.py
            python_cmd = [
                "python3", python_script,
                "--final-time", str(final_time),
                "--snapshots", "5",
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
    
    def _create_truth_comparison_plot(self, dt_values, final_energies, truth_energy, truth_deviations, output_dir, truth_dt):
        """
        Create a plot comparing final energy values to truth simulation (linear scale).
        
        Args:
            dt_values: List of time step values
            final_energies: List of final energy values
            truth_energy: Final energy from truth simulation
            truth_deviations: List of deviations from truth
            output_dir: Directory to save the plot
            truth_dt: Time step used for truth simulation
        """
        try:
            # Convert to numpy arrays
            dt_values = np.array(dt_values)
            final_energies = np.array(final_energies)
            truth_deviations = np.array(truth_deviations)
            
            # Create the plot
            plt.figure(figsize=(12, 8))
            
            # Plot final energy values
            plt.plot(dt_values, final_energies, 'r-', label='Final Energy')
            
            # Plot truth energy as horizontal line
            plt.axhline(y=truth_energy, color='g', linestyle='-', label=f'Truth Energy (dt={truth_dt:.6f})')
            
            # Use scientific notation for y-axis if values are very large
            if np.nanmax(np.abs(final_energies)) > 1000:
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            
            # Plot deviation on secondary axis if we have valid values
            valid_indices = ~np.isnan(truth_deviations)
            if np.any(valid_indices):
                ax2 = plt.gca().twinx()
                ax2.plot(dt_values[valid_indices], truth_deviations[valid_indices], 'b--', label='Deviation from Truth')
                ax2.set_ylabel('Relative Deviation from Truth')
                ax2.set_yscale('log')
                
                # Add reference lines for 0.1%, 1%, and 10% deviation thresholds
                ax2.axhline(y=0.001, color='g', linestyle=':', alpha=0.5, label='0.1% Deviation')
                ax2.axhline(y=0.01, color='y', linestyle=':', alpha=0.5, label='1% Deviation')
                ax2.axhline(y=0.1, color='r', linestyle=':', alpha=0.5, label='10% Deviation')
                
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
            plt.title('Energy vs Time Step Size - Comparison to Truth (Linear Scale)')
            
            # Combine legends from both axes
            lines1, labels1 = plt.gca().get_legend_handles_labels()
            plt.legend(lines1 + lines2, labels1 + labels2, loc='best')
            
            # Save the plot
            plt.tight_layout()
            plt.savefig(f"{output_dir}/energy_vs_timestep_truth_linear.png", dpi=300)
            plt.close()
            
            self.logger.info(f"Truth comparison plot saved to {output_dir}/energy_vs_timestep_truth_linear.png")
            
        except Exception as e:
            self.logger.error(f"Error creating truth comparison plot: {str(e)}")

    def _create_log_truth_comparison_plot(self, dt_values, final_energies, truth_energy, truth_deviations, output_dir, truth_dt):
        """
        Create a logarithmic plot comparing deviations from truth with gradient analysis.
        
        Args:
            dt_values: List of time step values
            final_energies: List of final energy values
            truth_energy: Final energy from truth simulation
            truth_deviations: List of deviations from truth
            output_dir: Directory to save the plot
            truth_dt: Time step used for truth simulation
        """
        try:
            # Convert to numpy arrays
            dt_values = np.array(dt_values)
            truth_deviations = np.array(truth_deviations)
            
            # Find valid data points
            valid_indices = ~np.isnan(truth_deviations)
            
            if not np.any(valid_indices):
                self.logger.warning("No valid data points for plotting")
                return
                
            # Filter to only valid values
            valid_dt = dt_values[valid_indices]
            valid_deviations = truth_deviations[valid_indices]
            
            # Create the plot with log-log scale
            plt.figure(figsize=(12, 8))
            plt.loglog(valid_dt, valid_deviations, 'b-o', label='Deviation from Truth')
            
            # Calculate gradient between consecutive points
            if len(valid_dt) > 1:
                # Calculate slope between points in log-log space
                log_dt = np.log10(valid_dt)
                log_dev = np.log10(valid_deviations)
                
                # Calculate gradient (first derivative)
                gradients = np.diff(log_dev) / np.diff(log_dt)
                
                # Plot gradient
                gradient_dt = valid_dt[:-1] + np.diff(valid_dt)/2  # Center points
                plt.loglog(gradient_dt, np.abs(gradients), 'r--', label='Slope (absolute value)')
                
                # Add average gradient annotation
                avg_gradient = np.mean(gradients)
                plt.text(0.05, 0.05, f'Average Slope: {avg_gradient:.2f}', 
                        transform=plt.gca().transAxes, bbox=dict(facecolor='white', alpha=0.7))
            
            # Add reference lines
            plt.axhline(y=0.001, color='g', linestyle=':', alpha=0.5, label='0.1% Deviation')
            plt.axhline(y=0.01, color='y', linestyle=':', alpha=0.5, label='1% Deviation')
            plt.axhline(y=0.1, color='r', linestyle=':', alpha=0.5, label='10% Deviation')
            
            # Add labels and grid
            plt.xlabel('Time Step Size (dt) - Log Scale')
            plt.ylabel('Relative Deviation from Truth - Log Scale')
            plt.grid(True, which="both", ls="-")
            plt.title('Deviation vs Time Step Size (Log-Log Scale)')
            
            # Add legend
            plt.legend(loc='best')
            
            # Save the plot
            plt.tight_layout()
            plt.savefig(f"{output_dir}/deviation_gradient_log.png", dpi=300)
            plt.close()
            
            self.logger.info(f"Log-log deviation plot with gradient saved to {output_dir}/deviation_gradient_log.png")
            
        except Exception as e:
            self.logger.error(f"Error creating logarithmic deviation plot: {str(e)}")

    def _create_log_truth_comparison_plot_old(self, dt_values, final_energies, truth_energy, truth_deviations, output_dir, truth_dt):
        """
        Create a logarithmic plot comparing final energy values to truth simulation.
        
        Args:
            dt_values: List of time step values
            final_energies: List of final energy values
            truth_energy: Final energy from truth simulation
            truth_deviations: List of deviations from truth
            output_dir: Directory to save the plot
            truth_dt: Time step used for truth simulation
        """
        try:
            # Convert to numpy arrays
            dt_values = np.array(dt_values)
            final_energies = np.array(final_energies)
            truth_deviations = np.array(truth_deviations)
            
            # Create the plot
            plt.figure(figsize=(12, 8))
            
            # Set up logarithmic x-axis
            plt.xscale('log')
            
            # Plot final energy values
            plt.plot(dt_values, final_energies, 'r-', label='Final Energy')
            
            # Plot truth energy as horizontal line
            plt.axhline(y=truth_energy, color='g', linestyle='-', label=f'Truth Energy (dt={truth_dt:.6f})')
            
            # Use scientific notation for y-axis if values are very large
            if np.nanmax(np.abs(final_energies)) > 1000:
                plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            
            # Plot deviation on secondary axis if we have valid values
            valid_indices = ~np.isnan(truth_deviations)
            if np.any(valid_indices):
                ax2 = plt.gca().twinx()
                ax2.plot(dt_values[valid_indices], truth_deviations[valid_indices], 'b--', label='Deviation from Truth')
                ax2.set_ylabel('Relative Deviation from Truth')
                ax2.set_yscale('log')
                
                # Add reference lines for error thresholds
                ax2.axhline(y=0.001, color='g', linestyle=':', alpha=0.5, label='0.1% Deviation')
                ax2.axhline(y=0.01, color='y', linestyle=':', alpha=0.5, label='1% Deviation')
                ax2.axhline(y=0.1, color='r', linestyle=':', alpha=0.5, label='10% Deviation')
                
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
            plt.xlabel('Time Step Size (dt) - Log Scale')
            plt.ylabel('Energy')
            plt.grid(True, which="both", ls="-")
            plt.title('Energy vs Time Step Size - Comparison to Truth (Logarithmic Scale)')
            
            # Combine legends from both axes
            lines1, labels1 = plt.gca().get_legend_handles_labels()
            plt.legend(lines1 + lines2, labels1 + labels2, loc='best')
            
            # Save the plot
            plt.tight_layout()
            plt.savefig(f"{output_dir}/energy_vs_timestep_truth_log.png", dpi=300)
            plt.close()
            
            self.logger.info(f"Logarithmic truth comparison plot saved to {output_dir}/energy_vs_timestep_truth_log.png")
            
        except Exception as e:
            self.logger.error(f"Error creating logarithmic truth comparison plot: {str(e)}")
    
    def _create_energy_plot(self, dt_values, initial_energies, final_energies, rel_errors, output_dir):
        """
        Create a plot of energy values vs time step size (linear scale).
        
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
            plt.title('Energy vs Time Step Size (Linear Scale)')
            
            # Combine legends from both axes
            lines1, labels1 = plt.gca().get_legend_handles_labels()
            plt.legend(lines1 + lines2, labels1 + labels2, loc='best')
            
            # Save the plot
            plt.tight_layout()
            plt.savefig(f"{output_dir}/energy_vs_timestep_linear.png", dpi=300)
            plt.close()
            
            self.logger.info(f"Linear scale energy plot saved to {output_dir}/energy_vs_timestep_linear.png")
            
        except Exception as e:
            self.logger.error(f"Error creating linear energy plot: {str(e)}")

    def _create_log_energy_plot(self, dt_values, initial_energies, final_energies, rel_errors, output_dir):
        """
        Create a logarithmic plot of energy values vs time step size.
        
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
            
            # Set up logarithmic x-axis
            plt.xscale('log')
            
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
                
                # Add reference lines for 0.1%, 1%, and 10% error thresholds
                ax2.axhline(y=0.001, color='g', linestyle=':', alpha=0.5, label='0.1% Error')
                ax2.axhline(y=0.01, color='y', linestyle=':', alpha=0.5, label='1% Error')
                ax2.axhline(y=0.1, color='r', linestyle=':', alpha=0.5, label='10% Error')
                
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
            plt.xlabel('Time Step Size (dt) - Log Scale')
            plt.ylabel('Energy')
            plt.grid(True, which="both", ls="-")
            plt.title('Energy vs Time Step Size (Logarithmic Scale)')
            
            # Combine legends from both axes
            lines1, labels1 = plt.gca().get_legend_handles_labels()
            plt.legend(lines1 + lines2, labels1 + labels2, loc='best')
            
            # Save the plot
            plt.tight_layout()
            plt.savefig(f"{output_dir}/energy_vs_timestep_log.png", dpi=300)
            plt.close()
            
            self.logger.info(f"Logarithmic scale energy plot saved to {output_dir}/energy_vs_timestep_log.png")
            
        except Exception as e:
            self.logger.error(f"Error creating logarithmic energy plot: {str(e)}")
    
    def _analyze_truth_stability(self, dt_values, truth_deviations):
        """
        Analyze the results to determine the largest stable time step based on truth comparison.
        
        Args:
            dt_values: List of time step values
            truth_deviations: List of deviations from truth simulation
        """
        try:
            # Convert to numpy arrays
            dt_values = np.array(dt_values)
            truth_deviations = np.array(truth_deviations)
            
            # Find valid data points
            valid_indices = ~np.isnan(truth_deviations)
            
            if not np.any(valid_indices):
                self.logger.warning("No valid data points found for stability analysis.")
                return
            
            # Find largest time step with deviation < 0.1%
            stable_indices = (truth_deviations < 0.001) & valid_indices
            if np.any(stable_indices):
                largest_stable_dt = np.max(dt_values[stable_indices])
                self.logger.info(f"Largest stable time step (deviation < 0.1%): {largest_stable_dt:.6f}")
            else:
                self.logger.info("No time steps meet the 0.1% deviation threshold.")
            
            # Find largest time step with deviation < 1%
            moderate_indices = (truth_deviations < 0.01) & valid_indices
            if np.any(moderate_indices):
                largest_moderate_dt = np.max(dt_values[moderate_indices])
                self.logger.info(f"Largest moderately stable time step (deviation < 1%): {largest_moderate_dt:.6f}")
            else:
                self.logger.info("No time steps meet the 1% deviation threshold.")
            
            # Find largest time step with deviation < 10%
            acceptable_indices = (truth_deviations < 0.1) & valid_indices
            if np.any(acceptable_indices):
                largest_acceptable_dt = np.max(dt_values[acceptable_indices])
                self.logger.info(f"Largest acceptable time step (deviation < 10%): {largest_acceptable_dt:.6f}")
            else:
                self.logger.info("No time steps meet the 10% deviation threshold.")
            
            # Find where instability (NaN values) begins
            nan_indices = np.isnan(truth_deviations)
            if np.any(nan_indices) and not np.all(nan_indices):
                first_nan_idx = np.where(nan_indices)[0][0]
                if first_nan_idx > 0:
                    instability_dt = dt_values[first_nan_idx]
                    self.logger.info(f"Instability occurs at time step: {instability_dt:.6f}")
                    
                    # Recommend a time step well below the instability threshold
                    if 'largest_moderate_dt' in locals():
                        recommended_dt = largest_moderate_dt
                    elif 'largest_acceptable_dt' in locals():
                        recommended_dt = largest_acceptable_dt * 0.5  # More conservative
                    else:
                        recommended_dt = instability_dt * 0.1  # Very conservative
                        
                    self.logger.info(f"Recommended safe time step: {recommended_dt:.6f}")
            
        except Exception as e:
            self.logger.error(f"Error analyzing stability: {str(e)}")


# Example usage as a script
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Test N-Body simulations with different time steps")
    parser.add_argument("--script", default="create_initial_conditions.py", 
                       help="Python script for generating initial conditions")
    parser.add_argument("--N", type=int, default=5, help="Number of bodies")
    parser.add_argument("--min-mass", type=float, default=1000.0, help="Minimum body mass")
    parser.add_argument("--max-mass", type=float, default=1000000.0, help="Maximum body mass")
    parser.add_argument("--truth-dt", type=float, default=0.0001, help="Time step for truth simulation")
    parser.add_argument("--start-dt", type=float, default=0.001, help="Starting time step")
    parser.add_argument("--end-dt", type=float, default=1.5, help="Ending time step")
    parser.add_argument("--dt-increment", type=float, default=0.0001, help="Time step increment")
    parser.add_argument("--final-time", type=float, default=10.0, help="Final simulation time")
    parser.add_argument("--executable", default="./step-4", help="Executable name")
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
        truth_dt=args.truth_dt,
        start_dt=args.start_dt,
        end_dt=args.end_dt,
        dt_increment=args.dt_increment,
        final_time=args.final_time,
        executable_name=args.executable,
        output_dir=args.output
    )
import subprocess
import os
import shutil
from pathlib import Path
import logging
import csv
import time
import numpy as np
from matplotlib import pyplot as plt

class NBodyScalingTester:
    """
    A class to test N-Body simulations with increasing N values and monitor execution time.
    """
    
    def __init__(self, log_level=logging.INFO):
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger('NBodyScalingTester')
    
    def run_scaling_tests(self, 
                          python_script="create_initial_conditions.py",
                          start_N=1,
                          max_N=20,
                          N_step=1,
                          min_mass=1.0,
                          max_mass=1000000.0,
                          dt=0.001,
                          final_time=600.0,
                          executable_name="./step-2",
                          output_dir="ScalingTest",
                          result_file="scaling_results.csv"):
        """
        Run simulations with increasing N and measure execution time.
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        results_path = output_path / result_file
        with open(results_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["N", "ExecutionTime"])
        
        N_values = list(range(start_N, max_N + 1, N_step))
        execution_times = []
        
        for N in N_values:
            self.logger.info(f"Running test for N = {N}")
            execution_time = self._run_single_test(python_script, N, min_mass, max_mass, dt, final_time, executable_name, output_dir)
            execution_times.append(execution_time)
            
            with open(results_path, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([N, execution_time])
            
            print(f"N: {N}, Execution Time: {execution_time:.6f} seconds")
        
        self._create_execution_plot(N_values, execution_times, output_dir)
    
    def _run_single_test(self, python_script, N, min_mass, max_mass, dt, final_time, executable_name, output_dir):
        """
        Run a single simulation test and measure execution time.
        """
        try:
            test_dir = f"{output_dir}/N_{N}"
            os.makedirs(test_dir, exist_ok=True)
            
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
            
            subprocess.run(python_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            shell_script = f"{os.path.basename(executable_name)}.sh"
            os.chmod(shell_script, 0o755)
            
            start_time = time.time()
            subprocess.run(f"bash ./{shell_script}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            execution_time = time.time() - start_time
            
            return execution_time
        except Exception as e:
            self.logger.error(f"Error in test with N={N}: {str(e)}")
            return None
    
    def _create_execution_plot(self, N_values, execution_times, output_dir):
        """
        Create a plot of execution time vs N.
        """
        try:
            plt.figure(figsize=(10, 6))
            plt.plot(N_values, execution_times, 'bo-', label='Execution Time')
            plt.xlabel('Number of Bodies (N)')
            plt.ylabel('Execution Time (seconds)')
            plt.title('Execution Time vs N')
            plt.grid(True)
            plt.legend()
            plt.savefig(f"{output_dir}/execution_time_vs_N.png", dpi=300)
            plt.close()
            
            self.logger.info(f"Execution time plot saved to {output_dir}/execution_time_vs_N.png")
        except Exception as e:
            self.logger.error(f"Error creating execution time plot: {str(e)}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Test N-Body simulations with increasing N values")
    parser.add_argument("--script", default="create_initial_conditions.py", help="Python script for generating initial conditions")
    parser.add_argument("--start-N", type=int, default=1, help="Starting number of bodies")
    parser.add_argument("--max-N", type=int, default=20, help="Maximum number of bodies")
    parser.add_argument("--N-step", type=int, default=1, help="Step increment for N")
    parser.add_argument("--min-mass", type=float, default=1.0, help="Minimum body mass")
    parser.add_argument("--max-mass", type=float, default=1000000.0, help="Maximum body mass")
    parser.add_argument("--dt", type=float, default=0.001, help="Time step size")
    parser.add_argument("--final-time", type=float, default=600.0, help="Final simulation time")
    parser.add_argument("--executable", default="./step-2", help="Executable name")
    parser.add_argument("--output", default="ScalingTest", help="Output directory")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    
    args = parser.parse_args()
    log_level = logging.DEBUG if args.debug else logging.INFO
    tester = NBodyScalingTester(log_level=log_level)
    
    tester.run_scaling_tests(
        python_script=args.script,
        start_N=args.start_N,
        max_N=args.max_N,
        N_step=args.N_step,
        min_mass=args.min_mass,
        max_mass=args.max_mass,
        dt=args.dt,
        final_time=args.final_time,
        executable_name=args.executable,
        output_dir=args.output
    )

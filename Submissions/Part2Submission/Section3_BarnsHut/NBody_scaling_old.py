# # import subprocess
# # import os
# # import shutil
# # from pathlib import Path
# # import logging
# # import csv
# # import time
# # import numpy as np
# # from matplotlib import pyplot as plt

# # class NBodyScalingTester:
# #     """
# #     A class to test N-Body simulations with increasing N values and monitor execution time.
# #     """
    
# #     def __init__(self, log_level=logging.INFO):
# #         logging.basicConfig(
# #             level=log_level,
# #             format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
# #         )
# #         self.logger = logging.getLogger('NBodyScalingTester')
    
# #     def run_scaling_tests(self, 
# #                           python_script="create_initial_conditions.py",
# #                           start_N=1,
# #                           max_N=10,
# #                           N_step=1,
# #                           min_mass=1000.0,
# #                           max_mass=1000000.0,
# #                           dt=0.01,
# #                           final_time=600.0,
# #                           executable_name="./step-4",
# #                           output_dir="ScalingTest",
# #                           result_file="scaling_results.csv"):
# #         """
# #         Run simulations with increasing N and measure execution time.
# #         """
# #         output_path = Path(output_dir)
# #         output_path.mkdir(parents=True, exist_ok=True)
        
# #         results_path = output_path / result_file
# #         with open(results_path, 'w', newline='') as f:
# #             writer = csv.writer(f)
# #             writer.writerow(["N", "ExecutionTime"])
        
# #         N_values = list(range(start_N, max_N + 1, N_step))
# #         execution_times = []
        
# #         for N in N_values:
# #             self.logger.info(f"Running test for N = {N}")
# #             execution_time = self._run_single_test(python_script, N, min_mass, max_mass, dt, final_time, executable_name, output_dir)
# #             execution_times.append(execution_time)
            
# #             with open(results_path, 'a', newline='') as f:
# #                 writer = csv.writer(f)
# #                 writer.writerow([N, execution_time])
            
# #             print(f"N: {N}, Execution Time: {execution_time:.6f} seconds")
        
# #         self._create_execution_plot(N_values, execution_times, output_dir)
    
# #     def _run_single_test(self, python_script, N, min_mass, max_mass, dt, final_time, executable_name, output_dir):
# #         """
# #         Run a single simulation test and measure execution time.
# #         """
# #         try:
# #             test_dir = f"{output_dir}/N_{N}"
# #             os.makedirs(test_dir, exist_ok=True)
            
# #             python_cmd = [
# #                 "python3", python_script,
# #                 "--final-time", str(final_time),
# #                 "--snapshots", "0",
# #                 "--executable-name", executable_name,
# #                 "--min-mass", str(min_mass),
# #                 "--max-mass", str(max_mass),
# #                 "--dt", str(dt),
# #                 "--N", str(N)
# #             ]
            
# #             subprocess.run(python_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
# #             shell_script = f"{os.path.basename(executable_name)}.sh"
# #             os.chmod(shell_script, 0o755)
            
# #             start_time = time.time()
# #             subprocess.run(f"bash ./{shell_script}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# #             execution_time = time.time() - start_time
            
# #             return execution_time
# #         except Exception as e:
# #             self.logger.error(f"Error in test with N={N}: {str(e)}")
# #             return None
    
# #     def _create_execution_plot(self, N_values, execution_times, output_dir):
# #         """
# #         Create a plot of execution time vs N.
# #         """
# #         try:
# #             plt.figure(figsize=(10, 6))
# #             plt.plot(N_values, execution_times, 'bo-', label='Execution Time')
# #             plt.xlabel('Number of Bodies (N)')
# #             plt.ylabel('Execution Time (seconds)')
# #             plt.title('Execution Time vs N')
# #             plt.grid(True)
# #             plt.legend()
# #             plt.savefig(f"{output_dir}/execution_time_vs_N.png", dpi=300)
# #             plt.close()
            
# #             self.logger.info(f"Execution time plot saved to {output_dir}/execution_time_vs_N.png")
# #         except Exception as e:
# #             self.logger.error(f"Error creating execution time plot: {str(e)}")

# # if __name__ == "__main__":
# #     import argparse
    
# #     parser = argparse.ArgumentParser(description="Test N-Body simulations with increasing N values")
# #     parser.add_argument("--script", default="create_initial_conditions.py", help="Python script for generating initial conditions")
# #     parser.add_argument("--start-N", type=int, default=1, help="Starting number of bodies")
# #     parser.add_argument("--max-N", type=int, default=20, help="Maximum number of bodies")
# #     parser.add_argument("--N-step", type=int, default=1, help="Step increment for N")
# #     parser.add_argument("--min-mass", type=float, default=1000.0, help="Minimum body mass")
# #     parser.add_argument("--max-mass", type=float, default=1000000.0, help="Maximum body mass")
# #     parser.add_argument("--dt", type=float, default=0.01, help="Time step size")
# #     parser.add_argument("--final-time", type=float, default=600.0, help="Final simulation time")
# #     parser.add_argument("--executable", default="./step-4", help="Executable name")
# #     parser.add_argument("--output", default="ScalingTest", help="Output directory")
# #     parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    
# #     args = parser.parse_args()
# #     log_level = logging.DEBUG if args.debug else logging.INFO
# #     tester = NBodyScalingTester(log_level=log_level)
    
# #     tester.run_scaling_tests(
# #         python_script=args.script,
# #         start_N=args.start_N,
# #         max_N=args.max_N,
# #         N_step=args.N_step,
# #         min_mass=args.min_mass,
# #         max_mass=args.max_mass,
# #         dt=args.dt,
# #         final_time=args.final_time,
# #         executable_name=args.executable,
# #         output_dir=args.output
# #     )

# import subprocess
# import os
# import shutil
# from pathlib import Path
# import logging
# import csv
# import time
# import numpy as np
# from matplotlib import pyplot as plt

# class NBodyScalingTester:
#     """
#     A class to test N-Body simulations with increasing N values and monitor execution time.
#     """
    
#     def __init__(self, log_level=logging.INFO):
#         logging.basicConfig(
#             level=log_level,
#             format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
#         )
#         self.logger = logging.getLogger('NBodyScalingTester')
    
#     def run_scaling_tests(self, 
#                           python_script="create_initial_conditions.py",
#                           start_N=1,
#                           max_N=10,
#                           N_step=1,
#                           min_mass=1000.0,
#                           max_mass=1000000.0,
#                           dt=0.01,
#                           final_time=600.0,
#                           executables=["./step-0", "./step-4"],
#                           output_dir="ScalingTest",
#                           result_file="scaling_results.csv"):
#         """
#         Run simulations with increasing N and measure execution time for multiple executables.
#         """
#         output_path = Path(output_dir)
#         output_path.mkdir(parents=True, exist_ok=True)
        
#         results_path = output_path / result_file
#         with open(results_path, 'w', newline='') as f:
#             writer = csv.writer(f)
#             header = ["N"] + [f"{os.path.basename(exec_name)}_ExecutionTime" for exec_name in executables]
#             writer.writerow(header)
        
#         N_values = list(range(start_N, max_N + 1, N_step))
#         execution_times = {exec_name: [] for exec_name in executables}
        
#         for N in N_values:
#             self.logger.info(f"Running test for N = {N}")
#             times_for_n = []
            
#             for executable in executables:
#                 exec_name = os.path.basename(executable)
#                 self.logger.info(f"Testing with executable: {exec_name}")
#                 execution_time = self._run_single_test(python_script, N, min_mass, max_mass, dt, 
#                                                       final_time, executable, f"{output_dir}/N_{N}_{exec_name}")
#                 execution_times[executable].append(execution_time)
#                 times_for_n.append(execution_time)
                
#                 print(f"N: {N}, Executable: {exec_name}, Execution Time: {execution_time:.6f} seconds")
            
#             with open(results_path, 'a', newline='') as f:
#                 writer = csv.writer(f)
#                 writer.writerow([N] + times_for_n)
        
#         self._create_execution_plots(N_values, execution_times, executables, output_dir)
    
#     def _run_single_test(self, python_script, N, min_mass, max_mass, dt, final_time, executable_name, output_dir):
#         """
#         Run a single simulation test and measure execution time.
#         """
#         try:
#             os.makedirs(output_dir, exist_ok=True)
            
#             # Generate initial conditions for this executable
#             python_cmd = [
#                 "python3", python_script,
#                 "--final-time", str(final_time),
#                 "--snapshots", "0",
#                 "--executable-name", executable_name,
#                 "--min-mass", str(min_mass),
#                 "--max-mass", str(max_mass),
#                 "--dt", str(dt),
#                 "--N", str(N)
#             ]
            
#             subprocess.run(python_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
#             shell_script = f"{os.path.basename(executable_name)}.sh"
#             os.chmod(shell_script, 0o755)
            
#             # Run the simulation and time it
#             start_time = time.time()
#             subprocess.run(f"bash ./{shell_script}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#             execution_time = time.time() - start_time
            
#             # Move output files to this test's directory
#             for file in os.listdir('.'):
#                 if file.startswith(os.path.basename(executable_name)) and file.endswith('.sh'):
#                     shutil.move(file, os.path.join(output_dir, file))
            
#             return execution_time
#         except Exception as e:
#             self.logger.error(f"Error in test with N={N}, executable={executable_name}: {str(e)}")
#             return None
    
#     def _create_execution_plots(self, N_values, execution_times, executables, output_dir):
#         """
#         Create plots comparing execution times between different executables.
#         """
#         try:
#             plt.figure(figsize=(12, 7))
            
#             # Line plot for each executable
#             for executable in executables:
#                 exec_name = os.path.basename(executable)
#                 plt.plot(N_values, execution_times[executable], 'o-', label=f'{exec_name}')
            
#             plt.xlabel('Number of Bodies (N)')
#             plt.ylabel('Execution Time (seconds)')
#             plt.title('Execution Time Comparison: step-0 vs step-4')
#             plt.grid(True)
#             plt.legend()
#             plt.savefig(f"{output_dir}/execution_time_comparison.png", dpi=300)
            
#             # Create bar chart for direct comparison
#             plt.figure(figsize=(12, 7))
#             bar_width = 0.35
#             index = np.arange(len(N_values))
            
#             for i, executable in enumerate(executables):
#                 exec_name = os.path.basename(executable)
#                 plt.bar(index + i*bar_width, execution_times[executable], bar_width, 
#                         label=f'{exec_name}')
            
#             plt.xlabel('Number of Bodies (N)')
#             plt.ylabel('Execution Time (seconds)')
#             plt.title('Execution Time Comparison: step-0 vs step-4')
#             plt.xticks(index + bar_width/2, N_values)
#             plt.legend()
#             plt.savefig(f"{output_dir}/execution_time_bar_comparison.png", dpi=300)
            
#             # Create speedup plot (ratio of step-0 to step-4)
#             if len(executables) >= 2:
#                 plt.figure(figsize=(12, 7))
#                 speedups = [t0/t4 if t4 and t0 else 0 
#                            for t0, t4 in zip(execution_times[executables[0]], 
#                                             execution_times[executables[1]])]
#                 plt.plot(N_values, speedups, 'go-', label='Speedup ratio')
#                 plt.xlabel('Number of Bodies (N)')
#                 plt.ylabel('Speedup (step-0 time / step-4 time)')
#                 plt.title('Performance Speedup: step-0 vs step-4')
#                 plt.grid(True)
#                 plt.legend()
#                 plt.savefig(f"{output_dir}/speedup_ratio.png", dpi=300)
                
#             plt.close('all')
            
#             self.logger.info(f"Execution time plots saved to {output_dir}/")
#         except Exception as e:
#             self.logger.error(f"Error creating execution time plots: {str(e)}")

# if __name__ == "__main__":
#     import argparse
    
#     parser = argparse.ArgumentParser(description="Test N-Body simulations with increasing N values")
#     parser.add_argument("--script", default="create_initial_conditions.py", help="Python script for generating initial conditions")
#     parser.add_argument("--start-N", type=int, default=1, help="Starting number of bodies")
#     parser.add_argument("--max-N", type=int, default=20, help="Maximum number of bodies")
#     parser.add_argument("--N-step", type=int, default=1, help="Step increment for N")
#     parser.add_argument("--min-mass", type=float, default=1000.0, help="Minimum body mass")
#     parser.add_argument("--max-mass", type=float, default=1000000.0, help="Maximum body mass")
#     parser.add_argument("--dt", type=float, default=0.01, help="Time step size")
#     parser.add_argument("--final-time", type=float, default=600.0, help="Final simulation time")
#     parser.add_argument("--output", default="ScalingTest", help="Output directory")
#     parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    
#     args = parser.parse_args()
#     log_level = logging.DEBUG if args.debug else logging.INFO
#     tester = NBodyScalingTester(log_level=log_level)
    
#     tester.run_scaling_tests(
#         python_script=args.script,
#         start_N=args.start_N,
#         max_N=args.max_N,
#         N_step=args.N_step,
#         min_mass=args.min_mass,
#         max_mass=args.max_mass,
#         dt=args.dt,
#         final_time=args.final_time,
#         executables=["./step-0", "./step-4"],
#         output_dir=args.output
#     )

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
                          max_N=10,
                          N_step=1,
                          min_mass=1000.0,
                          max_mass=1000000.0,
                          dt=0.01,
                          final_time=600.0,
                          executables=["./step-0", "./step-4"],
                          output_dir="ScalingTest",
                          result_file="scaling_results.csv",
                          show_output=True):
        """
        Run simulations with increasing N and measure execution time for multiple executables.
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        results_path = output_path / result_file
        with open(results_path, 'w', newline='') as f:
            writer = csv.writer(f)
            header = ["N"] + [f"{os.path.basename(exec_name)}_ExecutionTime" for exec_name in executables]
            writer.writerow(header)
        
        N_values = list(range(start_N, max_N + 1, N_step))
        execution_times = {exec_name: [] for exec_name in executables}
        
        for N in N_values:
            self.logger.info(f"Running test for N = {N}")
            times_for_n = []
            
            for executable in executables:
                exec_name = os.path.basename(executable)
                self.logger.info(f"Testing with executable: {exec_name}")
                execution_time = self._run_single_test(python_script, N, min_mass, max_mass, dt, 
                                                      final_time, executable, f"{output_dir}/N_{N}_{exec_name}")
                if execution_time is not None:
                    execution_times[executable].append(execution_time)
                    times_for_n.append(execution_time)
                    print(f"N: {N}, Executable: {exec_name}, Execution Time: {execution_time:.6f} seconds")
                else:
                    # Add a placeholder for failed runs to maintain data alignment
                    execution_times[executable].append(None)
                    times_for_n.append(None)
                    print(f"N: {N}, Executable: {exec_name}, Execution FAILED")
            
            with open(results_path, 'a', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([N] + times_for_n)
        
        self._create_execution_plots(N_values, execution_times, executables, output_dir)
    
    def _run_single_test(self, python_script, N, min_mass, max_mass, dt, final_time, executable_name, output_dir):
        """
        Run a single simulation test and measure execution time.
        """
        try:
            os.makedirs(output_dir, exist_ok=True)
            
            # Generate the initial conditions run output directory
            self.logger.info(f"Saving results to {output_dir}")
            
            # Generate initial conditions for this executable
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
            
            # Log the command being run
            cmd_str = " ".join(python_cmd)
            print(f"Running: {cmd_str}")
            
            # Save the command to a file
            with open(os.path.join(output_dir, "command.txt"), 'w') as f:
                f.write(f"Command: {cmd_str}\n")
            
            # Execute the initial conditions generation
            init_result = subprocess.run(python_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            # Save the initialization output
            with open(os.path.join(output_dir, "init_output.log"), 'w') as f:
                f.write("STDOUT:\n")
                f.write(init_result.stdout)
                f.write("\nSTDERR:\n")
                f.write(init_result.stderr)
            
            # Check if any errors occurred
            if init_result.returncode != 0:
                self.logger.error(f"Error generating initial conditions: {init_result.stderr}")
                return None
            shell_script = f"{os.path.basename(executable_name)}.sh"
            os.chmod(shell_script, 0o755)
            
            # Run the simulation and time it with timeout protection
            start_time = time.time()
            try:
                # Use a timeout value (adjust as needed)
                timeout_seconds = 3600  # 1 hour timeout
                result = subprocess.run(f"bash ./{shell_script}", shell=True, 
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                      timeout=timeout_seconds,
                                      text=True)
                
                # Save the output to a log file
                log_file = os.path.join(output_dir, f"{os.path.basename(executable_name)}_output.log")
                with open(log_file, 'w') as f:
                    f.write("STDOUT:\n")
                    f.write(result.stdout)
                    f.write("\nSTDERR:\n")
                    f.write(result.stderr)
                
                # Print the first 500 characters of output for quick reference
                print(f"Output (first 500 chars):\n{result.stdout[:500]}")
                if result.stderr:
                    print(f"Errors:\n{result.stderr[:500]}")
                
                # Check if the process completed successfully
                if result.returncode != 0:
                    self.logger.error(f"Process failed with return code {result.returncode}")
                    self.logger.error(f"STDERR: {result.stderr}")
                    return None
                    
                execution_time = time.time() - start_time
            except subprocess.TimeoutExpired as e:
                self.logger.error(f"Process timed out after {timeout_seconds} seconds")
                # Save any partial output that was captured
                log_file = os.path.join(output_dir, f"{os.path.basename(executable_name)}_timeout.log")
                with open(log_file, 'w') as f:
                    f.write("TIMEOUT ERROR\n")
                    if hasattr(e, 'stdout') and e.stdout:
                        f.write("\nPartial STDOUT:\n")
                        f.write(e.stdout)
                    if hasattr(e, 'stderr') and e.stderr:
                        f.write("\nPartial STDERR:\n")
                        f.write(e.stderr)
                return None
            
            # Move output files to this test's directory
            for file in os.listdir('.'):
                if file.startswith(os.path.basename(executable_name)) and file.endswith('.sh'):
                    shutil.move(file, os.path.join(output_dir, file))
            
            return execution_time
        except Exception as e:
            self.logger.error(f"Error in test with N={N}, executable={executable_name}: {str(e)}")
            return None
    
    def _create_execution_plots(self, N_values, execution_times, executables, output_dir):
        """
        Create plots comparing execution times between different executables.
        """
        try:
            plt.figure(figsize=(12, 7))
            
            # Line plot for each executable
            for executable in executables:
                exec_name = os.path.basename(executable)
                # Filter out None values for plotting
                valid_points = [(n, t) for n, t in zip(N_values, execution_times[executable]) if t is not None]
                if valid_points:
                    x_vals, y_vals = zip(*valid_points)
                    plt.plot(x_vals, y_vals, 'o-', label=f'{exec_name}')
            
            plt.xlabel('Number of Bodies (N)')
            plt.ylabel('Execution Time (seconds)')
            plt.title('Execution Time Comparison: step-0 vs step-4')
            plt.grid(True)
            plt.legend()
            plt.savefig(f"{output_dir}/execution_time_comparison.png", dpi=300)
            
            # Create bar chart for direct comparison
            plt.figure(figsize=(12, 7))
            bar_width = 0.35
            index = np.arange(len(N_values))
            
            for i, executable in enumerate(executables):
                exec_name = os.path.basename(executable)
                # Replace None values with zeros and create a mask for valid data
                times = execution_times[executable].copy()
                valid_mask = [t is not None for t in times]
                times = [t if t is not None else 0 for t in times]
                
                bars = plt.bar(index + i*bar_width, times, bar_width, label=f'{exec_name}')
                
                # Mark failed runs in the bar chart
                for j, valid in enumerate(valid_mask):
                    if not valid:
                        bars[j].set_hatch('//')
                        bars[j].set_facecolor('lightgray')
            
            plt.xlabel('Number of Bodies (N)')
            plt.ylabel('Execution Time (seconds)')
            plt.title('Execution Time Comparison: step-0 vs step-4')
            plt.xticks(index + bar_width/2, N_values)
            plt.legend()
            plt.savefig(f"{output_dir}/execution_time_bar_comparison.png", dpi=300)
            
            # Create speedup plot (ratio of step-0 to step-4)
            if len(executables) >= 2:
                plt.figure(figsize=(12, 7))
                # Calculate speedups only for valid data points
                speedups = []
                valid_n_values = []
                
                for i, (t0, t4) in enumerate(zip(execution_times[executables[0]], 
                                              execution_times[executables[1]])):
                    if t0 is not None and t4 is not None and t4 > 0:
                        speedups.append(t0/t4)
                        valid_n_values.append(N_values[i])
                
                if speedups:
                    plt.plot(valid_n_values, speedups, 'go-', label='Speedup ratio')
                plt.xlabel('Number of Bodies (N)')
                plt.ylabel('Speedup (step-0 time / step-4 time)')
                plt.title('Performance Speedup: step-0 vs step-4')
                plt.grid(True)
                plt.legend()
                plt.savefig(f"{output_dir}/speedup_ratio.png", dpi=300)
                
            plt.close('all')
            
            self.logger.info(f"Execution time plots saved to {output_dir}/")
        except Exception as e:
            self.logger.error(f"Error creating execution time plots: {str(e)}")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Test N-Body simulations with increasing N values")
    parser.add_argument("--script", default="create_initial_conditions.py", help="Python script for generating initial conditions")
    parser.add_argument("--start-N", type=int, default=1, help="Starting number of bodies")
    parser.add_argument("--max-N", type=int, default=20, help="Maximum number of bodies")
    parser.add_argument("--N-step", type=int, default=1, help="Step increment for N")
    parser.add_argument("--min-mass", type=float, default=1000.0, help="Minimum body mass")
    parser.add_argument("--max-mass", type=float, default=1000000.0, help="Maximum body mass")
    parser.add_argument("--dt", type=float, default=0.01, help="Time step size")
    parser.add_argument("--final-time", type=float, default=600.0, help="Final simulation time")
    parser.add_argument("--output", default="ScalingTest", help="Output directory")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    parser.add_argument("--executables", nargs='+', default=["./step-0", "./step-4"], 
                        help="List of executables to test")
    parser.add_argument("--hide-output", action="store_true", help="Hide stdout/stderr output")
    
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
        executables=args.executables,
        output_dir=args.output,
        show_output=not args.hide_output
    )
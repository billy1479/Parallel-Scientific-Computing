import re
import subprocess
import matplotlib.pyplot as plt
import numpy as np

class NBodyMetricsTracker:
    """
    Tracks and plots metrics from N-Body simulation output, specifically:
    - Maximum velocity (v_max)
    - Minimum distance (dx_min)
    - Time step size (dt)
    """
    
    def __init__(self):
        self.time_steps = []
        self.times = []
        self.dts = []
        self.v_maxs = []
        self.dx_mins = []
        
    def run_simulation(self, executable_path, script_path):
        """Run the simulation and capture output metrics"""
        try:
            result = subprocess.run(
                f"bash {script_path}",
                shell=True,
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse the output
            self._parse_output(result.stdout)
            return True
            
        except subprocess.CalledProcessError as e:
            print(f"Error running simulation: {e}")
            print(f"stderr: {e.stderr}")
            return False
    
    def extract_from_existing_output(self, output_file):
        """Extract metrics from an existing output file"""
        try:
            with open(output_file, 'r') as f:
                output = f.read()
            
            # Parse the output
            self._parse_output(output)
            return True
            
        except Exception as e:
            print(f"Error extracting metrics: {e}")
            return False
    
    def _parse_output(self, output):
        """Parse simulation output to extract metrics"""
        # Pattern to match the snapshot lines
        pattern = r"plot next snapshot,\s*time step=(\d+),\s*t=([0-9.]+),\s*dt=([0-9.]+),\s*v_max=([0-9.]+),\s*dx_min=([0-9.]+)"
        
        # Find all matches
        matches = re.findall(pattern, output)
        
        # Extract data
        for match in matches:
            self.time_steps.append(int(match[0]))
            self.times.append(float(match[1]))
            self.dts.append(float(match[2]))
            self.v_maxs.append(float(match[3]))
            self.dx_mins.append(float(match[4]))
        
        print(f"Extracted {len(self.time_steps)} data points")
    
    def calculate_cfl_condition(self):
        """
        Calculate the CFL condition (dx_min/v_max) which gives theoretical max dt
        """
        if not self.v_maxs or not self.dx_mins:
            print("No data available to calculate CFL condition")
            return []
        
        cfl_values = []
        for dx, v in zip(self.dx_mins, self.v_maxs):
            if v > 0:  # Avoid division by zero
                cfl_values.append(dx / v)
            else:
                cfl_values.append(float('nan'))
        
        return cfl_values
    
    def plot_metrics(self, output_file='nbody_metrics.png'):
        """Plot the tracked metrics"""
        if not self.time_steps:
            print("No data to plot")
            return False
        
        # Calculate CFL condition
        cfl_values = self.calculate_cfl_condition()
        
        # Create figure with subplots
        fig, axs = plt.subplots(4, 1, figsize=(10, 16), sharex=True)
        
        # Plot v_max
        axs[0].plot(self.times, self.v_maxs, 'r-', label='Maximum Velocity')
        axs[0].set_ylabel('v_max')
        axs[0].set_title('Maximum Velocity vs Time')
        axs[0].grid(True)
        
        # Plot dx_min
        axs[1].plot(self.times, self.dx_mins, 'b-', label='Minimum Distance')
        axs[1].set_ylabel('dx_min')
        axs[1].set_title('Minimum Distance vs Time')
        axs[1].grid(True)
        
        # Plot dt
        axs[2].plot(self.times, self.dts, 'g-', label='Time Step Size')
        axs[2].set_ylabel('dt')
        axs[2].set_title('Time Step Size vs Time')
        axs[2].grid(True)
        
        # Plot CFL condition (theoretical max dt)
        axs[3].plot(self.times, cfl_values, 'm-', label='CFL Condition (dx_min/v_max)')
        axs[3].plot(self.times, self.dts, 'g--', label='Actual dt', alpha=0.7)
        axs[3].set_ylabel('Time Step Size')
        axs[3].set_xlabel('Simulation Time')
        axs[3].set_title('CFL Condition vs Actual Time Step')
        axs[3].legend()
        axs[3].grid(True)
        
        # Set common x label
        fig.text(0.5, 0.04, 'Simulation Time', ha='center')
        
        # Adjust layout and save
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        plt.close()
        
        print(f"Plot saved to {output_file}")
        return True
    
    def generate_mock_data(self, num_points=100):
        """Generate mock data for testing (if no simulation output available)"""
        self.time_steps = list(range(1, num_points + 1))
        self.times = [0.01 * i for i in range(1, num_points + 1)]
        self.dts = [0.0038] * num_points  # Constant dt from example
        
        # Generate some varying v_max and dx_min
        self.v_maxs = [2000000 + 100000 * np.sin(i/10) for i in range(1, num_points + 1)]
        self.dx_mins = [400000 - 10000 * np.sin(i/5) for i in range(1, num_points + 1)]
        
        print(f"Generated {num_points} mock data points")
        return True

# Create a simple script to parse from file
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Track and plot N-Body simulation metrics")
    parser.add_argument("--input", help="Input file with simulation output")
    parser.add_argument("--output", default="nbody_metrics.png", help="Output plot file")
    parser.add_argument("--mock", action="store_true", help="Generate mock data if no input file")
    
    args = parser.parse_args()
    
    tracker = NBodyMetricsTracker()
    
    if args.input:
        tracker.extract_from_existing_output(args.input)
    elif args.mock:
        tracker.generate_mock_data()
    else:
        print("No input file provided. Use --input or --mock")
        exit(1)
    
    tracker.plot_metrics(args.output)
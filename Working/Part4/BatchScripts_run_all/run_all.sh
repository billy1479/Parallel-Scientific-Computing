#!/bin/bash
# This line is required to inform the Linux
#command line to parse the script using
#the bash shell

# Instructing SLURM to locate and assign resources
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p "cpu"
#SBATCH --qos="debug"
#SBATCH -t 00-03:00:00
#SBATCH --job-name=install_numactl
#SBATCH --output=run_all-%j.out
#SBATCH --error=run_all-%j.err

# Create data directory if it doesn't exist
mkdir -p data

# Generate initial conditions for different grid sizes
for N in 10 12 14 16; do
    echo "Generating initial conditions for N=$N..."
    python3 create_initial_conditions.py \
        --final-time 10.0 \
        --snapshots 0.1 \
        --dt 0.001 \
        --min-mass 0.1 \
        --max-mass 1.0 \
        --dim 3 \
        --N $N \
        --scenario random-grid \
        --executable-name step-0 > data/gen_N${N}.log
    
    # Move generated script to data directory and rename it
    mv step-0.sh data/input_N${N}.sh
    
    # Extract parameters from the script for direct use with step-0
    PARAMS=$(grep -oP '(?<=./step-0 ).*' data/input_N${N}.sh)
    echo "$PARAMS" > data/input_N${N}.params
done

# Run simulations for each set of initial conditions
echo "Running simulations for all generated initial conditions..."
for N in 10 12 14 16; do
    echo "Running simulation for N=$N..."
    PARAMS=$(cat data/input_N${N}.params)
    
    # Create a temporary script for this run
    cat > data/run_N${N}.sh << EOF
#!/bin/bash
./step-0 $PARAMS
EOF
    
    chmod +x data/run_N${N}.sh
    
    # Run the simulation
    data/run_N${N}.sh > data/output_N${N}.log 2>&1
    
    echo "Completed simulation for N=$N"
done

echo "All simulations completed."
#!/bin/bash
# This line is required to inform the Linux
#command line to parse the script using
#the bash shell

# Instructing SLURM to locate and assign resources
#SBATCH -N 1
#SBATCH -c 4
#SBATCH -p "cpu"
#SBATCH --qos="debug"
#SBATCH -t 00-00:30:00
#SBATCH --job-name=install_numactl
#SBATCH --output=numactl-install-%j.out
#SBATCH --error=numactl-install-%j.err

# Source the bash profile (required to use the module command)
source /etc/profile

# Set up installation directory in home
INSTALL_DIR=$HOME/local/numactl
BUILD_DIR=$HOME/build/numactl
NUMACTL_VERSION="2.0.16"  # You can change this to the latest version if needed

# Create necessary directories
mkdir -p $INSTALL_DIR
mkdir -p $BUILD_DIR

echo "Installing numactl version $NUMACTL_VERSION in $INSTALL_DIR"

# Navigate to build directory
cd $BUILD_DIR

# Download numactl source
echo "Downloading numactl source..."
wget https://github.com/numactl/numactl/releases/download/v$NUMACTL_VERSION/numactl-$NUMACTL_VERSION.tar.gz

# Extract the source
echo "Extracting source..."
tar -xzf numactl-$NUMACTL_VERSION.tar.gz
cd numactl-$NUMACTL_VERSION

# Configure, build and install numactl
echo "Configuring numactl..."
./configure --prefix=$INSTALL_DIR

echo "Building numactl..."
make

echo "Installing numactl..."
make install

# Add numactl to PATH if not already there
if ! grep -q "PATH=\$PATH:\$HOME/local/numactl/bin" $HOME/.bashrc; then
    echo 'export PATH=$PATH:$HOME/local/numactl/bin' >> $HOME/.bashrc
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/numactl/lib' >> $HOME/.bashrc
    echo "Added numactl to PATH in .bashrc"
fi

# Source .bashrc to update the current session
source $HOME/.bashrc

# Test the installation
echo "Testing numactl installation..."
$INSTALL_DIR/bin/numactl --version

echo "Getting hardware topology with numactl -H..."
$INSTALL_DIR/bin/numactl -H > numactl.out

echo "Installation complete. numactl hardware information saved to numactl.out"
echo "Please run 'source ~/.bashrc' or start a new session to use numactl"
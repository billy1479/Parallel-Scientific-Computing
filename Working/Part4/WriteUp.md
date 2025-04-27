### Timed methods

Need to run the instrumentation methods and attach to report. 

### Code performance and optimisations

- 2-D arrays used in step-0 were broken down into flat 1-D arrays as they provide better memory access patterns for GPU's. When data is offloaded to the GPU, it uses the flat arrays and it is merged back into the 2-D array when returned to the CPU. 
- OMP target data map is used to offload data from the CPU to the GPU for intense calculations like the force calculation method. This creates a persistent data region. Read-only maps was used for mass and time step size and read/write access was assigned for the flattened arrays. 
- OMP target teams distribute parallel is used to distribute the updateBody across GPU threads
- Reduction operations are used for minimum distance and maximum velcoity during the simulation
- Atomic opterations are implemented to ensure race conditions are prevented when threads are updating the same memory location
- Flat arrays are initialised in set-up then used throughout, instead of being created and destroyed in each time step like force0, force1, and force2 in step-0
- Position and velocity update in step-5-gpu were merged together into a single parallel loop and updates flat arrays on GPU to avoid transfer overhead
    - All operations are performaed on the same thread which improves cache locality
- Step-5-gpu allocates memory once during set-up instead of per time step in step-0


### 8,000 bodies (N=20)

Step-0: 
Step-5-gpu:

Graphs for both to show runtime improvement for both. Do a range from N=10 to N=20 -> use batch script for this

### Estimated speed up for N = 100,000

Pass...

### Batch scripts

Show results for each batch script

### PVM

Releasing an installation of the N-Body solver to all users on the PVM would require a centralised installation being created in the shared file system of the PVM. Initially, dependencies would need to be added to the head node via sudo (nvc++ and environment modules). Initialise a directory in the shared directory for the N-Body solver, copy the source code files to it and compile them using appropriate compiler flags. 

Following compilation of the N-Body solver, it can either be manually added to each users PATH manually which makes it harder to manage as it is in each users environment. Alternatively, an environment module could be created in the shared directory which could be loaded by users which allows for easier dependency management and cleaner user environments. This environment module will then need to be added to the system PATH.

A MPI file will then be created which allows the program to be used by more than one CPU core at a time along with distributed memory. The file can allow nodes with higher performance to allow more processes to be allocated to them, and particular processes in the simulation could be binded to certain nodes. In the case of the N-Body simulation, the head node should take the I/O processes and GPU-accelerated nodes



# Step-1 - Vectorisation
- Force update - reduced function calls by storing intermediate value in variables
- Force update - openmp simd for vectorisation with collapse(2) to break down the nested loop
- Force calculation - openmp declare simd vectorisation
- SIMD enabled version of square root, min and max function from std

# Step-2 - Parallelisation
- Brought SIMD enabled maths fucntions from step-1
- enabled parallel for schedule(dynamic) on force update loop
- enabled parallelisation for the subloop from above
- Ensuring variable visiblty via shared() and private() for loop definitions -> see lecture slides for more information in report -> variable i made private variable to do this
- All-to-one synch for the for loop updating forces for force0, force1, force2 

# Set-up loop  (parallelisation)
- Wasn't implemented because it would not bring much benefit, and there are shared variables

# Strong scaling from Claude
Let me explain strong scaling in parallel computing:

Strong scaling (also known as Speedup) refers to how the solution time varies with the number of processors for a fixed total problem size. In an ideal scenario, doubling the number of processors should halve the computation time.

The mathematical relationship is often expressed through Amdahl's Law:

Speedup(N) = T(1) / T(N)

Where:
- T(1) is the execution time on a single processor
- T(N) is the execution time on N processors
- N is the number of processors

Perfect strong scaling would mean that if you have a problem that takes 100 seconds on one processor, it should take:
- 50 seconds on 2 processors
- 25 seconds on 4 processors
- 12.5 seconds on 8 processors

However, perfect scaling is rarely achieved in practice due to several factors:

1. Communication overhead: As you add more processors, they need to communicate with each other more frequently, which adds time
2. Sequential portions: Some parts of the program cannot be parallelized (Amdahl's Law addresses this)
3. Load imbalance: Work may not be distributed equally among processors
4. Memory access patterns: Multiple processors accessing shared memory can create bottlenecks

This differs from weak scaling, where the problem size increases proportionally with the number of processors, keeping the work per processor constant.

Would you like me to elaborate on any of these aspects or discuss how to measure and optimize strong scaling?
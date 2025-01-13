# Step-1 - Vectorisation
- Force update - reduced function calls by storing intermediate value in variables
- Force update - openmp simd for vectorisation with collapse(2) to break down the nested loop
- Force calculation - openmp declare simd vectorisation
- SIMD enabled version of square root, min and max function from std

# Step-2 - Parallelisation
- Brought SIMD enabled maths fucntions from step-1
- enabled parallel for schedule(dynamic) on force update loop
- 
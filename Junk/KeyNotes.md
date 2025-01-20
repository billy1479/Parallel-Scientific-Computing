# Step-1 - Vectorisation
- Force update - reduced function calls by storing intermediate value in variables
- Force update - openmp simd for vectorisation with collapse(2) to break down the nested loop
- Force calculation - openmp declare simd vectorisation
- SIMD enabled version of square root, min and max function from std

# Step-2 - Parallelisation
- Brought SIMD enabled maths fucntions from step-1
- enabled parallel for schedule(static) on force update loop -> dynamic wasn't used as each iteration has the same computational load
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

# ISA Extensions
fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush mmx fxsr sse sse2 ht syscall nx mmxext fxsr_opt pdpe1gb rdtscp lm constant_tsc rep_good nopl nonstop_tsc cpuid extd_apicid aperfmperf pni pclmulqdq monitor ssse3   cx16 sse4_1 sse4_2 x2apic movbe popcnt aes xsave avx f16c rdrand lahf_lm cmp_legacy svm extapic cr8_legacy abm sse4a misalignsse 3dnowprefetch osvw ibs skinit wdt tce topoext perfctr_core perfctr_nb bpext perfctr_llc mwaitx cpb cat_l3 cdp_l3 hw_pstate ssbd mba ibrs ibpb stibp vmmcall fsgsbase bmi1 avx2 smep bmi2 cqm rdt_a rdseed adx smap clflushopt clwb sha_ni xsaveopt xsavec xgetbv1 xsaves cqm_llc cqm_occup_llc cqm_mbm_total cqm_mbm_local clzero irperf xsaveerptr wbnoinvd amd_ppin arat npt lbrv svm_lock nrip_save tsc_scale vmcb_clean flushbyasid decodeassists pausefilter pfthreshold avic v_vmsave_vmload vgif v_spec_ctrl umip rdpid overflow_recov succor smca sme sev sev_es
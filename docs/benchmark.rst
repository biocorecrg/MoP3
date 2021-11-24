*******************
Benchmark
*******************

We tested MoP on different datasets at the CRG's HPC where we can run up to 100 jobs in parallel (maximum 8 CPUs each) and using up to 10 GPU cards (GeForce RTX 2080 Ti).

MOP_PREPROCESS
-----------------

.. list-table:: Dataset
   
 * - 
   - Toy sample
   - Flongle
   - MinION
   - gridION
   - PromethION
 * - Fast5 files
   - 10 
   - 20 
   - 100 
   - 500 
   - 2,916 
 * - Reads
   - 40,000
   - 80,000
   - 400,000 
   - 2,000,000
   - 11,600,000
 * - Execution time
   - 10m
   - 18m
   - 1h 4m
   - 4h 19m
   - 10h 1m
 * - input folder space
   -   2 Gb
   -   5 Gb
   -  24 Gb
   - 118 Gb
   - 663 Gb
 * - work folder space
   - NOT MEASURED
   - NOT MEASURED
   - NOT MEASURED
   - 190 Gb
   - 1.1 Tb
 * - out folder space
   - 3 Gb
   - 7 Gb
   - 36 Gb
   - 178 Gb
   - 1.0 Tb
   
   
MOP_TAIL
-----------------

.. list-table:: Dataset
   
 * - 
   - Toy sample
   - Flongle
   - MinION
   - gridION
   - PromethION
 * - Fast5 files
   - 10 
   - 20 
   - 100 
   - 500 
   - 2,916 
 * - Reads
   - 40,000
   - 80,000
   - 400,000 
   - 2,000,000
   - 11,600,000
 * - Execution time
   - 20m
   - 20m
   - 47m
   - 5h 26m	
   - 2d 1h 34m
 * - work folder space
   - 153 Mb
   - 245 Mb
   - 1.2 Gb
   - 13 Gb
   - 9.2 Tb xx
 * - out folder space
   - 3 Mb 
   - 6 Mb
   - 29 Mb
   - 141 Mb
   - 1.8 Gb xx
 



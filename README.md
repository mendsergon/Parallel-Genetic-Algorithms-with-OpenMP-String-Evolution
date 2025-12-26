# **Parallel Genetic Algorithms with OpenMP: String Evolution**

### **Project Summary**

This C program implements a **parallel genetic algorithm benchmark** for solving string matching problems using **OpenMP parallelism**. The system compares three different parallelization strategies for evolutionary algorithms, measuring performance improvements across varying thread counts and migration patterns. The algorithm evolves a population of strings to match a randomly generated target sequence of DNA bases (A, C, G, T).

---

### **Core Features**

* **Three Parallel Implementation Strategies**:
  * **Version 1: Shared Population**: Single population with parallel fitness evaluation and generation processing
  * **Version 2: Independent Islands**: Multiple independent populations with final competition
  * **Version 3: Migration Islands**: Island model with periodic migration between populations

* **Complete Genetic Algorithm Framework**:
  * **Initialization**: Random population generation with DNA base characters
  * **Fitness Evaluation**: Character-by-character matching with target string
  * **Selection**: Roulette wheel (fitness-proportional) parent selection
  * **Crossover**: Single-point crossover with configurable point
  * **Mutation**: Random character replacement with adjustable rate
  * **Elitism**: Preservation of top-performing individuals

* **Comprehensive Benchmarking**:
  * Sequential baseline for comparison
  * Parallel execution time measurements
  * Generation count tracking
  * Progress reporting every 100 generations
  * Solution verification

* **Configurable Parameters**:
  * Chromosome length (default: 20)
  * Population size (default: 100)
  * Maximum generations (default: 10,000)
  * Mutation rate (default: 1%)
  * Elite count (default: 2)
  * Migration rate (default: 20%)
  * Migration interval (default: 20 generations)
  * Thread count (user-defined or system maximum)

* **Interactive Testing**: User can specify number of threads at runtime

---

### **Key Methods and Algorithms**

#### **1. Genetic Representation**
```c
typedef struct {
    char chromosome[CHROMOSOME_LENGTH + 1];
    int fitness;
} Individual;
```

#### **2. Fitness Calculation**
```c
int calculate_fitness(const char* chromosome, const char* target) {
    int fitness = 0;
    for (int i = 0; i < CHROMOSOME_LENGTH; i++) {
        if (chromosome[i] == target[i]) {
            fitness++;
        }
    }
    return fitness;
}
```

#### **3. Parallel Fitness Evaluation (Version 1)**
```c
#pragma omp parallel for
for (int i = 0; i < pop.size; i++) {
    pop.individuals[i].fitness = calculate_fitness(pop.individuals[i].chromosome, target);
}
```

#### **4. Island Model Parallelism (Version 2 & 3)**
```c
#pragma omp parallel
{
    int thread_id = omp_get_thread_num();
    Population* island = &island_pops[thread_id];
    // Each thread evolves its own population independently
    // ...
}
```

#### **5. Migration Mechanism (Version 3)**
```c
if (generation % MIGRATION_INTERVAL == 0) {
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int prev_island = (thread_id - 1 + num_threads) % num_threads;
        // Exchange top performers between islands in ring topology
        // ...
    }
}
```

---

### **File Overview**

| File | Description |
| :--- | :--- |
| **main_ga.c** | Main program with sequential baseline and parallel version execution |
| **ga_parallel.c** | Implementation of all genetic operators and three parallel strategies |
| **ga_parallel.h** | Data structures, constants, and function declarations |
| **Makefile_ga** | Build configuration with OpenMP compilation flags |

---

### **How to Compile and Run**

#### **1. Compilation**
```bash
make -f Makefile_ga
```

#### **2. Execution**
```bash
./ga_parallel
```

#### **3. Program Flow**
1. Generates random target DNA string
2. Runs sequential genetic algorithm as baseline
3. Prompts for number of threads to use
4. Executes three parallel versions with timing comparison
5. Reports results for each strategy

#### **4. Example Output**
```
Target string: ACGTGATCGTAGCTAGCTAG
Chromosome length: 20
Population size: 100
Max generations: 10000

=== SEQUENTIAL VERSION (Baseline) ===
Generation 1500: Best fitness = 20
Sequential result: 1500 generations, 0.452 seconds
Best: ACGTGATCGTAGCTAGCTAG (Fitness: 20/20)

=== PARALLEL VERSION 1: Shared Population ===
Solution found in 1400 generations!
Time: 0.127 seconds
Best individual: ACGTGATCGTAGCTAGCTAG (Fitness: 20/20)

=== PARALLEL VERSION 2: Independent Populations ===
Solution found in 1350 generations!
Time: 0.104 seconds

=== PARALLEL VERSION 3: Island Model with Migration ===
Solution found in 1200 generations!
Time: 0.089 seconds
```

---

### **Performance Characteristics**

#### **Version Comparison**
- **Version 1 (Shared)**: Best for shared-memory systems with low synchronization overhead
- **Version 2 (Islands)**: Good for maintaining diversity, less communication overhead
- **Version 3 (Migration)**: Best overall - combines diversity with information exchange

#### **Expected Speedups**
- 2-4x speedup on 4-core systems for Versions 2 & 3
- Version 1 may show super-linear speedup for large populations due to cache effects
- Migration (Version 3) typically finds solutions in fewer generations

#### **Optimal Configuration**
- Large populations benefit more from parallelism
- Migration interval of 10-50 generations works well
- 10-20% migration rate provides good balance
- Thread count should match available CPU cores

#### **Genetic Algorithm Parameters**
- **Mutation Rate**: Higher rates (0.01-0.05) help avoid local optima
- **Elite Count**: 1-5% of population preserves good solutions
- **Population Size**: Larger populations explore more but require more computation
- **Chromosome Length**: Longer strings make problem harder but show parallel benefits better

---

### **Technical Details**

#### **OpenMP Features Used**
- `#pragma omp parallel` - Creates parallel region
- `#pragma omp parallel for` - Parallel loop distribution
- `#pragma omp critical` - Thread-safe critical sections
- `#pragma omp single` - Single-threaded execution within parallel region
- `omp_get_wtime()` - High-precision timing
- `omp_get_thread_num()` - Thread identification

#### **Memory Management**
- Dynamic allocation for population arrays
- Proper freeing of memory between generations
- Stack allocation for individual chromosomes

#### **Random Number Generation**
- `srand(time(NULL))` for seeding
- Thread-safe random number usage in parallel sections
- Consistent random sequences for reproducibility

---

### **Applications and Extensions**

This implementation can be extended for:
- **Cryptography**: Password/encryption key searching
- **Bioinformatics**: DNA/protein sequence alignment
- **Machine Learning**: Hyperparameter optimization
- **Game AI**: Strategy evolution and optimization
- **Combinatorial Problems**: Traveling salesman, scheduling

To adapt for other problems:
1. Modify fitness function for new problem domain
2. Adjust genetic representation (int, float arrays instead of strings)
3. Tune mutation/crossover operators for problem constraints
4. Adjust parallel strategy based on problem characteristics

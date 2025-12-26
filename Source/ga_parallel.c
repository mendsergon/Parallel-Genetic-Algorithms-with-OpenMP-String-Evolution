#include "ga_parallel.h"

// Initialize target string with random genes (A, C, G, T)
void initialize_target(char* target) {
    const char genes[] = "ACGT";
    for (int i = 0; i < CHROMOSOME_LENGTH; i++) {
        target[i] = genes[rand() % 4];
    }
    target[CHROMOSOME_LENGTH] = '\0';
}

// Calculate fitness: number of matching characters with target
int calculate_fitness(const char* chromosome, const char* target) {
    int fitness = 0;
    for (int i = 0; i < CHROMOSOME_LENGTH; i++) {
        if (chromosome[i] == target[i]) {
            fitness++;
        }
    }
    return fitness;
}

// Initialize population with random individuals
void initialize_population(Population* pop, const char* target) {
    const char genes[] = "ACGT";
    
    for (int i = 0; i < pop->size; i++) {
        for (int j = 0; j < CHROMOSOME_LENGTH; j++) {
            pop->individuals[i].chromosome[j] = genes[rand() % 4];
        }
        pop->individuals[i].chromosome[CHROMOSOME_LENGTH] = '\0';
        pop->individuals[i].fitness = calculate_fitness(pop->individuals[i].chromosome, target);
    }
}

// Evaluate all individuals in population
void evaluate_population(Population* pop, const char* target) {
    for (int i = 0; i < pop->size; i++) {
        pop->individuals[i].fitness = calculate_fitness(pop->individuals[i].chromosome, target);
    }
}

// Select parent using roulette wheel selection
Individual select_parent(Population* pop) {
    int total_fitness = 0;
    for (int i = 0; i < pop->size; i++) {
        total_fitness += pop->individuals[i].fitness;
    }
    
    if (total_fitness == 0) {
        return pop->individuals[rand() % pop->size];
    }
    
    int random_point = rand() % total_fitness;
    int cumulative_fitness = 0;
    
    for (int i = 0; i < pop->size; i++) {
        cumulative_fitness += pop->individuals[i].fitness;
        if (cumulative_fitness >= random_point) {
            return pop->individuals[i];
        }
    }
    
    return pop->individuals[pop->size - 1];
}

// Single-point crossover
void crossover(const Individual* parent1, const Individual* parent2, Individual* child) {
    int crossover_point = rand() % CHROMOSOME_LENGTH;
    
    for (int i = 0; i < crossover_point; i++) {
        child->chromosome[i] = parent1->chromosome[i];
    }
    for (int i = crossover_point; i < CHROMOSOME_LENGTH; i++) {
        child->chromosome[i] = parent2->chromosome[i];
    }
    child->chromosome[CHROMOSOME_LENGTH] = '\0';
}

// Apply mutation with given probability
void mutate(Individual* individual) {
    const char genes[] = "ACGT";
    
    for (int i = 0; i < CHROMOSOME_LENGTH; i++) {
        if ((double)rand() / RAND_MAX < MUTATION_RATE) {
            individual->chromosome[i] = genes[rand() % 4];
        }
    }
}

// Copy individual
void copy_individual(const Individual* src, Individual* dest) {
    strcpy(dest->chromosome, src->chromosome);
    dest->fitness = src->fitness;
}

// Shuffle array for random migration
void shuffle_array(int* array, int size) {
    for (int i = size - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

// Print individual
void print_individual(const Individual* ind, int generation) {
    printf("Gen %4d: %s (Fitness: %2d/%2d)\n", 
           generation, ind->chromosome, ind->fitness, CHROMOSOME_LENGTH);
}

// VERSION 1: Shared population with parallel processing of stages
void genetic_algorithm_parallel_v1(const char* target, int num_threads) {
    printf("\n=== PARALLEL VERSION 1: Shared Population ===\n");
    printf("Target: %s\n", target);
    
    omp_set_num_threads(num_threads);
    
    // Allocate shared population
    Population pop;
    pop.size = POPULATION_SIZE;
    pop.individuals = (Individual*)malloc(pop.size * sizeof(Individual));
    
    // Initialize population
    #pragma omp parallel
    {
        #pragma omp single
        {
            initialize_population(&pop, target);
        }
    }
    
    Individual best_individual;
    int best_fitness = -1;
    int generation;
    
    double start_time = omp_get_wtime();
    
    for (generation = 0; generation < MAX_GENERATIONS; generation++) {
        // Parallel fitness evaluation
        #pragma omp parallel for
        for (int i = 0; i < pop.size; i++) {
            pop.individuals[i].fitness = calculate_fitness(pop.individuals[i].chromosome, target);
        }
        
        // Find best individual
        #pragma omp parallel
        {
            Individual local_best = pop.individuals[0];
            int local_best_fitness = local_best.fitness;
            
            #pragma omp for
            for (int i = 1; i < pop.size; i++) {
                if (pop.individuals[i].fitness > local_best_fitness) {
                    local_best = pop.individuals[i];
                    local_best_fitness = local_best.fitness;
                }
            }
            
            #pragma omp critical
            {
                if (local_best_fitness > best_fitness) {
                    best_individual = local_best;
                    best_fitness = local_best_fitness;
                }
            }
        }
        
        // Check for solution
        if (best_fitness == CHROMOSOME_LENGTH) {
            break;
        }
        
        // Create new generation in parallel
        Population new_pop;
        new_pop.size = pop.size;
        new_pop.individuals = (Individual*)malloc(pop.size * sizeof(Individual));
        
        // Elitism: keep best individuals
        #pragma omp parallel for
        for (int i = 0; i < ELITE_COUNT; i++) {
            copy_individual(&pop.individuals[i], &new_pop.individuals[i]);
        }
        
        // Generate rest of population in parallel
        #pragma omp parallel for
        for (int i = ELITE_COUNT; i < pop.size; i++) {
            Individual parent1 = select_parent(&pop);
            Individual parent2 = select_parent(&pop);
            crossover(&parent1, &parent2, &new_pop.individuals[i]);
            mutate(&new_pop.individuals[i]);
        }
        
        // Replace old population
        free(pop.individuals);
        pop = new_pop;
        
        // Print progress every 100 generations
        if (generation % 100 == 0) {
            printf("Generation %4d: Best fitness = %2d\n", generation, best_fitness);
        }
    }
    
    double end_time = omp_get_wtime();
    
    printf("\nSolution found in %d generations!\n", generation);
    printf("Time: %.3f seconds\n", end_time - start_time);
    printf("Best individual: %s (Fitness: %d/%d)\n", 
           best_individual.chromosome, best_fitness, CHROMOSOME_LENGTH);
    
    free(pop.individuals);
}

// VERSION 2: Independent populations (islands) with final competition
void genetic_algorithm_parallel_v2(const char* target, int num_threads) {
    printf("\n=== PARALLEL VERSION 2: Independent Populations ===\n");
    printf("Target: %s\n", target);
    
    omp_set_num_threads(num_threads);
    
    Population* island_pops = (Population*)malloc(num_threads * sizeof(Population));
    Individual* island_bests = (Individual*)malloc(num_threads * sizeof(Individual));
    int* island_best_fitnesses = (int*)malloc(num_threads * sizeof(int));
    
    double start_time = omp_get_wtime();
    
    // Initialize each island population
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        island_pops[thread_id].size = POPULATION_SIZE / num_threads;
        island_pops[thread_id].individuals = (Individual*)malloc(island_pops[thread_id].size * sizeof(Individual));
        
        initialize_population(&island_pops[thread_id], target);
        
        island_best_fitnesses[thread_id] = -1;
    }
    
    int global_best_fitness = -1;
    Individual global_best_individual;
    int generation;
    
    for (generation = 0; generation < MAX_GENERATIONS; generation++) {
        // Each island evolves independently
        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            Population* pop = &island_pops[thread_id];
            
            // Evaluate population
            evaluate_population(pop, target);
            
            // Find best in this island
            Individual island_best = pop->individuals[0];
            for (int i = 1; i < pop->size; i++) {
                if (pop->individuals[i].fitness > island_best.fitness) {
                    island_best = pop->individuals[i];
                }
            }
            
            island_bests[thread_id] = island_best;
            island_best_fitnesses[thread_id] = island_best.fitness;
            
            // Create new generation
            Population new_pop;
            new_pop.size = pop->size;
            new_pop.individuals = (Individual*)malloc(pop->size * sizeof(Individual));
            
            // Elitism
            for (int i = 0; i < ELITE_COUNT && i < pop->size; i++) {
                copy_individual(&pop->individuals[i], &new_pop.individuals[i]);
            }
            
            // Generate rest
            for (int i = ELITE_COUNT; i < pop->size; i++) {
                Individual parent1 = select_parent(pop);
                Individual parent2 = select_parent(pop);
                crossover(&parent1, &parent2, &new_pop.individuals[i]);
                mutate(&new_pop.individuals[i]);
            }
            
            // Replace population
            free(pop->individuals);
            *pop = new_pop;
        }
        
        // Check for global solution
        #pragma omp parallel for reduction(max:global_best_fitness)
        for (int i = 0; i < num_threads; i++) {
            if (island_best_fitnesses[i] > global_best_fitness) {
                global_best_fitness = island_best_fitnesses[i];
                #pragma omp critical
                {
                    global_best_individual = island_bests[i];
                }
            }
        }
        
        if (global_best_fitness == CHROMOSOME_LENGTH) {
            break;
        }
        
        // Print progress
        if (generation % 100 == 0) {
            printf("Generation %4d: Best fitness = %2d\n", generation, global_best_fitness);
        }
    }
    
    double end_time = omp_get_wtime();
    
    // Final competition between islands
    printf("\nFinal competition between islands:\n");
    for (int i = 0; i < num_threads; i++) {
        printf("Island %d: %s (Fitness: %2d)\n", 
               i, island_bests[i].chromosome, island_best_fitnesses[i]);
    }
    
    printf("\nSolution found in %d generations!\n", generation);
    printf("Time: %.3f seconds\n", end_time - start_time);
    printf("Best individual: %s (Fitness: %d/%d)\n", 
           global_best_individual.chromosome, global_best_fitness, CHROMOSOME_LENGTH);
    
    // Cleanup
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        free(island_pops[thread_id].individuals);
    }
    free(island_pops);
    free(island_bests);
    free(island_best_fitnesses);
}

// VERSION 3: Island model with migration
void genetic_algorithm_parallel_v3(const char* target, int num_threads) {
    printf("\n=== PARALLEL VERSION 3: Island Model with Migration ===\n");
    printf("Target: %s\n", target);
    printf("Migration: %d%% every %d generations\n", 
           (int)(MIGRATION_RATE * 100), MIGRATION_INTERVAL);
    
    omp_set_num_threads(num_threads);
    
    Population* island_pops = (Population*)malloc(num_threads * sizeof(Population));
    Individual* island_bests = (Individual*)malloc(num_threads * sizeof(Individual));
    int* island_best_fitnesses = (int*)malloc(num_threads * sizeof(int));
    
    double start_time = omp_get_wtime();
    
    // Initialize each island
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        island_pops[thread_id].size = POPULATION_SIZE / num_threads;
        island_pops[thread_id].individuals = (Individual*)malloc(island_pops[thread_id].size * sizeof(Individual));
        
        initialize_population(&island_pops[thread_id], target);
        
        island_best_fitnesses[thread_id] = -1;
    }
    
    int global_best_fitness = -1;
    Individual global_best_individual;
    int generation;
    
    for (generation = 0; generation < MAX_GENERATIONS; generation++) {
        // Each island evolves
        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            Population* pop = &island_pops[thread_id];
            
            evaluate_population(pop, target);
            
            // Find best in island
            Individual island_best = pop->individuals[0];
            for (int i = 1; i < pop->size; i++) {
                if (pop->individuals[i].fitness > island_best.fitness) {
                    island_best = pop->individuals[i];
                }
            }
            
            island_bests[thread_id] = island_best;
            island_best_fitnesses[thread_id] = island_best.fitness;
            
            // Create new generation
            Population new_pop;
            new_pop.size = pop->size;
            new_pop.individuals = (Individual*)malloc(pop->size * sizeof(Individual));
            
            // Elitism
            for (int i = 0; i < ELITE_COUNT && i < pop->size; i++) {
                copy_individual(&pop->individuals[i], &new_pop.individuals[i]);
            }
            
            for (int i = ELITE_COUNT; i < pop->size; i++) {
                Individual parent1 = select_parent(pop);
                Individual parent2 = select_parent(pop);
                crossover(&parent1, &parent2, &new_pop.individuals[i]);
                mutate(&new_pop.individuals[i]);
            }
            
            free(pop->individuals);
            *pop = new_pop;
        }
        
        // Migration phase (every MIGRATION_INTERVAL generations)
         // Fixed migration section in genetic_algorithm_parallel_v3
if (generation > 0 && generation % MIGRATION_INTERVAL == 0) {
    // Prepare migration buffers
    Individual** migration_pool = (Individual**)malloc(num_threads * sizeof(Individual*));
    int* migration_count = (int*)malloc(num_threads * sizeof(int));
    
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int migrate_n = island_pops[thread_id].size * MIGRATION_RATE;
        migration_count[thread_id] = migrate_n;
        migration_pool[thread_id] = (Individual*)malloc(migrate_n * sizeof(Individual));
        
        // Select individuals to migrate (top performers)
        for (int i = 0; i < migrate_n; i++) {
            migration_pool[thread_id][i] = island_pops[thread_id].individuals[i];
        }
    }
    
    // Perform migration in a ring topology
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int prev_island = (thread_id - 1 + num_threads) % num_threads;
        
        // Replace worst individuals with migrants from previous island
        for (int i = 0; i < migration_count[prev_island]; i++) {
            int replace_index = island_pops[thread_id].size - 1 - i;
            if (replace_index >= 0) {
                copy_individual(&migration_pool[prev_island][i], 
                              &island_pops[thread_id].individuals[replace_index]);
            }
        }
    }
    
    // Cleanup migration buffers
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        free(migration_pool[thread_id]);
    }
    free(migration_pool);
    free(migration_count);
}       
        // Check for global solution
        #pragma omp parallel for reduction(max:global_best_fitness)
        for (int i = 0; i < num_threads; i++) {
            if (island_best_fitnesses[i] > global_best_fitness) {
                global_best_fitness = island_best_fitnesses[i];
                #pragma omp critical
                {
                    global_best_individual = island_bests[i];
                }
            }
        }
        
        if (global_best_fitness == CHROMOSOME_LENGTH) {
            break;
        }
        
        // Print progress
        if (generation % 100 == 0) {
            printf("Generation %4d: Best fitness = %2d", generation, global_best_fitness);
            if (generation % MIGRATION_INTERVAL == 0) {
                printf(" [Migration occurred]");
            }
            printf("\n");
        }
    }
    
    double end_time = omp_get_wtime();
    
    printf("\nIsland results after migration:\n");
    for (int i = 0; i < num_threads; i++) {
        printf("Island %d: %s (Fitness: %2d)\n", 
               i, island_bests[i].chromosome, island_best_fitnesses[i]);
    }
    
    printf("\nSolution found in %d generations!\n", generation);
    printf("Time: %.3f seconds\n", end_time - start_time);
    printf("Best individual: %s (Fitness: %d/%d)\n", 
           global_best_individual.chromosome, global_best_fitness, CHROMOSOME_LENGTH);
    
    // Cleanup
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        free(island_pops[thread_id].individuals);
    }
    free(island_pops);
    free(island_bests);
    free(island_best_fitnesses);
}

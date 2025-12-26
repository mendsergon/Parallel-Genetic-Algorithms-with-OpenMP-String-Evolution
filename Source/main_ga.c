#include "ga_parallel.h"

int main() {
    srand(time(NULL));
    
    printf("===================================================\n");
    printf("PARALLEL GENETIC ALGORITHM - STRING GUESSING\n");
    printf("===================================================\n");
    
    // Generate random target string
    char target[CHROMOSOME_LENGTH + 1];
    initialize_target(target);
    
    printf("Target string: %s\n", target);
    printf("Chromosome length: %d\n", CHROMOSOME_LENGTH);
    printf("Population size: %d\n", POPULATION_SIZE);
    printf("Max generations: %d\n\n", MAX_GENERATIONS);
    
    // Get number of threads
    int num_threads = omp_get_max_threads();
    printf("Maximum threads available: %d\n", num_threads);
    
    // Test sequential version first (for comparison)
    printf("\n=== SEQUENTIAL VERSION (Baseline) ===\n");
    double seq_start = omp_get_wtime();
    
    Population seq_pop;
    seq_pop.size = POPULATION_SIZE;
    seq_pop.individuals = (Individual*)malloc(seq_pop.size * sizeof(Individual));
    
    initialize_population(&seq_pop, target);
    
    Individual seq_best;
    int seq_best_fitness = -1;
    int seq_generation;
    
    for (seq_generation = 0; seq_generation < MAX_GENERATIONS; seq_generation++) {
        evaluate_population(&seq_pop, target);
        
        // Find best
        for (int i = 0; i < seq_pop.size; i++) {
            if (seq_pop.individuals[i].fitness > seq_best_fitness) {
                seq_best = seq_pop.individuals[i];
                seq_best_fitness = seq_best.fitness;
            }
        }
        
        if (seq_best_fitness == CHROMOSOME_LENGTH) {
            break;
        }
        
        // Create new generation
        Population new_pop;
        new_pop.size = seq_pop.size;
        new_pop.individuals = (Individual*)malloc(seq_pop.size * sizeof(Individual));
        
        // Elitism
        for (int i = 0; i < ELITE_COUNT; i++) {
            copy_individual(&seq_pop.individuals[i], &new_pop.individuals[i]);
        }
        
        for (int i = ELITE_COUNT; i < seq_pop.size; i++) {
            Individual parent1 = select_parent(&seq_pop);
            Individual parent2 = select_parent(&seq_pop);
            crossover(&parent1, &parent2, &new_pop.individuals[i]);
            mutate(&new_pop.individuals[i]);
        }
        
        free(seq_pop.individuals);
        seq_pop = new_pop;
        
        if (seq_generation % 500 == 0) {
            printf("Generation %4d: Best fitness = %2d\n", seq_generation, seq_best_fitness);
        }
    }
    
    double seq_end = omp_get_wtime();
    printf("\nSequential result: %d generations, %.3f seconds\n", 
           seq_generation, seq_end - seq_start);
    printf("Best: %s (Fitness: %d/%d)\n", 
           seq_best.chromosome, seq_best_fitness, CHROMOSOME_LENGTH);
    
    free(seq_pop.individuals);
    
    // Run parallel versions
    printf("\n===================================================\n");
    printf("PARALLEL VERSIONS COMPARISON\n");
    printf("===================================================\n");
    
    printf("\nEnter number of threads to use (default: %d): ", num_threads);
    int user_threads;
    if (scanf("%d", &user_threads) == 1 && user_threads > 0) {
        num_threads = user_threads;
    }
    
    // Version 1
    genetic_algorithm_parallel_v1(target, num_threads);
    
    // Version 2  
    genetic_algorithm_parallel_v2(target, num_threads);
    
    // Version 3
    genetic_algorithm_parallel_v3(target, num_threads);
    
    // Performance comparison
    printf("\n===================================================\n");
    printf("PROGRAM FINISHED\n");
    printf("===================================================\n");
 
    
    return 0;
}

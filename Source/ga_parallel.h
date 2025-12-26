#ifndef GA_PARALLEL_H
#define GA_PARALLEL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <math.h>

// Configuration
#define MAX_GENERATIONS 10000
#define POPULATION_SIZE 100
#define CHROMOSOME_LENGTH 20
#define MUTATION_RATE 0.01
#define ELITE_COUNT 2
#define MIGRATION_RATE 0.2     // 20% migration
#define MIGRATION_INTERVAL 20  // Every 20 generations

// Genetic algorithm types
typedef struct {
    char chromosome[CHROMOSOME_LENGTH + 1];  // +1 for null terminator
    int fitness;
} Individual;

typedef struct {
    Individual* individuals;
    int size;
} Population;

// Function declarations
void initialize_target(char* target);
int calculate_fitness(const char* chromosome, const char* target);
void initialize_population(Population* pop, const char* target);
void evaluate_population(Population* pop, const char* target);
Individual select_parent(Population* pop);
void crossover(const Individual* parent1, const Individual* parent2, Individual* child);
void mutate(Individual* individual);
void copy_individual(const Individual* src, Individual* dest);

// Parallel versions
void genetic_algorithm_parallel_v1(const char* target, int num_threads);
void genetic_algorithm_parallel_v2(const char* target, int num_threads);
void genetic_algorithm_parallel_v3(const char* target, int num_threads);

// Helper functions
void print_individual(const Individual* ind, int generation);
void shuffle_array(int* array, int size);

#endif

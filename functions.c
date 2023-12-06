// Include everything necessary here
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "functions.h"

#define ELITISM_RATE 0.1

//?  my own stuff for internal usage

// function to copy an array from one destination to another
void copy_array(int NUM_VARIABLES, double data[NUM_VARIABLES], double target[NUM_VARIABLES])
{
    for (int i = 0; i < NUM_VARIABLES; i++)
    {
        target[i] = data[i];
    }
}

// binary search for when we find the parents through cumulative probabilities
// this just finds the index of the first element that is greater than or equal to the value
// binary search will keep halving the search space until it finds the value
int binary_search(double* array, int size, double value)
{
    int left = 0;
    int right = size - 1;
    while (left < right)
    {
        int mid = left + (right - left) / 2;
        if (array[mid] < value)
        {
            left = mid + 1;
        }
        else
        {
            right = mid;
        }
    }
    return (array[left] >= value) ? left : -1;
}

// compare function for qsort
int compare(const void* a, const void* b)
{
    double* x = (double*)a;
    double* y = (double*)b;
    return *x < *y;
}

// sorting function for the population
// we need to sort the fitness and the population in place so we can pick
// the best individuals for elitism
void sort_population(int POPULATION_SIZE, int NUM_VARIABLES, double fitness[POPULATION_SIZE], double population[POPULATION_SIZE][NUM_VARIABLES])
{
    // sort fitness
    qsort(fitness, POPULATION_SIZE, sizeof(double), compare);

    // sort population in place
    for (int i = 0; i < POPULATION_SIZE; ++i)
    {
        for (int j = 0; j < POPULATION_SIZE; ++j)
        {
            if (fitness[i] == 1 / (Objective_function(NUM_VARIABLES, population[j]) + 1e-6))
            {
                double temp[NUM_VARIABLES];
                copy_array(NUM_VARIABLES, population[i], temp);
                copy_array(NUM_VARIABLES, population[j], population[i]);
                copy_array(NUM_VARIABLES, temp, population[j]);
            }
        }
    }
}

//? everything else below is given

// generate random float between min and max because
// C doesn't have a built in function for this
double generate_random(double min, double max)
{
    double range = max - min;
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

// i literally never use this function
unsigned int generate_int()
{
    return (unsigned int)rand();
}

// Function to initialize a random population
void generate_population(int POPULATION_SIZE, int NUM_VARIABLES, double population[POPULATION_SIZE][NUM_VARIABLES], double Lbound[NUM_VARIABLES], double Ubound[NUM_VARIABLES])
{
    // just iterate through the population and generate a random number for each variable
    for (int i = 0; i < POPULATION_SIZE; ++i)
    {
        for (int j = 0; j < NUM_VARIABLES; ++j)
        {
            population[i][j] = generate_random(Lbound[j], Ubound[j]);
        }
    }
}

// Function to compute the objective function for each member of the population
void compute_objective_function(int POPULATION_SIZE, int NUM_VARIABLES, double population[POPULATION_SIZE][NUM_VARIABLES], double fitness[POPULATION_SIZE])
{

    /* compute "fitness[i]"" for each set of decision variables (individual) or each row in "population"
    by calling "Objective_function" */

    // we compute as 1 / (objective function + 1e-6) because we want
    // to get the global minimum, of the ackley function, which is 0
    for (int i = 0; i < POPULATION_SIZE; ++i)
    {
        fitness[i] = 1 / (Objective_function(NUM_VARIABLES, population[i]) + 1e-6);
    }
}

// main crossover function
void crossover(int POPULATION_SIZE, int NUM_VARIABLES, double fitness[POPULATION_SIZE], double new_population[POPULATION_SIZE][NUM_VARIABLES], double population[POPULATION_SIZE][NUM_VARIABLES], double crossover_rate)
{
    /* Implement the logic of crossover function here based on "fitness_probs" or each set
    of decision variables (individual) or each row in "population".
    And save the new population in "new_population"*/

    // compute the fitness probabilities and cumulative probabilities
    // beforehand so we don't have to do it for each pair of parents
    double fitness_probs[POPULATION_SIZE];
    double fitness_sum = 0.0;
    double cumulative_probs[POPULATION_SIZE];
    for (int i = 0; i < POPULATION_SIZE; ++i)
    {
        fitness_sum += fitness[i];
    }
    for (int i = 0; i < POPULATION_SIZE; ++i)
    {
        fitness_probs[i] = fitness[i] / fitness_sum;
    }
    cumulative_probs[0] = fitness_probs[0];
    for (int i = 1; i < POPULATION_SIZE; ++i)
    {
        cumulative_probs[i] = cumulative_probs[i - 1] + fitness_probs[i];
    }

    // initialize the parents and children
    double parent1[NUM_VARIABLES];
    double parent2[NUM_VARIABLES];
    double child1[NUM_VARIABLES];
    double child2[NUM_VARIABLES];

#ifdef ELITISM
    // sort the population
    sort_population(POPULATION_SIZE, NUM_VARIABLES, fitness, population);

    // copy the elites to the new population
    // (parents at the start of the population are the elites)
    // we want to keep the number of elites even
    int num_elites = (int)round(POPULATION_SIZE * ELITISM_RATE);
    num_elites = (num_elites % 2 == 0) ? num_elites : num_elites - 1;

    for (int i = 0; i < num_elites; ++i)
    {
        copy_array(NUM_VARIABLES, population[i], new_population[i]);
    }

    // start at the first non-elite
    int start = num_elites;
#else
    // start at the first individual (we do this because of the elitism stuff above)
    int start = 0;
#endif

    // select and crossover
    for (int i = start; i < POPULATION_SIZE; i += 2)
    {
        // select parents

        // we select parents based on a random number and the cumulative probabilities
        // the index is found using binary search

        // parent 1
        double r = generate_random(0, 1);
        int parent1_index = binary_search(cumulative_probs, POPULATION_SIZE, r);
        copy_array(NUM_VARIABLES, population[parent1_index], parent1);

        // parent 2
        r = generate_random(0, 1);
        int parent2_index = binary_search(cumulative_probs, POPULATION_SIZE, r);
        copy_array(NUM_VARIABLES, population[parent2_index], parent2);

        // crossover
        r = generate_random(0, 1);
        int crossover_point = generate_int() % NUM_VARIABLES;
        if (r < crossover_rate) // we only crossover if the random number is less than the crossover rate
        {
            // we copy the first part of the parent to the child
            for (int j = 0; j < crossover_point; ++j)
            {
                child1[j] = parent1[j];
                child2[j] = parent2[j];
            }

            // we copy the second part of the parent to the child
            for (int j = crossover_point; j < NUM_VARIABLES; ++j)
            {
                child1[j] = parent2[j];
                child2[j] = parent1[j];
            }
        }
        else
        {
            // otherwise we just copy the parents to the children
            copy_array(NUM_VARIABLES, parent1, child1);
            copy_array(NUM_VARIABLES, parent2, child2);
        }

        // save the new population
        copy_array(NUM_VARIABLES, child1, new_population[i]);
        copy_array(NUM_VARIABLES, child2, new_population[i + 1]);
    }
}

void mutate(int POPULATION_SIZE, int NUM_VARIABLES, double new_population[POPULATION_SIZE][NUM_VARIABLES], double population[POPULATION_SIZE][NUM_VARIABLES], double Lbound[NUM_VARIABLES], double Ubound[NUM_VARIABLES], double mutate_rate)
{
    // calculate the number of mutations
    int num_mutations = (int)round(POPULATION_SIZE * NUM_VARIABLES * mutate_rate);

    // and we iterate through the population and mutate that number of times
    for (int i = 0; i < num_mutations; ++i)
    {
        int row = generate_int() % POPULATION_SIZE;
        int col = generate_int() % NUM_VARIABLES;
        new_population[row][col] = generate_random(Lbound[col], Ubound[col]);
    }

    // copy the new population to the old one
    for (int i = 0; i < POPULATION_SIZE; ++i)
    {
        copy_array(NUM_VARIABLES, new_population[i], population[i]);
    }
}
// Include everything necessary here
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "functions.h"

// since we cannot change functions.h, ELITISM_RATE must be changed
// in both GA.c and functions.c
#define ELITISM_RATE 0.1

int main(int argc, char *argv[])
{
    // <YOUR CODE: Handle the possible errors in input data given by the user and say how to execute the code>

    // check if the number of arguments is correct
    if (argc != 6)
    {
        printf("Usage: %s <POPULATION_SIZE> <MAX_GENERATIONS> <crossover_rate> <mutate_rate> <stop_criteria>\n", argv[0]);
        printf("Example: %s 100 1000 0.8 0.1 0.001\n", argv[0]);
        exit(1);
    }

    // <YOUR CODE: Assign all inputs given by the user argv[i]> like
    // POPULATION_SIZE, MAX_GENERATIONS, crossover_rate, mutate_rate, stop_criteria

    int POPULATION_SIZE = atoi(argv[1]);
    int MAX_GENERATIONS = atoi(argv[2]);
    double crossover_rate = atof(argv[3]);
    double mutate_rate = atof(argv[4]);
    double stop_criteria = atof(argv[5]);

    // check if the inputs are valid
    if (POPULATION_SIZE <= 0 || MAX_GENERATIONS <= 0 || crossover_rate < 0 || crossover_rate > 1 || mutate_rate < 0 || mutate_rate > 1 || stop_criteria < 0)
    {
        printf("Invalid input\n");
        exit(1);
    }

    // seed the random number generator
    srand(time(NULL));

    // ###################################################################################
    // you dont need to change anything here
    // the number of variables
    int NUM_VARIABLES = 2;
    // the lower bounds of variables
    double Lbound[] = {-5.0, -5.0};
    // the upper bounds of variable
    double Ubound[] = {5.0, 5.0};
    // ###################################################################################

    // <YOUR CODE: Here make all the initial print outs>
    printf("Genetic algorithm start parameters\n");
    printf("-------------------------\n");
    printf("- Population size: %d\n", POPULATION_SIZE);
    printf("- Max generations: %d\n", MAX_GENERATIONS);
    printf("- Crossover rate: %f\n", crossover_rate);
    printf("- Mutation rate: %f\n", mutate_rate);
    printf("- Stop criteria: %.16f\n\n", stop_criteria);
    printf("- Number of variables: %d\n", NUM_VARIABLES);


    printf("- Lower bounds: [");
    for (int i = 0; i < NUM_VARIABLES; ++i)
    {
        printf("%f", Lbound[i]);
        if (i != NUM_VARIABLES - 1)
        {
            printf(", ");
        }
    }
    printf("]\n");
    printf("- Upper bounds: [");
    for (int i = 0; i < NUM_VARIABLES; ++i)
    {
        printf("%f", Ubound[i]);
        if (i != NUM_VARIABLES - 1)
        {
            printf(", ");
        }
    }
    printf("]\n\n");
#ifdef ELITISM
    printf("- Elitism: ON\n");
    printf("- Elitism rate: %f\n\n\n", ELITISM_RATE);
#else
    printf("- Elitism: OFF\n\n");
#endif

    // wall clock time instead of cpu time
    // i was doing testing with multithreading and cpu time was not accurate
    struct timeval start, end;
    double wall_time_used;
    gettimeofday(&start, NULL);

    // <YOUR CODE: Declare all the arrays you need here>
    double population[POPULATION_SIZE][NUM_VARIABLES];
    double new_population[POPULATION_SIZE][NUM_VARIABLES];
    double fitness[POPULATION_SIZE];

    // these variables are all for calculating the best solution
    double best_solution[NUM_VARIABLES];
    double best_fitness = 0.0;
    double prev_best_fitness = 0.0;
    int generations_taken = 0;
    int gens_low_fitness = 0;

    // <YOUR CODE: Call generate_population function to initialize the "population"> like:
    // generate_population(POPULATION_SIZE, NUM_VARIABLES, population, Lbound, Ubound);
    generate_population(POPULATION_SIZE, NUM_VARIABLES, population, Lbound, Ubound);

    // iteration starts here. The loop continues until MAX_GENERATIONS is reached
    // Or stopping criteria is met
    for (int generation = 0; generation < MAX_GENERATIONS; generation++)
    {
        // <YOUR CODE: Compute the fitness values using objective function for
        // each row in "population" (each set of variables)> like:
        // compute_objective_function(POPULATION_SIZE, NUM_VARIABLES, population, fitness);
        compute_objective_function(POPULATION_SIZE, NUM_VARIABLES, population, fitness);

        // <YOUR CODE: Here implement the logic of finding best solution with minimum fitness value
        // and the stopping criteria>

        // find the best solution in this generation
        best_fitness = fitness[0];
        int best_solution_index = 0;
        for (int i = 0; i < POPULATION_SIZE; ++i)
        {
            if (fitness[i] > best_fitness)
            {
                best_fitness = fitness[i];
                best_solution_index = i;
            }
        }

        // check if we have reached the stopping criteria for 5 consecutive generations
        // if so, break out of the loop
        if (gens_low_fitness > (int)round(MAX_GENERATIONS * 0.05))
        {
            break;
        }

        if (fabs(best_fitness - prev_best_fitness) < stop_criteria)
        {
            gens_low_fitness++;
        }
        else
        {
            // we only want to reset counter and update best solution if 
            // the best fitness has improved
            if (best_fitness > prev_best_fitness)
            {
                // update the best fitness
                prev_best_fitness = best_fitness;
                
                // update the best solution
                for (int i = 0; i < NUM_VARIABLES; ++i)
                {
                    best_solution[i] = population[best_solution_index][i];
                }

                // reset the number of generations with low fitness
                gens_low_fitness = 0;
            }
        }

        // <YOUR CODE: Here call the crossover function>
        crossover(POPULATION_SIZE, NUM_VARIABLES, fitness, new_population, population, crossover_rate);

        // <YOUR CODE: Here call the mutation function>
        mutate(POPULATION_SIZE, NUM_VARIABLES, new_population, population, Lbound, Ubound, mutate_rate);

        // Now you have the a new population, and it goes to the beginning of loop to re-compute all again
        ++generations_taken;
    }

    // <YOUR CODE: Jump to this part of code if the stopping criteria is met before MAX_GENERATIONS is met>

    // ###################################################################################
    // You dont need to change anything here
    // Here we print the CPU time taken for your code
    gettimeofday(&end, NULL);
    wall_time_used = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    // ###################################################################################

    // <Here print out the best solution and objective function value for the best solution like the format>
    printf("Results\n");
    printf("-------------------------\n");
    printf("- Best solution: [");
    for (int i = 0; i < NUM_VARIABLES; ++i)
    {
        printf("%.16f", best_solution[i]);
        if (i != NUM_VARIABLES - 1)
        {
            printf(", ");
        }
    }
    printf("]\n");
    printf("- Best fitness: %f\n", best_fitness);
    printf("- Generations taken: %d\n", generations_taken);
    printf("- Time used: %.16f s\n", wall_time_used);

    return 0;
}

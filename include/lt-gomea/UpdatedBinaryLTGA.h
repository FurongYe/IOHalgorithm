#ifndef _LTGA_H_
#define _LTGA_H_

#include <stdint.h>

/*-=-=-=-=-=-=-=-=-=-=-=- Section Structs -=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-*/

const double FITNESS_EPSILON = 0.0000000001;
typedef unsigned int uint;

struct InputParametersC {
  uint M, k, o, b;
};


struct LTGAParameters {
  short     write_generational_statistics, 
            write_generational_solutions, 
            print_verbose_overview,
            print_FOSs_contents,
            use_ilse,
            use_random_linkage,
            use_value_to_reach,
            use_fitness_variance_tolerance,
            use_premature_stopping;
  long      maximum_number_of_evaluations,
            maximum_number_of_milliseconds;
  uint      m, k, o, b;
  double    global_optima_score,
            fitness_variance_tolerance;
};

struct LTGAResultC {
  char   *elitist_solution;
  double  elitist_fitness;
  long    elitist_evaluations;
  long    elitist_time;
  short   global_optima_found;
  int     pop_size_glob_opt; 
};

/**
 * All global variables in a struct
 */
struct Environment {
  char       *terminated;                                             /* Whether a specific GOMEA with the restart scheme has terminated. */
  char       *elitist_solution,                                       /* The very best solution ever evaluated. */
            **sostr,                                                  /* Set of solutions to reach. */
           ***populations,                                            /* The populations containing the solutions. */
           ***offsprings;                                             /* Offspring solutions (one set per population). */
  short       halt_execution,                                         /* Whether the execution of the library should be halted (due to the desired value having been reached or the max number of evaluations has been exceeded) */
              write_generational_statistics,                          /* Whether to compute and write statistics every generation (0 = no). */
              write_generational_solutions,                           /* Whether to write the population every generation (0 = no). */
              print_verbose_overview,                                 /* Whether to print a overview of settings (0 = no). */
              print_FOSs_contents,                                    /* Whether to print the contents of the FOS structure each generation (0 = no). */
              use_ilse,                                               /* Whether incremental linkage subset evaluation is to be used. */
              use_random_linkage,                                     /* Whether to use random linkage during FOS model construction */
              use_fitness_variance_tolerance,                         /* Whether to use fitness variance to terminate a population */
              use_premature_stopping,                                 /* Whether to use premature stopping; when a population has worse avg fitness than a bigger population, stop the smaller population */
              vosostr_hit_status,                                     /* Whether the vosostr hit has happened yet: a solution has been evaluated and a value >= VTR or a STR has been found (0 = no, 1 = yes, 2 = yes, but this is no longer the first time. */
              vtr_exists,                                             /* Whether a vtr exists. */
              sostr_exists;                                           /* Whether a sostr exists. */
  int         problem_index,                                          /* The index of the optimization problem. */
              FOSs_structure_index,                                   /* The index of the FOS structure. */
              number_of_parameters,                                   /* The number of parameters to be optimized. */
              number_of_solutions_in_sostr,                           /* The number of solutions in the set of solutions to reach. */
              number_of_generations,                                  /* The current generation count. */
             *population_sizes,                                       /* The number of solutions in each population. */
              base_population_size,                                   /* The minimum population size used in the smallest GOMEA instance. */
             *number_of_subgenerations_per_GOMEA,                     /* The number of subgenerations per GOMEA instance. */
              number_of_subgenerations_per_GOMEA_factor,              /* The factor by which the number of subgenerations increases with every new population. */
             *no_improvement_stretchs,                                /* The number of subsequent generations without an improvement for every GOMEA. */
              number_of_GOMEAs,                                       /* The number of GOMEAs currently running in multipop configuration. */
              maximum_number_of_GOMEAs,                               /* The maximum number of GOMEAs running in multipop configuration. */
           ***FOSs,                                                   /* The family of subsets linkage struture. */
            **FOSs_number_of_indices,                                 /* The number of variables in each linkage subset. */
             *FOSs_length,                                            /* The number of linkage subsets. */
              minimum_GOMEA_index;                                    /* The minimum GOMEA index that corresponds to the GOMEA that is still allowed to run (lower ones should be stopped because of average fitness being lower than that of a higher one). */
  long        maximum_number_of_evaluations,                          /* The maximum number of evaluations. */
              maximum_number_of_milliseconds,                         /* The maximum number of milliseconds. */
              timestamp_start,                                        /* The time stamp in milliseconds for when the program was started. */
              timestamp_start_after_init,                             /* The time stamp in milliseconds for when the algorithm was started (after problem initialization). */
              number_of_evaluations,                                  /* The current number of times a function evaluation was performed. */
              elitist_solution_number_of_evaluations,                 /* The number of evaluations until the elitist solution. */
              elitist_solution_hitting_time,                          /* The hitting time of the elitist solution. */
              vosostr_number_of_evaluations,                          /* The number of evaluations until a solution that was to be reached (either vtr or in sostr). */
              vosostr_hitting_time;                                   /* The hitting time of a solution that was to be reached (either vtr or in sostr). */
  long long   number_of_bit_flip_evaluations,                         /* The number of bit-flip evaluations. */
              elitist_solution_number_of_bit_flip_evaluations,        /* The number of bit-flip evaluations until the elitist solution. */
              vosostr_number_of_bit_flip_evaluations;                 /* The number of bit-flip evaluations until a solution that was to be reached (either vtr or in sostr). */
  double      elitist_solution_objective_value,                       /* The objective value of the elitist solution. */
              elitist_solution_constraint_value,                      /* The constraint value of the elitist solution. */
              vtr,                                                    /* The value to reach (fitness of best solution that is feasible). */
            **objective_values,                                       /* Objective values for population members. */
            **constraint_values,                                      /* Sum of all constraint violations for population members. */
            **objective_values_offsprings,                            /* Objective values of selected solutions. */
            **constraint_values_offsprings,                           /* Sum of all constraint violations of selected solutions. */
             *objective_values_best_of_generation,
             *constraint_values_best_of_generation,
             *average_objective_values,
             *average_constraint_values,
           ***dependency_matrices;                                    /* Measure of dependency between any two variables (higher is more dependent). */
  int64_t     random_seed,                                            /* The seed used for the random-number generator. */
              random_seed_changing;                                   /* Internally used variable for randomly setting a random seed. */

  // int         njobs = 0,                                              /* Number of jobs. */
  //             nmach = 0,                                              /* Number of machines. */
  //           **ptime;                                                  /* Processing time of each jobs on each machine */

  struct InputParametersC  input_parameters;                           /* The Mk Landscape's input parameters passed using a file, CLI, or ABI. */
  uint                  **cliques;                                    /* The Mk Landscape's generated clique tree's cliques. Stores per clique the problem indices */
  double                **codomain_values;                            /* The Mk Landscape's generated or read codomain values. Stores the subfunction per clique */
  double                  global_optima_score;                        /* The Mk Landscape's global optima score/fitness, defined by the clique tree and the codomain values */
  double                  fitness_variance_tolerance;                 /* Fitness variance tolerance for the populations; below it, a population is terminated. */

  struct LTGAResultC       ltga_result;                                /*The LTGA results */
};
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
/* Global variables */
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*-=-=-=-=-=-=-=-=-=-=-=-= Section Header Functions -=-=-=-=-=-=-=-=-=-=-=-=*/
void *Malloc( long size );
double randomRealUniform01( struct Environment *env );
int randomInt( struct Environment *env, int maximum );
int *randomPermutation( struct Environment *env, int n );

void problemEvaluation( struct Environment *env, int index, char *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, char *parameters_before, double objective_value_before, double constraint_value_before, int GOMEA_index );
void installedProblemEvaluation( int index, double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *parameters_before, double objective_value_before, double constraint_value_before );

void initializeNewGOMEA( struct Environment *env );
void initializeNewGOMEAMemory( struct Environment *env );
void initializeNewGOMEAPopulationAndFitnessValues( struct Environment *env );
void initializeValueAndSetOfSolutionsToReach( struct Environment *env );
short initializeValueToReach( struct Environment *env );
short initializeSetOfSolutionsToReach( struct Environment *env );
void initializeRandomNumberGenerator( struct Environment *env );

void selectFinalSurvivorsSpecificGOMEA( struct Environment *env, int gomea_index );
char betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );
char equalFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y );

void writeGenerationalStatistics( struct Environment *env );
void writeGenerationalSolutions( struct Environment *env, char is_final_generation );
void writeRunningTime( struct Environment *env, char *filename );
void writeElitistSolution( struct Environment *env );

char checkTermination( struct Environment *env );
char checkNumberOfEvaluationsTerminationCondition( struct Environment *env );
char checkVOSOSTRTerminationCondition( struct Environment *env );
char checkNumberOfMilliSecondsTerminationCondition( struct Environment *env );
char checkFitnessVarianceTermination(struct Environment *env, int GOMEA_index );

void generationalStepAllGOMEAs( struct Environment *env );
void makeOffspringSpecificGOMEA( struct Environment *env, int gomea_index );
void learnFOSSpecificGOMEA( struct Environment *env, int gomea_index );
void selectForLearningLT(int gomea_index);
int **learnLTFOSSpecificGOMEA( struct Environment *env, int gomea_index, short compute_dependency_matrices, short compute_parent_child_relations, int *number_of_parent_child_relations );
int determineNearestNeighbour( int index, double **S_matrix, int *mpm_number_of_indices, int mpm_length );
void computeDependencyMatrixSpecificGOMEA( int gomea_index );
void computeDependencyMatrixMutualInformationSpecificGOMEA( struct Environment *env, int gomea_index );
double *estimateParametersForSingleBinaryMarginal( struct Environment *env, int gomea_index, int *indices, int number_of_indices, int *factor_size );
void uniquifyFOSSpecificGOMEA( int gomea_index );
void printFOSContentsSpecificGOMEA( struct Environment *env, int gomea_index );
double log2( double x );
void generateAndEvaluateNewSolutionsToFillOffspringSpecificGOMEA( struct Environment *env, int gomea_index );
char *generateAndEvaluateNewSolutionBinarySpecificGOMEA( struct Environment *env, int gomea_index, int parent_index, double *obj, double *con );
void shuffleFOSSpecificGOMEA( struct Environment *env, int gomea_index );
void shuffleFOSSubsetsSpecificGOMEA( struct Environment *env, int gomea_index );

void ezilaitiniAllGOMEAs( struct Environment *env );
void ezilaitiniSpecificGOMEA( struct Environment *env, int gomea_index );
void ezilaitiniSpecificGOMEAMemoryForPopulationAndOffspring( struct Environment *env, int gomea_index );
void ezilaitiniValueAndSetOfSolutionsToReach( struct Environment *env );
void ezilaitiniProblem( int index );

long getMilliSecondsRunning( struct Environment *env );
long getMilliSecondsRunningAfterInit( struct Environment *env );
long getMilliSecondsRunningSinceTimeStamp( long timestamp );
long getCurrentTimeStampInMilliSeconds();

void run( struct Environment *env );
void multiPopGOMEA( struct Environment *env );
int main( int argc, char **argv );

void calculateFitness( struct Environment *env, char *parameters, double *objective_value, double *constraint_value);
int setParameters(struct Environment *env, struct LTGAParameters *run_parameters, uint **cliques_memory, double **codomain_memory );
extern struct LTGAResultC run_with_parameters(struct LTGAParameters run_parameters, uint **cliques_memory, double **codomain_memory );
struct Environment getNewEnvironment();

#endif
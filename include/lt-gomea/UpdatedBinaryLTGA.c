/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Header -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * PermutationGOMEA.c
 *
 * Copyright (c) Peter A.N. Bosman
 *
 * The software in this file is the proprietary information of
 * Peter A.N. Bosman.
 *
 * IN NO EVENT WILL THE AUTHOR OF THIS SOFTWARE BE LIABLE TO YOU FOR ANY
 * DAMAGES, INCLUDING BUT NOT LIMITED TO LOST PROFITS, LOST SAVINGS, OR OTHER
 * INCIDENTIAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR THE INABILITY
 * TO USE SUCH PROGRAM, EVEN IF THE AUTHOR HAS BEEN ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGES, OR FOR ANY CLAIM BY ANY OTHER PARTY. THE AUTHOR MAKES NO
 * REPRESENTATIONS OR WARRANTIES ABOUT THE SUITABILITY OF THE SOFTWARE, EITHER
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE
 * AUTHOR SHALL NOT BE LIABLE FOR ANY DAMAGES SUFFERED BY ANYONE AS A RESULT OF
 * USING, MODIFYING OR DISTRIBUTING THIS SOFTWARE OR ITS DERIVATIVES.
 */
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=- Section Global Variables -=-=-=-=-=-=-=-=-=-=-=-=-*/
#if defined(_WIN32) || defined(WIN32) || defined(__CYGWIN__) || defined(__MINGW32__) || defined(__BORLANDC__)
#define OS_WIN
#endif
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Includes -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef OS_WIN
#include <stdint.h>
#endif
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include "UpdatedBinaryLTGA.h"
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


#define MAX(a,b) \
  ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Constants -=-=-=-=-=-=-=-=-=-=-=-=-=-*/
#define FALSE 0
#define TRUE 1
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*-=-=-=-=-=-=-=-=-=-=-= Section Elementary Operations -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Allocates memory and exits the program in case of a memory allocation failure.
 */
void *Malloc( long size )
{
  void *result;

  result = (void *) malloc( size );
  if( !result )
  {
    printf( "\n" );
    printf( "Error while allocating memory in Malloc( %ld ), aborting program.", size );
    printf( "\n" );

    exit( 0 );
  }

  return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Merge Sort -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Sorts an array of integers and returns the sort-order (large to small).
 */
int *mergeSortIntegersDecreasing( int *array, int array_size )
{
  int i, *sorted, *tosort;

  sorted = (int *) Malloc( array_size * sizeof( int ) );
  tosort = (int *) Malloc( array_size * sizeof( int ) );
  for( i = 0; i < array_size; i++ )
    tosort[i] = i;

  if( array_size == 1 )
    sorted[0] = 0;
  else
    mergeSortIntegersDecreasingWithinBounds( array, sorted, tosort, 0, array_size-1 );

  free( tosort );

  return( sorted );
}

/**
 * Subroutine of merge sort, sorts the part of the array between p and q.
 */
void mergeSortIntegersDecreasingWithinBounds( int *array, int *sorted, int *tosort, int p, int q )
{
  int r;

  if( p < q )
  {
    r = (p + q) / 2;
    mergeSortIntegersDecreasingWithinBounds( array, sorted, tosort, p, r );
    mergeSortIntegersDecreasingWithinBounds( array, sorted, tosort, r+1, q );
    mergeSortIntegersDecreasingMerge( array, sorted, tosort, p, r+1, q );
  }
}

/**
 * Subroutine of merge sort, merges the results of two sorted parts.
 */
void mergeSortIntegersDecreasingMerge( int *array, int *sorted, int *tosort, int p, int r, int q )
{
  int i, j, k, first;

  i = p;
  j = r;
  for( k = p; k <= q; k++ )
  {
    first = 0;
    if( j <= q )
    {
      if( i < r )
      {
        if( array[tosort[i]] > array[tosort[j]] )
          first = 1;
      }
    }
    else
      first = 1;

    if( first )
    {
      sorted[k] = tosort[i];
      i++;
    }
    else
    {
      sorted[k] = tosort[j];
      j++;
    }
  }

  for( k = p; k <= q; k++ )
    tosort[k] = sorted[k];
}

/**
 * Sorts an array of doubles and returns the sort-order (large to small).
 */
int *mergeSortDoublesIncreasing( double *array, int array_size )
{
  int i, *sorted, *tosort;

  sorted = (int *) Malloc( array_size * sizeof( int ) );
  tosort = (int *) Malloc( array_size * sizeof( int ) );
  for( i = 0; i < array_size; i++ )
    tosort[i] = i;

  if( array_size == 1 )
    sorted[0] = 0;
  else
    mergeSortDoublesIncreasingWithinBounds( array, sorted, tosort, 0, array_size-1 );

  free( tosort );

  return( sorted );
}

/**
 * Subroutine of merge sort, sorts the part of the array between p and q.
 */
void mergeSortDoublesIncreasingWithinBounds( double *array, int *sorted, int *tosort, int p, int q )
{
  int r;

  if( p < q )
  {
    r = (p + q) / 2;
    mergeSortDoublesIncreasingWithinBounds( array, sorted, tosort, p, r );
    mergeSortDoublesIncreasingWithinBounds( array, sorted, tosort, r+1, q );
    mergeSortDoublesIncreasingMerge( array, sorted, tosort, p, r+1, q );
  }
}

/**
 * Subroutine of merge sort, merges the results of two sorted parts.
 */
void mergeSortDoublesIncreasingMerge( double *array, int *sorted, int *tosort, int p, int r, int q )
{
  int i, j, k, first;

  i = p;
  j = r;
  for( k = p; k <= q; k++ )
  {
    first = 0;
    if( j <= q )
    {
      if( i < r )
      {
        if( array[tosort[i]] < array[tosort[j]] )
          first = 1;
      }
    }
    else
      first = 1;

    if( first )
    {
      sorted[k] = tosort[i];
      i++;
    }
    else
    {
      sorted[k] = tosort[j];
      j++;
    }
  }

  for( k = p; k <= q; k++ )
    tosort[k] = sorted[k];
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/




/*-=-=-=-=-=-=-=-=-=-=- Section Interpret Command Line -=-=-=-=-=-=-=-=-=-=-*/
/**
 * Parses and checks the command line.
 */
void interpretCommandLine( int argc, char **argv )
{
  parseCommandLine( argc, argv );
  
  checkOptions();
}

/**
 * Parses the command line.
 * For options, see printUsage.
 */
void parseCommandLine( int argc, char **argv )
{
  int index;

  index = 1;

  parseOptions( argc, argv, &index );
  
  parseParameters( argc, argv, &index );
}

/**
 * Parses only the options from the command line.
 */
// void parseOptions( int argc, char **argv, int *index )
// {
//   double dummy;

//   write_generational_statistics = 0;
//   write_generational_solutions  = 0;
//   print_verbose_overview        = 0;
//   print_FOSs_contents           = 0;
//   use_ilse                      = 0;

//   for( ; (*index) < argc; (*index)++ )
//   {
//     if( argv[*index][0] == '-' )
//     {
//       /* If it is a negative number, the option part is over */
//       if( sscanf( argv[*index], "%lf", &dummy ) && argv[*index][1] != '\0' )
//         break;

//       if( argv[*index][1] == '\0' )
//         optionError( argv, *index );
//       else if( argv[*index][2] != '\0' )
//         optionError( argv, *index );
//       else
//       {
//         switch( argv[*index][1] )
//         {
//           case '?': printUsage(); break;
//           case 'P': printAllInstalledProblems(); break;
//           case 'F': printAllInstalledFOSs(); break;
//           case 's': write_generational_statistics = 1; break;
//           case 'w': write_generational_solutions  = 1; break;
//           case 'v': print_verbose_overview        = 1; break;
//           case 'f': print_FOSs_contents           = 1; break;
//           case 'i': use_ilse                      = 1; break;
//           default : optionError( argv, *index );
//         }
//       }
//     }
//     else /* Argument is not an option, so option part is over */
//      break;
//   }
// }

/**
 * Writes the names of all installed problems to the standard output.
 */
void printAllInstalledProblems()
{
  int i, n;

  n = numberOfInstalledProblems();
  printf( "Installed optimization problems:\n" );
  for( i = 0; i < n; i++ )
    printf( "%3d: %s\n", i, installedProblemName( i ) );

  exit( 0 );
}

/**
 * Writes the names of all installed FOS structures to the standard output.
 */
void printAllInstalledFOSs()
{
  int i, n;

  n = numberOfInstalledFOSStructures();
  printf( "Installed FOS structures:\n" );
  for( i = 0; i < n; i++ )
    printf( "%3d: %s\n", i, installedFOSStructureName( i ) );

  exit( 0 );
}

/**
 * Informs the user of an illegal option and exits the program.
 */
void optionError( char **argv, int index )
{
  printf( "Illegal option: %s\n\n", argv[index] );

  printUsage();
}

/**
 * Parses only the EA parameters from the command line.
 */
// void parseParameters( int argc, char **argv, int *index )
// {
//   int noError;

//   if( (argc - *index) != 5 )
//   {
//     printf( "Number of parameters is incorrect, require 5 parameters (you provided %d).\n\n", (argc - *index) );

//     printUsage();
//   }

//   noError = 1;
//   noError = noError && sscanf( argv[*index+0], "%d", &problem_index );
//   noError = noError && sscanf( argv[*index+1], "%d", &number_of_parameters );
//   noError = noError && sscanf( argv[*index+2], "%d", &FOSs_structure_index );
//   noError = noError && sscanf( argv[*index+3], "%ld", &maximum_number_of_evaluations );
//   noError = noError && sscanf( argv[*index+4], "%ld", &maximum_number_of_milliseconds );

//   if( !noError )
//   {
//     printf( "Error parsing parameters.\n\n" );

//     printUsage();
//   }
// }

/**
 * Prints usage information and exits the program.
 */
void printUsage()
{
  printf("Usage: GOMEA [-?] [-P] [-F] [-s] [-w] [-v] [-f] [-i] pro dim fos eva mil\n");
  printf("   -?: Prints out this usage information.\n");
  printf("   -P: Prints out a list of all installed optimization problems.\n");
  printf("   -F: Prints out a list of all installed FOS structures.\n");
  printf("   -s: Enables computing and writing of statistics every generation.\n");
  printf("   -w: Enables writing of solutions and their fitnesses every generation.\n");
  printf("   -v: Enables verbose mode. Prints the settings before starting the run.\n");
  printf("   -f: Enables printing the contents of the FOS every generation.\n");
  printf("   -i: Enables Incremental Linkage Subset Evaluation (ILSE).\n");
  printf("\n");
  printf("  pro: Index of optimization problem to be solved (minimization).\n");
  printf("  dim: Number of parameters.\n");
  printf("  fos: Index of FOS structure to use.\n");
  printf("  eva: Maximum number of evaluations allowed (-1 for no limit).\n");
  printf("  mil: Maximum number of milliseconds allowed (-1 for no limit).\n");
  exit( 0 );
}

/**
 * Checks whether the selected options are feasible.
 */
// void checkOptions()
// {
//   if( number_of_parameters < 1 )
//   {
//     printf("\n");
//     printf("Error: number of parameters < 1 (read: %d). Require number of parameters >= 1.", number_of_parameters);
//     printf("\n\n");

//     exit( 0 );
//   }

//   if( installedProblemName( problem_index ) == NULL )
//   {
//     printf("\n");
//     printf("Error: unknown index for problem (read index %d).", problem_index );
//     printf("\n\n");

//     exit( 0 );
//   }

//   if( installedFOSStructureName( FOSs_structure_index ) == NULL )
//   {
//     printf("\n");
//     printf("Error: unknown index for FOS structure (read index %d).", FOSs_structure_index );
//     printf("\n\n");

//     exit( 0 );
//   }
// }

/**
 * Prints the settings as read from the command line.
 */
void printVerboseOverview(struct Environment *env)
{
  printf("### Settings ######################################\n");
  printf("#\n");
  printf("# Statistics writing every generation: %s\n", env->write_generational_statistics ? "enabled" : "disabled");
  printf("# Population file writing            : %s\n", env->write_generational_solutions ? "enabled" : "disabled");
  printf("# FOS contents printing              : %s\n", env->print_FOSs_contents ? "enabled" : "disabled");
  printf("# ILSE                               : %s\n", env->use_ilse ? "enabled" : "disabled" );
  printf("#\n");
  printf("###################################################\n");
  printf("#\n");
  printf("# Problem                            = %s\n", installedProblemName( env->problem_index ));
  printf("# Number of parameters               = %d\n", env->number_of_parameters);
  printf("# FOS structure                      = %s\n", installedFOSStructureName( env->FOSs_structure_index ));
  printf("# Maximum number of evaluations      = %ld\n", env->maximum_number_of_evaluations);
  printf("# Maximum number of milliseconds     = %ld\n", env->maximum_number_of_milliseconds);
  printf("# Random seed                        = %ld\n", (long) env->random_seed);
  printf("#\n");
  printf("###################################################\n");
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Random Numbers -=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns a random double, distributed uniformly between 0 and 1.
 */
double randomRealUniform01(struct Environment *env)
{
  int64_t n26, n27;
  double  result;

  env->random_seed_changing = (env->random_seed_changing * 0x5DEECE66DLLU + 0xBLLU) & ((1LLU << 48) - 1);
  n26                  = (int64_t)(env->random_seed_changing >> (48 - 26));
  env->random_seed_changing = (env->random_seed_changing * 0x5DEECE66DLLU + 0xBLLU) & ((1LLU << 48) - 1);
  n27                  = (int64_t)(env->random_seed_changing >> (48 - 27));
  result               = (((int64_t)n26 << 27) + n27) / ((double) (1LLU << 53));

  return( result );
}
        
/**
 * Returns a random integer, distributed uniformly between 0 and maximum.
 */
int randomInt(struct Environment *env, int maximum )
{
  int result;
  
  result = (int) (((double) maximum)*randomRealUniform01(env));
  
  return( result );
}

/**
 * Returns a random compact (using integers 0,1,...,n-1) permutation
 * of length n using the Fisher-Yates shuffle.
 */
int *randomPermutation(struct Environment *env, int n )
{
  int i, j, dummy, *result;

  result = (int *) Malloc( n*sizeof( int ) );
  for( i = 0; i < n; i++ )
    result[i] = i;

  for( i = n-1; i > 0; i-- )
  {
    j         = randomInt(env, i+1 );
    dummy     = result[j];
    result[j] = result[i];
    result[i] = dummy;
  }

  return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-= Section Problems -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns the name of an installed problem.
 */
char *installedProblemName( int index )
{
  switch( index )
  {
    case  0: return( (char *) "Sort" );
    default: 
      if (index <= 120)
        return ((char*) "Taillard Flowshop");                  
  }

  return( NULL );
}

/**
 * Returns the number of problems installed.
 */
int numberOfInstalledProblems()
{
  static int result = -1;
  
  if( result == -1 )
  {
    result = 0;
    while( installedProblemName( result ) != NULL )
      result++;
  }
  
  return( result );
}


/**
 * Calculate the fitness using the codomain values passed as input
 */
void calculateFitness(struct Environment *env, char *parameters, double *objective_value, double *constraint_value) 
{
  for(uint i = 0; i < env->input_parameters.M; i++) {
    //Calculate for each clique the solution substring for this clique, as an index into an array of these substrings
    int clique_substring_as_index = 0;

    //Store pointer to the clique, for convenience
    uint* clique = env->cliques[i];
    for(int j = env->input_parameters.k - 1; j >= 0; j--) {
      //TODO: maybe store parameters as uint values? Then we can skip this cast
      //Go over each variable index in the clique and for each one, take the bit value from the solution string and add it to the clique substring.
      // char value = parameters[clique[j]];
      // printf("value of parameter %d: %d", clique[j], value);
      clique_substring_as_index += (uint)(parameters[clique[j]]) << (env->input_parameters.k - j - 1);
    }

    //Add the fitness contribution of this clique/subfunction
    *objective_value += env->codomain_values[i][clique_substring_as_index];
  }
}

/**
 * Computes the value of the single objective and the sum of all
 * constraint violations, using the clique tree read from a file or passed.
 */
void problemEvaluation(struct Environment *env, int index, char *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, char *parameters_before, double objective_value_before, double constraint_value_before, int GOMEA_index )
{
  short same;
  int i, j;

  /* Count the evaluation */
  env->number_of_evaluations++;
  if( number_of_touched_parameters == 0 )
    env->number_of_bit_flip_evaluations += env->number_of_parameters;
  else
  {
    for( i = 0; i < number_of_touched_parameters; i++ )
    {
      if( parameters[touched_parameters_indices[i]] != parameters_before[touched_parameters_indices[i]] )
        env->number_of_bit_flip_evaluations++;
    }
  }

  /* Do the actual evaluation */
  *objective_value  = 0.0;
  *constraint_value = 0.0;

  calculateFitness(env, parameters, objective_value, constraint_value);
  // switch( index )
  // {
  //   case  0: sortFunctionProblemEvaluation( parameters_as_integer_permutation, objective_value, constraint_value ); break;
  //   default:
  //     taillardFlowshopProblemEvaluation(parameters_as_integer_permutation, objective_value, constraint_value); break;
  // }

  // free( parameters_as_integer_permutation );

  /* Check the VTR */
  if( !env->vosostr_hit_status )
  {
    if( env->vtr_exists > 0 )
    {
      if( ((*constraint_value) == 0) && (equalFitness((*objective_value), 0.0, env->global_optima_score, 0.0) || betterFitness((*objective_value), 0.0, env->global_optima_score, 0.0) ) ) //((*objective_value) >= env->vtr)  )
      {
        env->vosostr_hit_status = 1;
      }
    }
  }

  /* Check the SOSTR */
  if( !env->vosostr_hit_status )
  {
    if( env->sostr_exists )
    {
      for( i = 0; i < env->number_of_solutions_in_sostr; i++ )
      {
        same = 1;
        for( j = 0; j < env->number_of_parameters; j++ )
        {
          if( parameters[j] != env->sostr[i][j] )
          {
            same = 0;
            break;
          }
        }
        if( same )
        {
          env->vosostr_hit_status = 1;
          break;
        }
      }
    }
  }

  /* Check the VOSOSTR */
  if( env->vosostr_hit_status == 1 )
  {
    env->vosostr_hit_status                     = 2;
    env->vosostr_hitting_time                   = getMilliSecondsRunningAfterInit( env );
    env->vosostr_number_of_evaluations          = env->number_of_evaluations;
    env->vosostr_number_of_bit_flip_evaluations = env->number_of_bit_flip_evaluations;

    env->ltga_result.global_optima_found       = 1;
    env->ltga_result.pop_size_glob_opt          = env->population_sizes[GOMEA_index];
    // writeRunningTime( env, (char *) "vosostr_hitting_time.dat" );
  }

  /* Update elitist solution */
  if( (env->number_of_evaluations == 1) || betterFitness( *objective_value, *constraint_value, env->elitist_solution_objective_value, env->elitist_solution_constraint_value ) )
  {
    for( i = 0; i < env->number_of_parameters; i++ )
      env->elitist_solution[i] = parameters[i];

    env->elitist_solution_objective_value                = *objective_value;
    env->elitist_solution_constraint_value               = *constraint_value;
    env->elitist_solution_hitting_time                   = getMilliSecondsRunningAfterInit( env );
    env->elitist_solution_number_of_evaluations          = env->number_of_evaluations;
    env->elitist_solution_number_of_bit_flip_evaluations = env->number_of_bit_flip_evaluations;


    //Update LTGA elitist solution stats
    for( i = 0; i < env->number_of_parameters; i++ )
      env->ltga_result.elitist_solution[i] = parameters[i];

    env->ltga_result.elitist_fitness                     = *objective_value;
    env->ltga_result.elitist_time                        = env->elitist_solution_hitting_time;
    env->ltga_result.elitist_evaluations                 = env->number_of_evaluations;



    // writeRunningTime( env,  (char *) "elitist_solution_hitting_time.dat" );

    // writeElitistSolution( env );

    //Update Elitist result here
  }
  /* Exit early, depending on VOSOSTR status */
  if( env->vosostr_hit_status != 0 )
    env->halt_execution = 1;
    // exit(0);
  if( (env->number_of_evaluations >= env->maximum_number_of_evaluations) )
    env->halt_execution = 1;
    // exit(0);
}


// /**
//  * Computes the value of the single objective
//  * and the sum of all constraint violations
//  * function.
//  */
// void installedProblemEvaluation( int index, double *parameters, double *objective_value, double *constraint_value, int number_of_touched_parameters, int *touched_parameters_indices, double *parameters_before, double objective_value_before, double constraint_value_before )
// {
//   short same;
//   int   i, j, *parameters_as_integer_permutation;

//   parameters_as_integer_permutation = randomKeysToIntegerPermutation( parameters );

//   /* Count the evaluation */
//   number_of_evaluations++;
//   if( number_of_touched_parameters == 0 )
//     number_of_bit_flip_evaluations += number_of_parameters;
//   else
//   {
//     for( i = 0; i < number_of_touched_parameters; i++ )
//     {
//       if( parameters[touched_parameters_indices[i]] != parameters_before[touched_parameters_indices[i]] )
//         number_of_bit_flip_evaluations++;
//     }
//   }

//   /* Do the actual evaluation */
//   *objective_value  = 0.0;
//   *constraint_value = 0.0;
//   switch( index )
//   {
//     case  0: sortFunctionProblemEvaluation( parameters_as_integer_permutation, objective_value, constraint_value ); break;
//     default:
//       taillardFlowshopProblemEvaluation(parameters_as_integer_permutation, objective_value, constraint_value); break;
//   }

//   free( parameters_as_integer_permutation );

//   /* Check the VTR */
//   if( !vosostr_hit_status )
//   {
//     if( vtr_exists > 0 )
//     {
//       if( ((*constraint_value) == 0) && ((*objective_value) >= vtr)  )
//       {
//         vosostr_hit_status = 1;
//       }
//     }
//   }

//   /* Check the SOSTR */
//   if( !vosostr_hit_status )
//   {
//     if( sostr_exists )
//     {
//       for( i = 0; i < number_of_solutions_in_sostr; i++ )
//       {
//         same = 1;
//         for( j = 0; j < number_of_parameters; j++ )
//         {
//           if( parameters[j] != sostr[i][j] )
//           {
//             same = 0;
//             break;
//           }
//         }
//         if( same )
//         {
//           vosostr_hit_status = 1;
//           break;
//         }
//       }
//     }
//   }

//   /* Check the VOSOSTR */
//   if( vosostr_hit_status == 1 )
//   {
//     vosostr_hit_status                     = 2;
//     vosostr_hitting_time                   = getMilliSecondsRunningAfterInit();
//     vosostr_number_of_evaluations          = number_of_evaluations;
//     vosostr_number_of_bit_flip_evaluations = number_of_bit_flip_evaluations;

//     writeRunningTime( (char *) "vosostr_hitting_time.dat" );
//   }

//   /* Update elitist solution */
//   if( (number_of_evaluations == 1) || betterFitness( *objective_value, *constraint_value, elitist_solution_objective_value, elitist_solution_constraint_value ) )
//   {
//     for( i = 0; i < number_of_parameters; i++ )
//       elitist_solution[i] = parameters[i];

//     elitist_solution_objective_value                = *objective_value;
//     elitist_solution_constraint_value               = *constraint_value;
//     elitist_solution_hitting_time                   = getMilliSecondsRunningAfterInit();
//     elitist_solution_number_of_evaluations          = number_of_evaluations;
//     elitist_solution_number_of_bit_flip_evaluations = number_of_bit_flip_evaluations;

//     writeRunningTime( (char *) "elitist_solution_hitting_time.dat" );

//     writeElitistSolution();
//   }

//   /* Exit early, depending on VOSOSTR status */
//   if( vosostr_hit_status != 0 )
//     exit(0);
//   if( (number_of_evaluations >= maximum_number_of_evaluations) )
//     exit(0);
// }

// int *randomKeysToIntegerPermutation( double *parameters )
// {
//   int *result;

//   result = mergeSortDoublesIncreasing( parameters, number_of_parameters );
  
//   return( result );
// }

// void sortFunctionProblemEvaluation( int *parameters, double *objective_value, double *constraint_value )
// {
//   int    i, j;
//   double result;

//   result = 0.0;
//   for( i = 0; i < number_of_parameters; i++ )
//     for( j = i+1; j < number_of_parameters; j++ )
//       result += parameters[i] > parameters[j] ? 0 : 1;

//   *objective_value  = result;
//   *constraint_value = 0;
// }

// void taillardFlowshopProblemEvaluation( int *parameters, double *objective_value, double *constraint_value)
// {
//   int c[njobs][nmach];
//   int i, m;

//   c[0][0] = ptime[parameters[0]][0];
//   for (m = 1; m < nmach; m++)
//     c[0][m] = c[0][m - 1] + ptime[parameters[0]][m];

//   for (i = 1; i < njobs; i++) {
//     c[i][0] = c[i - 1][0] + ptime[parameters[i]][0];
//     for (m = 1; m < nmach; m++)
//       c[i][m] = MAX(c[i - 1][m], c[i][m - 1]) + ptime[parameters[i]][m];
//   }

//   *objective_value = -c[njobs - 1][nmach - 1];
//   *objective_value = 0;
//   for (i = 0; i < njobs; i++)
//     *objective_value += -c[i][nmach - 1];
//   *constraint_value = 0;
// }

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Initialize -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Performs initialization for a single GOMEA.
 */
void initializeNewGOMEA(struct Environment *env)
{
  int gomea_index;

  if( env->number_of_GOMEAs == 0 )
  {
    env->populations                  = (char ***) Malloc( env->maximum_number_of_GOMEAs*sizeof( char ** ) );
    env->objective_values             = (double **) Malloc( env->maximum_number_of_GOMEAs*sizeof( double * ) );
    env->constraint_values            = (double **) Malloc( env->maximum_number_of_GOMEAs*sizeof( double * ) );
    env->offsprings                   = (char ***) Malloc( env->maximum_number_of_GOMEAs*sizeof( char ** ) );
    env->objective_values_offsprings  = (double **) Malloc( env->maximum_number_of_GOMEAs*sizeof( double * ) );
    env->constraint_values_offsprings = (double **) Malloc( env->maximum_number_of_GOMEAs*sizeof( double * ) );
    env->average_objective_values     = (double *) Malloc( env->maximum_number_of_GOMEAs*sizeof( double ) );
    env->average_constraint_values    = (double *) Malloc( env->maximum_number_of_GOMEAs*sizeof( double ) );
    env->terminated                   = (char *) Malloc( env->maximum_number_of_GOMEAs*sizeof( char ) );
    env->dependency_matrices          = (double ***) Malloc( env->maximum_number_of_GOMEAs*sizeof( double ** ) );
    env->FOSs                         = (int ***) Malloc( env->maximum_number_of_GOMEAs*sizeof( int ** ) );
    env->FOSs_number_of_indices       = (int **) Malloc( env->maximum_number_of_GOMEAs*sizeof( int * ) );
    env->FOSs_length                  = (int *) Malloc( env->maximum_number_of_GOMEAs*sizeof( int ) );
    env->objective_values_best_of_generation  = (double *) Malloc( env->maximum_number_of_GOMEAs*sizeof( double ) );
    env->constraint_values_best_of_generation = (double *) Malloc( env->maximum_number_of_GOMEAs*sizeof( double ) );

    env->population_sizes[env->number_of_GOMEAs]                   = env->base_population_size;
    env->number_of_subgenerations_per_GOMEA[env->number_of_GOMEAs] = 1;
  }
  else
  {
    env->population_sizes[env->number_of_GOMEAs]                   = 2*env->population_sizes[env->number_of_GOMEAs-1];
    env->number_of_subgenerations_per_GOMEA[env->number_of_GOMEAs] = 1;
    for( gomea_index = 0; gomea_index < env->number_of_GOMEAs; gomea_index++ )
      env->number_of_subgenerations_per_GOMEA[gomea_index] *= env->number_of_subgenerations_per_GOMEA_factor;
  }

  env->terminated[env->number_of_GOMEAs]              = 0;
  env->no_improvement_stretchs[env->number_of_GOMEAs] = 0;
  env->FOSs[env->number_of_GOMEAs]                    = NULL;

  initializeNewGOMEAMemory(env);

  initializeNewGOMEAPopulationAndFitnessValues(env);

  env->number_of_GOMEAs++;
}

/**
 * Initializes the memory for a single GOMEA.
 */
void initializeNewGOMEAMemory(struct Environment *env)
{
  int i;

  env->populations[env->number_of_GOMEAs]                  = (char **) Malloc( env->population_sizes[env->number_of_GOMEAs]*sizeof( char * ) );
  env->objective_values[env->number_of_GOMEAs]             = (double *) Malloc( env->population_sizes[env->number_of_GOMEAs]*sizeof( double ) );
  env->constraint_values[env->number_of_GOMEAs]            = (double *) Malloc( env->population_sizes[env->number_of_GOMEAs]*sizeof( double ) );
  env->offsprings[env->number_of_GOMEAs]                   = (char **) Malloc( env->population_sizes[env->number_of_GOMEAs]*sizeof( char * ) );
  env->objective_values_offsprings[env->number_of_GOMEAs]  = (double *) Malloc( env->population_sizes[env->number_of_GOMEAs]*sizeof( double ) );
  env->constraint_values_offsprings[env->number_of_GOMEAs] = (double *) Malloc( env->population_sizes[env->number_of_GOMEAs]*sizeof( double ) );

  for( i = 0; i < env->population_sizes[env->number_of_GOMEAs]; i++ )
    env->populations[env->number_of_GOMEAs][i] = (char *) Malloc( env->number_of_parameters*sizeof( char ) );

  for( i = 0; i < env->population_sizes[env->number_of_GOMEAs]; i++ )
    env->offsprings[env->number_of_GOMEAs][i] = (char *) Malloc( env->number_of_parameters*sizeof( char ) );

  env->dependency_matrices[env->number_of_GOMEAs] = (double **) Malloc( env->number_of_parameters*sizeof( double * ) );
  for( i = 0; i < env->number_of_parameters; i++ )
    (env->dependency_matrices[env->number_of_GOMEAs])[i] = (double *) Malloc( env->number_of_parameters*sizeof( double ) );

  env->FOSs[env->number_of_GOMEAs] = NULL;
}

/**
 * Initializes the population and the objective values by randomly
 * generation n solutions.
 */
void initializeNewGOMEAPopulationAndFitnessValues(struct Environment *env)
{
  int    i, j;
  double obj, con;

  env->objective_values_best_of_generation[env->number_of_GOMEAs]  = -1e+308;
  env->constraint_values_best_of_generation[env->number_of_GOMEAs] = 1e+308;
  for( i = 0; i < env->population_sizes[env->number_of_GOMEAs]; i++ )
  {
    for( j = 0; j < env->number_of_parameters; j++ )
      env->populations[env->number_of_GOMEAs][i][j] = (randomInt(env, 2 ) == 1) ? TRUE : FALSE;
      //populations[number_of_GOMEAs][i][j] = randomRealUniform01();

    problemEvaluation(env, env->problem_index, env->populations[env->number_of_GOMEAs][i], &obj, &con, 0, NULL, NULL, 0, 0, env->number_of_GOMEAs );
    env->objective_values[env->number_of_GOMEAs][i]  = obj;
    env->constraint_values[env->number_of_GOMEAs][i] = con;

    if( betterFitness( env->objective_values[env->number_of_GOMEAs][i], env->constraint_values[env->number_of_GOMEAs][i], env->objective_values_best_of_generation[env->number_of_GOMEAs], env->constraint_values_best_of_generation[env->number_of_GOMEAs] ) )
    {
      env->objective_values_best_of_generation[env->number_of_GOMEAs]  = env->objective_values[env->number_of_GOMEAs][i];
      env->constraint_values_best_of_generation[env->number_of_GOMEAs] = env->constraint_values[env->number_of_GOMEAs][i];
    }
  }
}

// }
/**
 * Checks to see if files exists with values and solutions to reach.
 */
void initializeValueAndSetOfSolutionsToReach(struct Environment *env)
{
  env->vtr_exists = 0;
  if( initializeValueToReach(env) )
    env->vtr_exists = 1;

  env->sostr_exists = 0;
  if( initializeSetOfSolutionsToReach(env) )
    env->sostr_exists = 1;
}

/**
 * Attempts to read the value to reach.
 */
short initializeValueToReach(struct Environment *env)
{
  char  c, string[100000], filename[1000];
  int   i;
  FILE *file;

  sprintf( filename, "vtr.txt" );
  file = fopen( filename, "r" );
  if( file == NULL )
    return( 0 );

  i    = 0;
  c    = fgetc( file );
  while( (c != '\n') && (c != EOF) && (c != ' ') )
  {
    string[i] = (char) c;
    c         = fgetc( file );
    i++;
  }
  string[i] = '\0';
  sscanf( string, "%le", &(env->vtr));

  fclose( file );

  return( 1 );
}

/**
 * Attempts to read assets of solutions to reach.
 */
short initializeSetOfSolutionsToReach(struct Environment *env)
{
  char  c, string[100000], filename[1000];
  int   i, j;
  FILE *file;

  env->number_of_solutions_in_sostr = 0;

  sprintf( filename, "sostr.txt" );
  file = fopen( filename, "r" );
  if( file == NULL )
    return( 0 );

  do
  {
    c = fgetc( file );
    if( (c == '0') || (c == '1') )
    {
      env->number_of_solutions_in_sostr++;
      do
      {
        c = fgetc( file );
      } while( (c == '0') || (c == '1') );
    }
  }
  while( c != EOF );
  fclose( file );

  env->sostr = (char **) Malloc( env->number_of_solutions_in_sostr*sizeof( char * ) );
  for( i = 0; i < env->number_of_solutions_in_sostr; i++ )
    env->sostr[i] = (char *) Malloc( env->number_of_parameters*sizeof( char ) );

  file = fopen( filename, "r" );
  i    = 0;
  do
  {
    c = fgetc( file );
    if( c == EOF )
      break;
    if( !((c == '0') || (c == '1')) )
      continue;
    j = 0;
    while( (c != '\n') && (c != EOF) && (c != ' ') )
    {
      string[j] = (char) c;
      c         = fgetc( file );
      j++;
    }
    if( j != env->number_of_parameters )
    {
      printf("Error while reading %s: the number of parameters in at least one of the solutions does not match the number of parameters on the commandline (%d != %d)\n",filename,j,env->number_of_parameters);
      exit(0);
    }
    for( j = 0; j < env->number_of_parameters; j++ )
      env->sostr[i][j] = string[j] == '0' ? 0 : 1;
  }
  while( c != EOF );

  fclose( file );

  return( 1 );
}

/**
 * Initializes the pseudo-random number generator.
 */
void initializeRandomNumberGenerator(struct Environment *env)
{
  struct timeval tv;
  struct tm *timep;

  while( env->random_seed_changing == 0 )
  {
    gettimeofday( &tv, NULL );
    timep = localtime (&tv.tv_sec);
    env->random_seed_changing = timep->tm_hour * 3600 * 1000 + timep->tm_min * 60 * 1000 + timep->tm_sec * 1000 + tv.tv_usec / 1000;
  }

  env->random_seed = env->random_seed_changing;
}


// void initializeProblem( int index )
// {
//   char fn[100];

//   switch( index )
//   {
//     case  0: break;
//     default:
//       switch((index - 1) / 10) {
//         case 0: 
//           sprintf(fn, "tai20_5.txt");
//           loadTaillardFlowshop(20, 5, fn, index); break;
//         case 1: 
//           sprintf(fn, "tai20_10.txt");
//           loadTaillardFlowshop(20, 10, fn, index - 10); break;
//         case 2: 
//           sprintf(fn, "tai20_20.txt");
//           loadTaillardFlowshop(20, 20, fn, index - 20); break;
//         case 3: 
//           sprintf(fn, "tai50_5.txt");
//           loadTaillardFlowshop(50, 5, fn, index - 30); break;
//         case 4: 
//           sprintf(fn, "tai50_10.txt");
//           loadTaillardFlowshop(50, 10, fn, index - 40); break;
//         case 5: 
//           sprintf(fn, "tai50_20.txt");
//           loadTaillardFlowshop(50, 20, fn, index - 50); break;
//         case 6: 
//           sprintf(fn, "tai100_5.txt");
//           loadTaillardFlowshop(100, 5, fn, index - 60); break;
//         case 7: 
//           sprintf(fn, "tai100_10.txt");
//           loadTaillardFlowshop(100, 10, fn, index - 70); break;
//         case 8: 
//           sprintf(fn, "tai100_20.txt");
//           loadTaillardFlowshop(100, 20, fn, index - 80); break;
//         case 9: 
//           sprintf(fn, "tai200_10.txt");
//           loadTaillardFlowshop(200, 10, fn, index - 90); break;
//         case 10: 
//           sprintf(fn, "tai200_20.txt");
//           loadTaillardFlowshop(200, 20, fn, index - 100); break;
//         case 11: 
//           sprintf(fn, "tai500_20.txt");
//           loadTaillardFlowshop(500, 20, fn, index - 110); break;
//      }
//   }
// }

// void loadTaillardFlowshop(int J, int M, char *fn, int inst) {
//   int i, j, c, m, t;
//   FILE *fp;
//   char *line;
//   size_t maxlinelen = 2000;

//   njobs = J;
//   nmach = M;
//   ptime = (int **) Malloc(sizeof(int *) * njobs);
//   for (j = 0; j < njobs; j++)
//     ptime[j] = (int *) Malloc(sizeof(int) * nmach);

//   fp = fopen(fn, "r");
//   line = (char *) Malloc(sizeof(char) * (maxlinelen + 1));

//   for (i = 1; i < inst; i++) 
//     for (c = 0; c < (3 + nmach); c++)
//       if (getline(&line, &maxlinelen, fp) == -1) {
//         printf("%s\n", line);
//         printf("Error: Something wrong with the format of file %s", fn);
//         exit(0);
//       }

//   for (c = 0; c < 3; c++)
//     if (getline(&line, &maxlinelen, fp) == -1) {
//       printf("Error: Something wrong with the format of file %s", fn);
//       exit(0);
//     }

//   for (m = 0; m < nmach; m++) 
//     for (j = 0; j < njobs; j++) {
//       fscanf(fp, "%d", &t);
//       ptime[j][m] = t;
//     }
  
//   fclose(fp);
//   free(line);
// }

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*=-=-=-=-=-=-=-=-=-=-= Section Survivor Selection =-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Determines the solutions that finally survive the generation (offspring only).
 */
void selectFinalSurvivorsSpecificGOMEA(struct Environment *env, int gomea_index )
{
  int    i, j;
  double objective_values_best_of_generation_before, constraint_values_best_of_generation_before;

  objective_values_best_of_generation_before  = env->objective_values_best_of_generation[gomea_index];
  constraint_values_best_of_generation_before = env->constraint_values_best_of_generation[gomea_index];

  for( i = 0; i < env->population_sizes[gomea_index]; i++ )
  {
    for( j = 0; j < env->number_of_parameters; j++ )
      env->populations[gomea_index][i][j] = env->offsprings[gomea_index][i][j];
    env->objective_values[gomea_index][i]  = env->objective_values_offsprings[gomea_index][i];
    env->constraint_values[gomea_index][i] = env->constraint_values_offsprings[gomea_index][i];

    if( betterFitness( env->objective_values[gomea_index][i], env->constraint_values[gomea_index][i], env->objective_values_best_of_generation[gomea_index], env->constraint_values_best_of_generation[gomea_index] ) )
    {
      env->objective_values_best_of_generation[gomea_index]  = env->objective_values[gomea_index][i];
      env->constraint_values_best_of_generation[gomea_index] = env->constraint_values[gomea_index][i];
    }
  }

  if( !betterFitness( env->objective_values_best_of_generation[gomea_index], env->constraint_values_best_of_generation[gomea_index], objective_values_best_of_generation_before, constraint_values_best_of_generation_before ) )
    env->no_improvement_stretchs[gomea_index]++;
  else
    env->no_improvement_stretchs[gomea_index] = 0;
}

/**
 * Returns 1 if x is better than y, 0 otherwise.
 * x is not better than y unless:
 * - x and y are both infeasible and x has a smaller sum of constraint violations, or
 * - x is feasible and y is not, or
 * - x and y are both feasible and x has a larger objective value than y
 */
char betterFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y )
{
  char result;
  
  result = 0;

  if( constraint_value_x > 0 ) /* x is infeasible */
  {
    if( constraint_value_y > 0 ) /* Both are infeasible */
    {
      if( constraint_value_x < constraint_value_y )
       result = 1;
    }
  }
  else /* x is feasible */
  {
    if( constraint_value_y > 0 ) /* x is feasible and y is not */
      result = 1;
    else /* Both are feasible */
    {
      if( objective_value_x > objective_value_y && fabs(objective_value_x - objective_value_y) >= FITNESS_EPSILON )
        result = 1;
    }
  }

  return( result );
}

/**
 * Returns 1 if x is equally preferable to y, 0 otherwise.
 */
char equalFitness( double objective_value_x, double constraint_value_x, double objective_value_y, double constraint_value_y )
{
  char result;
  
  result = 0;

  if( (constraint_value_x == constraint_value_y) && ((objective_value_x == objective_value_y) || fabs(objective_value_x - objective_value_y) < FITNESS_EPSILON ) )
    result = 1;

  return( result );
}

void computeAverageFitnessSpecificGOMEA( struct Environment *env, int GOMEA_index )
{
  int i;

  env->average_objective_values[GOMEA_index]  = 0;
  env->average_constraint_values[GOMEA_index] = 0;
  for( i = 0; i < env->population_sizes[GOMEA_index]; i++ )
  {
    env->average_objective_values[GOMEA_index]  += env->objective_values[GOMEA_index][i];
    env->average_constraint_values[GOMEA_index] += env->constraint_values[GOMEA_index][i];
  }
  env->average_objective_values[GOMEA_index] /= (double) (env->population_sizes[GOMEA_index]);
  env->average_constraint_values[GOMEA_index] /= (double) (env->population_sizes[GOMEA_index]);
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Output =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Writes (appends) statistics about the current generation to a
 * file named "statistics.dat".
 */
void writeGenerationalStatistics(struct Environment *env)
{
  char    filename[1000], string[10000];
  int     i, n, gomea_index;
  double  objective_avg, objective_var, objective_best, objective_worst,
          constraint_avg, constraint_var, constraint_best, constraint_worst;
  FILE   *file;

  /* First compute the statistics */
  /* Average, best and worst */ 
  objective_avg    = 0.0;
  constraint_avg   = 0.0;
  objective_best   = env->objective_values[0][0];
  objective_worst  = env->objective_values[0][0];
  constraint_best  = env->constraint_values[0][0];
  constraint_worst = env->constraint_values[0][0];
  n                = 0;
  for( gomea_index = 0; gomea_index < env->number_of_GOMEAs; gomea_index++ )
  {
    for( i = 0; i < env->population_sizes[gomea_index]; i++ )
    {
      objective_avg += env->objective_values[gomea_index][i];
      constraint_avg += env->constraint_values[gomea_index][i];
      if( betterFitness( env->objective_values[gomea_index][i], env->constraint_values[gomea_index][i], objective_best, constraint_best ) )
      {
        objective_best  = env->objective_values[gomea_index][i];
        constraint_best = env->constraint_values[gomea_index][i];
      }
      if( betterFitness( objective_worst, constraint_worst, env->objective_values[gomea_index][i], env->constraint_values[gomea_index][i] ) )
      {
        objective_worst  = env->objective_values[gomea_index][i];
        constraint_worst = env->constraint_values[gomea_index][i];
      }
      n++;
    }
  }
  objective_avg = objective_avg / ((double) n);
  constraint_avg = constraint_avg / ((double) n);

  /* Variance */
  objective_var  = 0.0;
  constraint_var = 0.0;
  for( gomea_index = 0; gomea_index < env->number_of_GOMEAs; gomea_index++ )
  {
    for( i = 0; i < env->population_sizes[gomea_index]; i++ )
    {
      objective_var += (env->objective_values[gomea_index][i] - objective_avg)*(env->objective_values[gomea_index][i] - objective_avg);
      constraint_var += (env->constraint_values[gomea_index][i] - constraint_avg)*(env->constraint_values[gomea_index][i] - constraint_avg);
    }
  }
  objective_var = objective_var / ((double) n);
  constraint_var = constraint_var / ((double) n);

  if( objective_var <= 0.0 )
     objective_var = 0.0;
  if( constraint_var <= 0.0 )
     constraint_var = 0.0;

  /* Then write them */
  sprintf( filename, "statistics.dat" );
  file = NULL;
  if( env->number_of_generations == 0 )
  {
    file = fopen( filename, "w" );

    sprintf( string, "# Generation Evaluations      Average-obj.     Variance-obj.         Best-obj.        Worst-obj.        Elite-obj.      Average-con.     Variance-con.         Best-con.        Worst-con.        Elite-con.\n");
    fputs( string, file );
  }
  else
    file = fopen( filename, "a" );

  sprintf( string, "  %10d %11ld %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e %17.10e\n", env->number_of_generations, env->number_of_evaluations, objective_avg, objective_var, objective_best, objective_worst, env->elitist_solution_objective_value, constraint_avg, constraint_var, constraint_best, constraint_worst, env->elitist_solution_constraint_value );
  fputs( string, file );

  fclose( file );
}

/**
 * Writes the solutions to various files. The filenames
 * contain the generation. If the flag is_final_generation
 * is set the generation number in the filename
 * is replaced with the word "final".
 *
 * populations_xxxxx_generation_xxxxx.dat: the populations
 * offsprings_xxxxx_generation_xxxxx.dat : the offsprings
 */
void writeGenerationalSolutions( struct Environment *env, char is_final_generation )
{
  char  filename[1000], string[10000];
  int   i, j, gomea_index;
  FILE *file;

  /* Populations */
  if( is_final_generation )
    sprintf( filename, "populations_generation_final.dat" );
  else
    sprintf( filename, "populations_generation_%05d.dat", env->number_of_generations );
  file = fopen( filename, "w" );

  for( gomea_index = 0; gomea_index < env->number_of_GOMEAs; gomea_index++ )
  {
    for( i = 0; i < env->population_sizes[gomea_index]; i++ )
    {
      for( j = 0; j < env->number_of_parameters; j++ )
      {
        sprintf( string, "%u ", env->populations[gomea_index][i][j] );
        fputs( string, file );
      }
      sprintf( string, "     %17.10e %17.10e\n", env->objective_values[gomea_index][i], env->constraint_values[gomea_index][i] );
      fputs( string, file );
    }
  }
  fclose( file );
  
  /* Offsprings */
  if( (env->number_of_generations > 0) && (!is_final_generation) )
  {
    sprintf( filename, "offsprings_generation_%05d.dat", env->number_of_generations-1 );
    file = fopen( filename, "w" );

    for( gomea_index = 0; gomea_index < env->number_of_GOMEAs; gomea_index++ )
    {
      for( i = 0; i < env->population_sizes[gomea_index]; i++ )
      {
        for( j = 0; j < env->number_of_parameters; j++ )
        {
          sprintf( string, "%u", env->offsprings[gomea_index][i][j] );
          fputs( string, file );
        }
        sprintf( string, "     %17.10e %17.10e\n", env->objective_values_offsprings[gomea_index][i], env->constraint_values_offsprings[gomea_index][i] );
        fputs( string, file );
      }
    }

    fclose( file );
  }
}

void writeRunningTime(struct Environment *env, char *filename )
{
  char  string[10000];
  FILE *file;

  file = fopen( filename, "w" );
  sprintf( string, "# Column 1: Total number of milliseconds.\n");
  fputs( string, file );
  sprintf( string, "# Column 2: Total number of milliseconds after initialization.\n");
  fputs( string, file );
  sprintf( string, "# Column 3: Total number of evaluations.\n"); 
  fputs( string, file );
  sprintf( string, "# Column 4: Total number of bit-flip evaluations.\n"); 
  fputs( string, file );
#ifdef OS_WIN
  sprintf( string, "%ld %ld %ld %I64d\n", getMilliSecondsRunning(), getMilliSecondsRunningAfterInit(), number_of_evaluations, number_of_bit_flip_evaluations );
#else
  sprintf( string, "%ld %ld %ld %lld\n", getMilliSecondsRunning( env ), getMilliSecondsRunningAfterInit( env ), env->number_of_evaluations, env->number_of_bit_flip_evaluations );
#endif
  fputs( string, file );
  fclose( file );
}

void writeElitistSolution(struct Environment *env)
{
  char       string[10000];
  int        i;//, *elitist_solution_as_integer_permutation;
  FILE      *file;

  //elitist_solution_as_integer_permutation = randomKeysToIntegerPermutation( elitist_solution );
  file = fopen( "elitist_solution.dat", "w" );
  sprintf( string, "# Column 1: Solution.\n");
  fputs( string, file );
  sprintf( string, "# Column 2: Objective value.\n");
  fputs( string, file );
  sprintf( string, "# Column 3: Constraint value.\n");
  fputs( string, file );
  for( i = 0; i < env->number_of_parameters; i++ )
  {
    // sprintf( string, "%d ", elitist_solution_as_integer_permutation[i] );
    sprintf( string, "%d ", env->elitist_solution[i]);
    fputs( string, file );
  }
  sprintf( string, "     %17.10e %17.10e\n", env->elitist_solution_objective_value, env->elitist_solution_constraint_value );
  fputs( string, file );
  fclose( file );
  // free( elitist_solution_as_integer_permutation );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Termination -=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Returns 1 if termination should be enforced at the end of a generation, 0 otherwise.
 */
char checkTermination(struct Environment *env)
{
  if( checkNumberOfEvaluationsTerminationCondition(env) )
    return( 1 );
  
  if( checkNumberOfMilliSecondsTerminationCondition(env) )
    return( 1 );

  if( checkVOSOSTRTerminationCondition(env) )
      return( 1 );

  return( 0 );
}

/**
 * Returns 1 if the maximum number of evaluations
 * has been reached, 0 otherwise.
 */
char checkNumberOfEvaluationsTerminationCondition(struct Environment *env)
{
  if( (env->maximum_number_of_evaluations >= 0) && (env->number_of_evaluations >= env->maximum_number_of_evaluations) )
    return( 1 );

  return( 0 );
}

/**
 * Returns 1 if the value-to-reach has been reached.
 */
char checkVOSOSTRTerminationCondition(struct Environment *env)
{
  if( env->vosostr_hit_status > 0 )
    return( 1 );

  return( 0 );
}

/**
 * Returns 1 if the maximum number of milliseconds
 * has passed, 0 otherwise.
 */
char checkNumberOfMilliSecondsTerminationCondition(struct Environment *env)
{
  if( (env->maximum_number_of_milliseconds >= 0) && (getMilliSecondsRunning( env ) > env->maximum_number_of_milliseconds) )
    return( 1 );

  return( 0 );
}

char checkTerminationSpecificGOMEA( struct Environment *env, int GOMEA_index )
{
  int i, j;

  // EDIT: NOT PREMATURELY STOPPING SMALLER POPULATIONS (PETER)
  // Tobias Edit: Do allow use of permature stopping

  if(env->use_premature_stopping) {
    for( i = GOMEA_index+1; i < env->number_of_GOMEAs; i++ )
    {
      if( betterFitness( env->average_objective_values[i], env->average_constraint_values[i], env->average_objective_values[GOMEA_index], env->average_constraint_values[GOMEA_index] ) )
      {
        env->minimum_GOMEA_index = GOMEA_index+1;

        return( 1 );
      }
    }
  }


  if(env->use_fitness_variance_tolerance && checkFitnessVarianceTermination(env, GOMEA_index) ) {
    return( 1 );
  }

  for( i = 1; i < env->population_sizes[GOMEA_index]; i++ )
  {
    for( j = 0; j < env->number_of_parameters; j++ )
    {
      if( env->populations[GOMEA_index][i][j] != env->populations[GOMEA_index][0][j] )
        return( 0 );
    }
  }

  return( 1 );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/**
 * Checks whether the fitness variance
 * has become too small (user-defined tolerance).
 */
char checkFitnessVarianceTermination(struct Environment *env, int GOMEA_index )
{
  int    i;
  double objective_avg, objective_var;
  
  objective_avg = 0.0;
  for( i = 0; i < env->population_sizes[GOMEA_index]; i++ )
    objective_avg  += env->objective_values[GOMEA_index][i];
  objective_avg  = objective_avg / ((double) env->population_sizes[GOMEA_index]);

  objective_var = 0.0;
  for( i = 0; i < env->population_sizes[GOMEA_index]; i++ )
    objective_var  += (env->objective_values[GOMEA_index][i]-objective_avg)*(env->objective_values[GOMEA_index][i]-objective_avg);
  objective_var  = objective_var / ((double) env->population_sizes[GOMEA_index]);

  if( objective_var <= 0.0 )
     objective_var = 0.0;

  if( objective_var <= env->fitness_variance_tolerance )
    return( 1 );

  return( 0 );
}



/*-=-=-=-=-=-=-=-=-=-=-=-=-= Section Variation -==-=-=-=-=-=-=-=-=-=-=-=-=-=*/
void generationalStepAllGOMEAsRecursiveFold( struct Environment *env, int GOMEA_index_smallest, int GOMEA_index_biggest );
void generationalStepAllGOMEAs(struct Environment *env)
{
  int GOMEA_index_smallest, GOMEA_index_biggest;

  GOMEA_index_biggest  = env->number_of_GOMEAs-1;
  GOMEA_index_smallest = 0;
  while( GOMEA_index_smallest <= GOMEA_index_biggest )
  {
    if( !env->terminated[GOMEA_index_smallest] )
      break;

    GOMEA_index_smallest++;
  }

  generationalStepAllGOMEAsRecursiveFold( env, GOMEA_index_smallest, GOMEA_index_biggest );
}

void generationalStepAllGOMEAsRecursiveFold( struct Environment *env, int GOMEA_index_smallest, int GOMEA_index_biggest )
{
  int i, GOMEA_index;

  //Do this step 4 times (value of number_of_subgenerations_per_GOMEA_factor)
  for( i = 0; i < env->number_of_subgenerations_per_GOMEA_factor-1; i++ )
  {
    //For each population between index_smalles and index_biggest (incl.), 
    for( GOMEA_index = GOMEA_index_smallest; GOMEA_index <= GOMEA_index_biggest; GOMEA_index++ )
    {
      //Check if it should be terminated
      if( !env->terminated[GOMEA_index] )
      {
        env->terminated[GOMEA_index] = checkTerminationSpecificGOMEA( env, GOMEA_index );
      }

      //If it shouldn't, then perform one generational step
      if( (!env->terminated[GOMEA_index]) && (GOMEA_index >= env->minimum_GOMEA_index) )
      {
        makeOffspringSpecificGOMEA( env, GOMEA_index );

        if ( env->halt_execution )
          return;

        selectFinalSurvivorsSpecificGOMEA( env, GOMEA_index );

        if ( env->halt_execution )
          return;

        computeAverageFitnessSpecificGOMEA( env, GOMEA_index );

        if ( env->halt_execution )
          return;

      }
    }

    //For all population indices between index_smallest and index_biggest (excl.)
    for( GOMEA_index = GOMEA_index_smallest; GOMEA_index < GOMEA_index_biggest; GOMEA_index++ )
    {
      //Recursively call this function again, but now with GOMEA_index as index_biggest. 
      //This call right here makes sure that population P{i} gets 4 generations for every generation of P{i + 1}
      generationalStepAllGOMEAsRecursiveFold( env, GOMEA_index_smallest, GOMEA_index );

      if ( env->halt_execution )
        return;

    }
  }
}

void makeOffspringSpecificGOMEA( struct Environment *env, int gomea_index )
{
  //Commented away re-encoding, which is applicable only to permutation (random keys)
  // // EDIT: RE-ENCODING OF SOLUTIONS
  // double *surrogate = (double *) Malloc( number_of_parameters*sizeof( double ) );
  // for( int p = 0; p < population_sizes[gomea_index]; p++ )
  // {
  //   for( int i = 0; i < number_of_parameters; i++ )
  //     surrogate[i] = randomRealUniform01();
  //   int *order1 = mergeSortDoublesIncreasing( populations[gomea_index][p], number_of_parameters );
  //   int *order2 = mergeSortDoublesIncreasing( surrogate, number_of_parameters );
  //   for( int i = 0; i < number_of_parameters; i++ )
  //     populations[gomea_index][p][order1[i]] = surrogate[order2[i]];
  //   free( order1 );
  //   free( order2 );
  // }
  // for( int i = 0; i < number_of_parameters; i++ )
  //   surrogate[i] = randomRealUniform01();
  // int *order1 = mergeSortDoublesIncreasing( elitist_solution, number_of_parameters );
  // int *order2 = mergeSortDoublesIncreasing( surrogate, number_of_parameters );
  // for( int i = 0; i < number_of_parameters; i++ )
  //   elitist_solution[order1[i]] = surrogate[order2[i]];
  // free( surrogate );

  //Assign solutions to selection, removed as this would increase selection pressure and decrease variation
  //selectForLearningLT(gomea_index);

  //Learn LT from selection
  learnFOSSpecificGOMEA( env, gomea_index );

  if( env->print_FOSs_contents )
  {
    printf("### FOS contents for GOMEA #%02d in generation #%03d\n", gomea_index, env->number_of_generations);
    printFOSContentsSpecificGOMEA( env, gomea_index );
    printf( "###################################################\n" );
  }

  generateAndEvaluateNewSolutionsToFillOffspringSpecificGOMEA( env, gomea_index );
}

/**
 * Returns the name of an installed problem.
 */
char *installedFOSStructureName( int index )
{
  switch( index )
  {
    case  0: return( (char *) "Univariate" );
    case  1: return( (char *) "Linkage Tree" );
    case  2: return( (char *) "Multiscale Linkage Neighbors" );
    case  3: return( (char *) "Linkage Trees and Neighbors" );
    case  4: return( (char *) "Filtered Linkage Tree" );
    case  5: return( (char *) "Filtered Multiscale Linkage Neighbors" );
    case  6: return( (char *) "Filtered Linkage Trees and Neighbors" );
  }

  return( NULL );
}

/**
 * Returns the number of FOS structures installed.
 */
int numberOfInstalledFOSStructures()
{
  static int result = -1;
  
  if( result == -1 )
  {
    result = 0;
    while( installedFOSStructureName( result ) != NULL )
      result++;
  }
  
  return( result );
}

/*
 * Performs tournament selection to create the dataset
 * upon which to perform structure learning.
 */
// void selectForLearningLT(int gomea_index)
// {
//   int i, j, a, b, c, *order, from, to, index_best, tournament_size;
//   int population_size = population_sizes[gomea_index];
//   int selection_size = population_sizes[gomea_index];

//   tournament_size = 2;

//   order = (int *) Malloc( population_size*sizeof( int ) );
//   for( i = 0; i < population_size; i++ )
//     order[i] = i;
//   for( i = 0; i < population_size; i++ )
//   {
//     a        = randomInt( population_size );
//     b        = randomInt( population_size );
//     c        = order[a];
//     order[a] = order[b];
//     order[b] = c;
//   }

//   from = 0;
//   to   = tournament_size-1;
//   for( i = 0; i < selection_size; i++ )
//   {
//     index_best = from;
//     for( j = from+1; j <= to; j++ )
//       if( betterFitness( objective_values[gomea_index][order[j%population_size]], constraint_values[gomea_index][order[j%population_size]], objective_values[gomea_index][order[index_best]], constraint_values[gomea_index][order[index_best]] ) )
//         index_best = j%population_size;

//       selections[gomea_index][i][j] = populations[gomea_index][order[index_best]][j];

//     from += tournament_size;
//     to   += tournament_size;
//     if( to >= population_size )
//     {
//       for( j = 0; j < population_size; j++ )
//       {
//         a        = randomInt( population_size );
//         b        = randomInt( population_size );
//         c        = order[a];
//         order[a] = order[b];
//         order[b] = c;
//       }
//       from = 0;
//       to   = tournament_size-1;
//     }
//   }

//   free( order );
// }

/**
 * Selects the FOS to be learned and calls the appropriate function to do
 * the learning.
 */
void learnFOSSpecificGOMEA( struct Environment *env, int gomea_index )
{
  int i;

  if( env->FOSs[gomea_index] != NULL )
  {
    for( i = 0; i < env->FOSs_length[gomea_index]; i++ )
      free( env->FOSs[gomea_index][i] );
    free( env->FOSs[gomea_index] );
    free( env->FOSs_number_of_indices[gomea_index] );
  }

  learnLTFOSSpecificGOMEA( env, gomea_index, 1, 0, NULL ); 

  // switch( env->FOSs_structure_index )
  // {
  //   case  0: learnUnivariateFOSSpecificGOMEA( gomea_index ); break;
  //   case  1: learnLTFOSSpecificGOMEA( env, gomea_index, 1, 0, NULL ); break;
  //   case  2: learnMLNFOSSpecificGOMEA( gomea_index, 1, 0, NULL ); break;
  //   case  3: learnLTNFOSSpecificGOMEA( gomea_index ); break;
  //   case  4: learnFilteredLTFOSSpecificGOMEA( gomea_index, 1 ); break;
  //   case  5: learnFilteredMLNFOSSpecificGOMEA( gomea_index, 1 ); break;
  //   case  6: learnFilteredLTNFOSSpecificGOMEA( gomea_index ); break;
  // }

  // switch( env->FOSs_structure_index )
  // {
  //   case  0: break;
  //   case  1: break;
  //   case  2: uniquifyFOSSpecificGOMEA( gomea_index ); break;
  //   case  3: uniquifyFOSSpecificGOMEA( gomea_index ); break;
  //   case  4: break;
  //   case  5: uniquifyFOSSpecificGOMEA( gomea_index ); break;
  //   case  6: uniquifyFOSSpecificGOMEA( gomea_index ); break;
  //   case  7: break;
  // }
}

/**
 * Learns a univariate FOS (randomized ordering of the singletons).
 */
// void learnUnivariateFOSSpecificGOMEA( int gomea_index )
// {
//   int i, FOSs_index, *order;

//   order                               = randomPermutation( number_of_parameters );
//   FOSs_length[gomea_index]            = number_of_parameters;
//   FOSs[gomea_index]                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
//   FOSs_number_of_indices[gomea_index] = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
//   FOSs_index                           = 0;
//   for( i = 0; i < number_of_parameters; i++ )
//   {
//     FOSs[gomea_index][FOSs_index]                   = (int *) Malloc( 1*sizeof( int ) );
//     FOSs[gomea_index][FOSs_index][0]                = order[FOSs_index];
//     FOSs_number_of_indices[gomea_index][FOSs_index] = 1;
//     FOSs_index++;
//   }

//   free( order );
// }

/**
 * Learns a linkage tree FOS by means of hierarchical clustering.
 * This implementation follows the reciprocal nearest neighbor approach.
 */
int **learnLTFOSSpecificGOMEA( struct Environment *env, int gomea_index, short compute_dependency_matrices, short compute_parent_child_relations, int *number_of_parent_child_relations )
{
  char     done;
  int      i, j, r0, r1, rswap, *indices, *order,
           FOSs_index, **mpm, *mpm_number_of_indices, mpm_length,
         **mpm_new, *mpm_new_number_of_indices, mpm_new_length,
          *NN_chain, NN_chain_length, **parent_child_relations,
           PCR_index, *FOSs_index_of_mpm_element;
  double **S_matrix, mul0, mul1;

  parent_child_relations   = NULL; /* Only needed to prevent compiler warnings. */
  PCR_index                = 0;    /* Only needed to prevent compiler warnings. */
  FOSs_index_of_mpm_element = NULL; /* Only needed to prevent compiler warnings. */
  if( compute_parent_child_relations )
  {
    *number_of_parent_child_relations = env->number_of_parameters-1;
    parent_child_relations = (int **) Malloc( (*number_of_parent_child_relations)*sizeof( int * ) );
    for( i = 0; i < (*number_of_parent_child_relations); i++ )
      parent_child_relations[i] = (int *) Malloc( 3*sizeof( int ) );
    FOSs_index_of_mpm_element = (int *) Malloc( env->number_of_parameters*sizeof( int ) );
    for( i = 0; i < env->number_of_parameters; i++ )
      FOSs_index_of_mpm_element[i] = i;
  }

  /* Compute dependency matrix using mutual information*/
  if( compute_dependency_matrices )
    computeDependencyMatrixMutualInformationSpecificGOMEA( env, gomea_index );

  /* Initialize MPM to the univariate factorization */
  order                 = randomPermutation( env, env->number_of_parameters );
  mpm                   = (int **) Malloc( env->number_of_parameters*sizeof( int * ) );
  mpm_number_of_indices = (int *) Malloc( env->number_of_parameters*sizeof( int ) );
  mpm_length            = env->number_of_parameters;
  for( i = 0; i < env->number_of_parameters; i++ )
  {
    indices                  = (int *) Malloc( 1*sizeof( int ) );
    indices[0]               = order[i];
    mpm[i]                   = indices;
    mpm_number_of_indices[i] = 1;
  }
  free( order );

  /* Initialize LT to the initial MPM */
  env->FOSs_length[gomea_index]            = env->number_of_parameters+env->number_of_parameters-1;
  env->FOSs[gomea_index]                   = (int **) Malloc( env->FOSs_length[gomea_index]*sizeof( int * ) );
  env->FOSs_number_of_indices[gomea_index] = (int *) Malloc( env->FOSs_length[gomea_index]*sizeof( int ) );
  FOSs_index                                             = 0;
  for( i = 0; i < mpm_length; i++ )
  {
    env->FOSs[gomea_index][FOSs_index]                   = mpm[i];
    env->FOSs_number_of_indices[gomea_index][FOSs_index] = mpm_number_of_indices[i];
    FOSs_index++;
  }

  /* Initialize similarity matrix */
  S_matrix = (double **) Malloc( env->number_of_parameters*sizeof( double * ) );
  for( i = 0; i < env->number_of_parameters; i++ )
    S_matrix[i] = (double *) Malloc( env->number_of_parameters*sizeof( double ) );
  for( i = 0; i < mpm_length; i++ )
    for( j = 0; j < mpm_length; j++ )
      S_matrix[i][j] = env->dependency_matrices[gomea_index][mpm[i][0]][mpm[j][0]];
  for( i = 0; i < mpm_length; i++ )
    S_matrix[i][i] = 0;

  NN_chain        = (int *) Malloc( (env->number_of_parameters+2)*sizeof( int ) );
  NN_chain_length = 0;
  done            = 0;
  while( !done )
  {
    if( NN_chain_length == 0 )
    {
      NN_chain[NN_chain_length] = randomInt( env, mpm_length );
      NN_chain_length++;
    }

    while( NN_chain_length < 3 )
    {
      NN_chain[NN_chain_length] = determineNearestNeighbour( NN_chain[NN_chain_length-1], S_matrix, mpm_number_of_indices, mpm_length );
      NN_chain_length++;
    }

    while( NN_chain[NN_chain_length-3] != NN_chain[NN_chain_length-1] )
    {
      NN_chain[NN_chain_length] = determineNearestNeighbour( NN_chain[NN_chain_length-1], S_matrix, mpm_number_of_indices, mpm_length );
      if( ((S_matrix[NN_chain[NN_chain_length-1]][NN_chain[NN_chain_length]] == S_matrix[NN_chain[NN_chain_length-1]][NN_chain[NN_chain_length-2]])) && (NN_chain[NN_chain_length] != NN_chain[NN_chain_length-2]) )
        NN_chain[NN_chain_length] = NN_chain[NN_chain_length-2];
      NN_chain_length++;
      if( NN_chain_length > env->number_of_parameters )
        break;
    }
    r0 = NN_chain[NN_chain_length-2];
    r1 = NN_chain[NN_chain_length-1];
    if( r0 > r1 )
    {
      rswap = r0;
      r0    = r1;
      r1    = rswap;
    }
    NN_chain_length -= 3;

    if( r1 < mpm_length ) /* This test is required for exceptional cases in which the nearest-neighbor ordering has changed within the chain while merging within that chain */
    {
      indices = (int *) Malloc( (mpm_number_of_indices[r0]+mpm_number_of_indices[r1])*sizeof( int ) );
  
      i = 0;
      for( j = 0; j < mpm_number_of_indices[r0]; j++ )
      {
        indices[i] = mpm[r0][j];
        i++;
      }
      for( j = 0; j < mpm_number_of_indices[r1]; j++ )
      {
        indices[i] = mpm[r1][j];
        i++;
      }
    
      if( compute_parent_child_relations )
      {
        parent_child_relations[PCR_index][0] = FOSs_index;
        parent_child_relations[PCR_index][1] = FOSs_index_of_mpm_element[r0];
        parent_child_relations[PCR_index][2] = FOSs_index_of_mpm_element[r1];
        FOSs_index_of_mpm_element[r0]         = FOSs_index;
        FOSs_index_of_mpm_element[r1]         = FOSs_index_of_mpm_element[mpm_length-1];
        PCR_index++;
      }
      env->FOSs[gomea_index][FOSs_index]                   = indices;
      env->FOSs_number_of_indices[gomea_index][FOSs_index] = mpm_number_of_indices[r0]+mpm_number_of_indices[r1];
      FOSs_index++;
  
      mul0 = ((double) mpm_number_of_indices[r0])/((double) mpm_number_of_indices[r0]+mpm_number_of_indices[r1]);
      mul1 = ((double) mpm_number_of_indices[r1])/((double) mpm_number_of_indices[r0]+mpm_number_of_indices[r1]);
      for( i = 0; i < mpm_length; i++ )
      {
        if( (i != r0) && (i != r1) )
        {
          S_matrix[i][r0] = mul0*S_matrix[i][r0] + mul1*S_matrix[i][r1];
          S_matrix[r0][i] = S_matrix[i][r0];
        }
      }
  
      mpm_new                   = (int **) Malloc( (mpm_length-1)*sizeof( int * ) );
      mpm_new_number_of_indices = (int *) Malloc( (mpm_length-1)*sizeof( int ) );
      mpm_new_length            = mpm_length-1;
      for( i = 0; i < mpm_new_length; i++ )
      {
        mpm_new[i]                   = mpm[i];
        mpm_new_number_of_indices[i] = mpm_number_of_indices[i];
      }
  
      mpm_new[r0]                   = indices;
      mpm_new_number_of_indices[r0] = mpm_number_of_indices[r0]+mpm_number_of_indices[r1];
      if( r1 < mpm_length-1 )
      {
        mpm_new[r1]                   = mpm[mpm_length-1];
        mpm_new_number_of_indices[r1] = mpm_number_of_indices[mpm_length-1];
  
        for( i = 0; i < r1; i++ )
        {
          S_matrix[i][r1] = S_matrix[i][mpm_length-1];
          S_matrix[r1][i] = S_matrix[i][r1];
        }
  
        for( j = r1+1; j < mpm_new_length; j++ )
        {
          S_matrix[r1][j] = S_matrix[j][mpm_length-1];
          S_matrix[j][r1] = S_matrix[r1][j];
        }
      }
  
      for( i = 0; i < NN_chain_length; i++ )
      {
        if( NN_chain[i] == mpm_length-1 )
        {
          NN_chain[i] = r1;
          break;
        }
      }
  
      free( mpm );
      free( mpm_number_of_indices );
      mpm                   = mpm_new;
      mpm_number_of_indices = mpm_new_number_of_indices;
      mpm_length            = mpm_new_length;
  
      if( mpm_length == 1 )
        done = 1;
    }
  }

  free( NN_chain );

  free( mpm_new );
  free( mpm_number_of_indices );

  for( i = 0; i < env->number_of_parameters; i++ )
    free( S_matrix[i] );
  free( S_matrix );

  free( FOSs_index_of_mpm_element );

  return( parent_child_relations );
}

/**
 * Estimates the cumulative probability distribution of a
 * single binary marginal (EDA legacy code).
 * Edited by Tobias to use the population as-is, instead of performing selection first.
 */
double *estimateParametersForSingleBinaryMarginal( struct Environment *env, int gomea_index, int *indices, int number_of_indices, int *factor_size )
{
  int     i, j, index, power_of_two;
  double *result;

  *factor_size = (int) pow( 2, number_of_indices );
  result       = (double *) Malloc( (*factor_size)*sizeof( double ) );

  for( i = 0; i < (*factor_size); i++ )
    result[i] = 0.0;

  for( i = 0; i < env->population_sizes[gomea_index]; i++ )
  {
    index        = 0;
    power_of_two = 1;
    for( j = number_of_indices-1; j >= 0; j-- )
    {
      index += (env->populations[gomea_index][i][indices[j]] == TRUE) ? power_of_two : 0;
      power_of_two *= 2;
    }

    result[index] += 1.0;
  }

  for( i = 0; i < (*factor_size); i++ )
    result[i] /= (double) env->population_sizes[gomea_index];

  for( i = 1; i < (*factor_size); i++ )
    result[i] += result[i-1];

  result[(*factor_size)-1] = 1.0;

  return( result );
}

/**
 * Determines nearest neighbour according to similarity values.
 */
int determineNearestNeighbour( int index, double **S_matrix, int *mpm_number_of_indices, int mpm_length )
{
  int i, result;

  result = 0;
  if( result == index )
    result++;
  for( i = 1; i < mpm_length; i++ )
  {
    if( ((S_matrix[index][i] > S_matrix[index][result]) || ((S_matrix[index][i] == S_matrix[index][result]) && (mpm_number_of_indices[i] < mpm_number_of_indices[result]))) && (i != index) )
      result = i;
  }

  return( result );
}

/**
 * Learns a multiscale linkage neighbors FOS.
 */
// int **learnMLNFOSSpecificGOMEA( int gomea_index, short compute_dependency_matrices, short compute_parent_child_relations, int *number_of_parent_child_relations )
// {
//   int    i, j, k, k2, k3, q, q2, **neighbors, ***buckets, **bucket_sizes, bucket_index, number_of_buckets,
//          PCR_index, **parent_child_relations, **parent_child_relations_temp, number_of_parent_child_relations_temp, parent_child_relations_temp_size;
//   double MI_max, MI_ratio;

//   parent_child_relations           = NULL; /* Only needed to prevent compiler warnings. */
//   parent_child_relations_temp      = NULL; /* Only needed to prevent compiler warnings. */
//   PCR_index                        = 0;    /* Only needed to prevent compiler warnings. */
//   parent_child_relations_temp_size = 0;    /* Only needed to prevent compiler warnings. */
//   number_of_buckets                = 1+sqrt( population_sizes[gomea_index] );
//   if( compute_parent_child_relations )
//   {
//     number_of_parent_child_relations_temp = 0;
//     parent_child_relations_temp_size      = number_of_parameters*number_of_buckets;
//     parent_child_relations_temp           = (int **) Malloc( (parent_child_relations_temp_size)*sizeof( int * ) );
//     for( i = 0; i < parent_child_relations_temp_size; i++ )
//       parent_child_relations_temp[i] = (int *) Malloc( 3*sizeof( int ) );
//   }

//   /* Compute dependency matrix */
//   if( compute_dependency_matrices )
//     computeDependencyMatrixSpecificGOMEA( gomea_index );

//   /* Create a random ordering of neighbors for each variable */
//   neighbors = (int **) Malloc( number_of_parameters*sizeof( int * ) );
//   for( i = 0; i < number_of_parameters; i++ )
//     neighbors[i] = randomPermutation( number_of_parameters );

//   /* Determine for each variable i a particular ordering: a bucket sort without sorted buckets */
//   buckets           = (int ***) Malloc( number_of_parameters*sizeof( int ** ) );
//   bucket_sizes      = (int **) Malloc( number_of_parameters*sizeof( int * ) );
//   for( i = 0; i < number_of_parameters; i++ )
//   {
//     buckets[i]      = NULL;
//     bucket_sizes[i] = NULL;

//     buckets[i]      = (int **) Malloc( number_of_buckets*sizeof( int * ) );
//     bucket_sizes[i] = (int *) Malloc( number_of_buckets*sizeof( int ) );

//     MI_max = 0;
//     for( j = 0; j < number_of_parameters; j++ )
//     {
//       if( neighbors[i][j] != i )
//       {
//         if( dependency_matrices[gomea_index][i][neighbors[i][j]] > MI_max )
//           MI_max = dependency_matrices[gomea_index][i][neighbors[i][j]];
//       }
//     }
//     MI_max = MI_max > 1 ? 1 : MI_max;

//     for( bucket_index = 0; bucket_index < number_of_buckets; bucket_index++ )
//       bucket_sizes[i][bucket_index] = 0;
//     for( j = 0; j < number_of_parameters; j++ )
//     {
//       if( neighbors[i][j] == i )
//         bucket_index = number_of_buckets-1;
//       else
//       {
//         MI_ratio = 1;
//         if( MI_max > 0 )
//           MI_ratio = dependency_matrices[gomea_index][i][neighbors[i][j]]/MI_max;
//         bucket_index = MI_ratio*((double) (number_of_buckets-1));
//         bucket_index = bucket_index >= (number_of_buckets-1) ? (number_of_buckets-1)-1 : bucket_index;
//         bucket_index = bucket_index <= 0 ? 0 : bucket_index;
//       }
//       bucket_sizes[i][bucket_index]++;
//     }

//     for( bucket_index = 0; bucket_index < number_of_buckets; bucket_index++ )
//     {
//       buckets[i][bucket_index] = NULL;
//       if( bucket_sizes[i][bucket_index] > 0 )
//         buckets[i][bucket_index] = (int *) Malloc( bucket_sizes[i][bucket_index]*sizeof( int ) );
//     }

//     for( bucket_index = 0; bucket_index < number_of_buckets; bucket_index++ )
//       bucket_sizes[i][bucket_index] = 0;
//     for( j = 0; j < number_of_parameters; j++ )
//     {
//       if( neighbors[i][j] == i )
//         bucket_index = number_of_buckets-1;
//       else
//       {
//         MI_ratio = 1;
//         if( MI_max > 0 )
//           MI_ratio = dependency_matrices[gomea_index][i][neighbors[i][j]]/MI_max;
//         bucket_index = MI_ratio*((double) (number_of_buckets-1));
//         bucket_index = bucket_index >= (number_of_buckets-1) ? (number_of_buckets-1)-1 : bucket_index;
//         bucket_index = bucket_index <= 0 ? 0 : bucket_index;
//       }
//       buckets[i][bucket_index][bucket_sizes[i][bucket_index]] = neighbors[i][j];
//       bucket_sizes[i][bucket_index]++;
//     }
//   }

//   FOSs_length[gomea_index] = 0;
//   for( i = 0; i < number_of_parameters; i++ )
//   {
//     if( number_of_parameters == 1 )
//       FOSs_length[gomea_index]++;
//     else
//     {
//       if( bucket_sizes[i][number_of_buckets-1] > 1 )
//         FOSs_length[gomea_index]++;
//       q  = bucket_sizes[i][number_of_buckets-1];
//       q2 = q;
//       for( j = number_of_buckets-2; j >= 0; j-- )
//       {
//         q2 += bucket_sizes[i][j];
//         if( ((q2 >= (2*q)) || (q2 == number_of_parameters)) && (q2 <= (number_of_parameters/2)) )
//         {
//           q = q2;
//           if( compute_parent_child_relations )
//           {
//             parent_child_relations_temp[PCR_index][0] = FOSs_length[gomea_index];
//             parent_child_relations_temp[PCR_index][1] = FOSs_length[gomea_index]-1;
//             parent_child_relations_temp[PCR_index][2] = -1;
//             PCR_index++;
//           }
//           FOSs_length[gomea_index]++;
//         }
//         if( q2 == number_of_parameters )
//           break;
//       }
//     }
//   }

//   /* Create the multiscale set of neighbors for each variable */
//   FOSs[gomea_index]                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
//   FOSs_number_of_indices[gomea_index] = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
  
//   FOSs_length[gomea_index] = 0;
//   for( i = 0; i < number_of_parameters; i++ )
//   {
//     if( number_of_parameters == 1 )
//     {
//       FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]] = 1;
//       FOSs[gomea_index][FOSs_length[gomea_index]]                   = (int *) Malloc( FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]]*sizeof( int ) );
//       FOSs[gomea_index][FOSs_length[gomea_index]][0]                = neighbors[i][0];
//       FOSs_length[gomea_index]++;
//     }
//     else
//     {
//       if( bucket_sizes[i][number_of_buckets-1] > 1 )
//       {
//         FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]] = bucket_sizes[i][number_of_buckets-1];
//         FOSs[gomea_index][FOSs_length[gomea_index]]                   = (int *) Malloc( FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]]*sizeof( int ) );
//         for( k = 0; k < bucket_sizes[i][number_of_buckets-1]; k++ )
//           FOSs[gomea_index][FOSs_length[gomea_index]][k] = buckets[i][number_of_buckets-1][k];
//         FOSs_length[gomea_index]++;
//       }
//       q  = bucket_sizes[i][number_of_buckets-1];
//       q2 = q;
//       for( j = number_of_buckets-2; j >= 0; j-- )
//       {
//         q2 += bucket_sizes[i][j];
//         if( ((q2 >= (2*q)) || (q2 == number_of_parameters)) && (q2 <= (number_of_parameters/2)) )
//         {
//           q = q2;
//           FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]] = q;
//           FOSs[gomea_index][FOSs_length[gomea_index]]                   = (int *) Malloc( FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]]*sizeof( int ) );
//           k                                 = 0;
//           for( k2 = number_of_buckets-1; k2 >= j; k2-- )
//           {
//             for( k3 = 0; k3 < bucket_sizes[i][k2]; k3++ )
//             {
//               FOSs[gomea_index][FOSs_length[gomea_index]][k] = buckets[i][k2][k3];
//               k++;
//             }
//           }
//           FOSs_length[gomea_index]++;
//         }
//         if( q2 == number_of_parameters )
//           break;
//       }
//     }
//   }

//   for( i = 0; i < number_of_parameters; i++ )
//   {
//     if( buckets[i] != NULL )
//     {
//       for( j = 0; j < number_of_buckets; j++ )
//       {
//         if( buckets[i][j] != NULL )
//         {
//           free( buckets[i][j] );
//         }
//       }
//       free( buckets[i] );
//       free( bucket_sizes[i] );
//     }
//   }
//   free( bucket_sizes );
//   free( buckets );

//   for( i = 0; i < number_of_parameters; i++ )
//   {
//     if( neighbors[i] != NULL )
//     {
//       free( neighbors[i] );
//     }
//   }
//   free( neighbors );

//   if( compute_parent_child_relations )
//   {
//     number_of_parent_child_relations_temp = PCR_index;
//     *number_of_parent_child_relations     = number_of_parent_child_relations_temp;

//     parent_child_relations = (int **) Malloc( (*number_of_parent_child_relations)*sizeof( int * ) );
//     for( i = 0; i < (*number_of_parent_child_relations); i++ )
//     {
//       parent_child_relations[i] = (int *) Malloc( 3*sizeof( int ) );
//       for( j = 0; j < 3; j++ )
//         parent_child_relations[i][j] = parent_child_relations_temp[i][j];
//     }
//     for( i = 0; i < parent_child_relations_temp_size; i++ )
//       free( parent_child_relations_temp[i] );
//     free( parent_child_relations_temp );
//   }
//   return( parent_child_relations );
// }

/**
 * Learns a multiscale linkage neighbors FOS.
 */
// void learnLTNFOSSpecificGOMEA( int gomea_index )
// {
//   learnLTNFOSWithOrWithoutFilteringSpecificGOMEA( gomea_index, 0 );
// }

/**
 * Learns a multiscale linkage neighbors FOS with or without filtering. The actual
 * work is done here.
 */
// void learnLTNFOSWithOrWithoutFilteringSpecificGOMEA( int gomea_index, short use_filtering )
// {
//   int i, j, **LT_FOS, **MLN_FOS, *LT_FOS_number_of_indices, *MLN_FOS_number_of_indices, LT_FOS_length, MLN_FOS_length;

//   /* Learn the LT FOS and create a backup copy */
//   if( use_filtering )
//     learnFilteredLTFOSSpecificGOMEA( gomea_index, 1 );
//   else
//     learnLTFOSSpecificGOMEA( gomea_index, 1, 0, NULL );
  
//   LT_FOS_length            = FOSs_length[gomea_index];
//   LT_FOS_number_of_indices = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
//   LT_FOS                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//   {
//     LT_FOS_number_of_indices[i] = FOSs_number_of_indices[gomea_index][i];
//     LT_FOS[i] = (int *) Malloc( FOSs_number_of_indices[gomea_index][i]*sizeof( int ) );
//     for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//       LT_FOS[i][j] = FOSs[gomea_index][i][j];
//   }

//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//     free( FOSs[gomea_index][i] );
//   free( FOSs[gomea_index] );
//   free( FOSs_number_of_indices[gomea_index] );

//   /* Learn the MLN FOS and create a backup copy */
//   if( use_filtering )
//     learnFilteredMLNFOSSpecificGOMEA( gomea_index, 0 );
//   else
//     learnMLNFOSSpecificGOMEA( gomea_index, 0, 0, NULL );

//   MLN_FOS_length            = FOSs_length[gomea_index];
//   MLN_FOS_number_of_indices = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
//   MLN_FOS                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//   {
//     MLN_FOS_number_of_indices[i] = FOSs_number_of_indices[gomea_index][i];
//     MLN_FOS[i] = (int *) Malloc( FOSs_number_of_indices[gomea_index][i]*sizeof( int ) );
//     for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//       MLN_FOS[i][j] = FOSs[gomea_index][i][j];
//   }

//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//     free( FOSs[gomea_index][i] );
//   free( FOSs[gomea_index] );
//   free( FOSs_number_of_indices[gomea_index] );

//   /* Construct the LTN FOS: join the LT FOS and the MLN FOS */
//   FOSs_length[gomea_index]            = LT_FOS_length + MLN_FOS_length; 
//   FOSs_number_of_indices[gomea_index] = (int *) Malloc( (LT_FOS_length + MLN_FOS_length)*sizeof( int ) );
//   FOSs[gomea_index]                   = (int **) Malloc( (LT_FOS_length + MLN_FOS_length)*sizeof( int * ) );
//   for( i = 0; i < LT_FOS_length; i++ )
//   {
//     FOSs_number_of_indices[gomea_index][i] = LT_FOS_number_of_indices[i];
//     FOSs[gomea_index][i] = (int *) Malloc( LT_FOS_number_of_indices[i]*sizeof( int ) );
//     for( j = 0; j < LT_FOS_number_of_indices[i]; j++ )
//       FOSs[gomea_index][i][j] = LT_FOS[i][j];
//   }
//   for( i = 0; i < MLN_FOS_length; i++ )
//   {
//     FOSs_number_of_indices[gomea_index][LT_FOS_length+i] = MLN_FOS_number_of_indices[i];
//     FOSs[gomea_index][LT_FOS_length+i] = (int *) Malloc( MLN_FOS_number_of_indices[i]*sizeof( int ) );
//     for( j = 0; j < MLN_FOS_number_of_indices[i]; j++ )
//       FOSs[gomea_index][LT_FOS_length+i][j] = MLN_FOS[i][j];
//   }

//   /* Free up backup memory */
//   for( i = 0; i < LT_FOS_length; i++ )
//     free( LT_FOS[i] );
//   free( LT_FOS_number_of_indices );
//   free( LT_FOS );
//   for( i = 0; i < MLN_FOS_length; i++ )
//     free( MLN_FOS[i] );
//   free( MLN_FOS_number_of_indices );
//   free( MLN_FOS );
// }

// void learnFilteredLTFOSSpecificGOMEA( int gomea_index, short compute_dependency_matrices )
// {
//   int i, **parent_child_relations, number_of_parent_child_relations;

//   parent_child_relations = learnLTFOSSpecificGOMEA( gomea_index, compute_dependency_matrices, 1, &number_of_parent_child_relations );

//   filterParentChildRelationsAndCreateNewFOSSpecificGOMEA( gomea_index, parent_child_relations, number_of_parent_child_relations );

//   for( i = 0; i < number_of_parent_child_relations; i++ )
//     free( parent_child_relations[i] );
//   free( parent_child_relations );
// }

// void learnFilteredMLNFOSSpecificGOMEA( int gomea_index, short compute_dependency_matrices )
// {
//   int i, **parent_child_relations, number_of_parent_child_relations;

//   parent_child_relations = learnMLNFOSSpecificGOMEA( gomea_index, compute_dependency_matrices, 1, &number_of_parent_child_relations );

//   filterParentChildRelationsAndCreateNewFOSSpecificGOMEA( gomea_index, parent_child_relations, number_of_parent_child_relations );

//   for( i = 0; i < number_of_parent_child_relations; i++ )
//     free( parent_child_relations[i] );
//   free( parent_child_relations );
// }

/**
 * Learns a multiscale linkage neighbors FOS.
 */
// void learnFilteredLTNFOSSpecificGOMEA( int gomea_index )
// {
//   learnLTNFOSWithOrWithoutFilteringSpecificGOMEA( gomea_index, 1 );
// }

// void filterParentChildRelationsAndCreateNewFOSSpecificGOMEA( int gomea_index, int **parent_child_relations, int number_of_parent_child_relations )
// {
//   char   *FOS_element_accepted, *linkage_strength_computed;
//   int     i, j, parent_index, child0_index, child1_index, **FOS, *FOS_number_of_indices, FOS_length;
//   double *linkage_strength;

//   FOS_element_accepted = (char *) Malloc( FOSs_length[gomea_index]*sizeof( char ) );
//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//     FOS_element_accepted[i] = 1;

//   linkage_strength_computed = (char *) Malloc( FOSs_length[gomea_index]*sizeof( char ) );
//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//     linkage_strength_computed[i] = 0;

//   linkage_strength = (double *) Malloc( FOSs_length[gomea_index]*sizeof( double ) );
//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//     linkage_strength[i] = 0;

//   for( i = 0; i < number_of_parent_child_relations; i++ )
//   {
//     parent_index = parent_child_relations[i][0];
//     child0_index = parent_child_relations[i][1];
//     child1_index = parent_child_relations[i][2];

//     if( !linkage_strength_computed[parent_index] )
//     {
//       linkage_strength[parent_index]          = computeLinkageStrengthSpecificGOMEA( gomea_index, FOSs[gomea_index][parent_index], FOSs_number_of_indices[gomea_index][parent_index] );
//       linkage_strength_computed[parent_index] = 1;
//     }
//     if( child0_index != -1 )
//     {
//       if( !linkage_strength_computed[child0_index] )
//       {
//         linkage_strength[child0_index]          = computeLinkageStrengthSpecificGOMEA( gomea_index, FOSs[gomea_index][child0_index], FOSs_number_of_indices[gomea_index][child0_index] );
//         linkage_strength_computed[child0_index] = 1;
//       }
//     }
//     if( child1_index != -1 )
//     {
//       if( !linkage_strength_computed[child1_index] )
//       {
//         linkage_strength[child1_index]          = computeLinkageStrengthSpecificGOMEA( gomea_index, FOSs[gomea_index][child1_index], FOSs_number_of_indices[gomea_index][child1_index] );
//         linkage_strength_computed[child1_index] = 1;
//       }
//     }

//     /* Remove each child if it has the same linkage strength as its parent */
//     if( child0_index != -1 )
//     {
//       if( ((linkage_strength[parent_index] >= (1.0-1e-100)*linkage_strength[child0_index])) && ((linkage_strength[parent_index] <= (1.0+1e-100)*linkage_strength[child0_index])) )
//         FOS_element_accepted[child0_index] = 0;
//     }
//     if( child1_index != -1 )
//     {
//       if( ((linkage_strength[parent_index] >= (1.0-1e-100)*linkage_strength[child1_index])) && ((linkage_strength[parent_index] <= (1.0+1e-100)*linkage_strength[child1_index])) )
//         FOS_element_accepted[child1_index] = 0;
//     }
//   }

//   /* Create a backup copy of the FOS */
//   FOS_length            = FOSs_length[gomea_index];
//   FOS_number_of_indices = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
//   FOS                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//   {
//     FOS_number_of_indices[i] = FOSs_number_of_indices[gomea_index][i];
//     FOS[i] = (int *) Malloc( FOSs_number_of_indices[gomea_index][i]*sizeof( int ) );
//     for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//       FOS[i][j] = FOSs[gomea_index][i][j];
//   }

//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//     free( FOSs[gomea_index][i] );
//   free( FOSs[gomea_index] );
//   free( FOSs_number_of_indices[gomea_index] );

//   /* Apply filters */
//   FOSs_length[gomea_index] = 0;
//   for( i = 0; i < FOS_length; i++ )
//     if( FOS_element_accepted[i] )
//       FOSs_length[gomea_index]++;
 
//   FOSs_number_of_indices[gomea_index] = (int *) Malloc( FOSs_length[gomea_index]*sizeof( int ) );
//   FOSs[gomea_index]                   = (int **) Malloc( FOSs_length[gomea_index]*sizeof( int * ) );
//   FOSs_length[gomea_index]            = 0;
//   for( i = 0; i < FOS_length; i++ )
//   {
//     if( FOS_element_accepted[i] )
//     {
//       FOSs_number_of_indices[gomea_index][FOSs_length[gomea_index]] = FOS_number_of_indices[i];
//       FOSs[gomea_index][FOSs_length[gomea_index]]                   = (int *) Malloc( FOS_number_of_indices[i]*sizeof( int ) );
//       for( j = 0; j < FOS_number_of_indices[i]; j++ )
//         FOSs[gomea_index][FOSs_length[gomea_index]][j] = FOS[i][j];
//       FOSs_length[gomea_index]++;
//     }
//   }

//   /* Free up backup memory */
//   for( i = 0; i < FOS_length; i++ )
//     free( FOS[i] );
//   free( FOS_number_of_indices );
//   free( FOS );

//   free( linkage_strength );
//   free( linkage_strength_computed );
//   free( FOS_element_accepted );
// }

// double computeLinkageStrengthSpecificGOMEA( int gomea_index, int *variables, int number_of_variables )
// {
//   int    i, j, n;
//   double result;

//   result = dependency_matrices[gomea_index][variables[0]][variables[0]];
//   if( number_of_variables > 1 )
//   {
//     result = 0;
//     n      = 0;
//     for( i = 0; i < number_of_variables; i++ )
//     {
//       for( j = i+1; j < number_of_variables; j++ )
//       {
//         result += dependency_matrices[gomea_index][variables[i]][variables[j]];
//         n++;
//       }
//     }
//     result /= (double) n;
//   }

//   return( result );
// }

/**
 * Compute the mutual information dependency matrices for the binary format
 */
void computeDependencyMatrixMutualInformationSpecificGOMEA( struct Environment *env, int gomea_index )
{
  int i, j, k, *indices, factor_size;
  double p, *cumulative_probabilities;
    /* Compute joint entropy matrix */
  if (env->use_random_linkage){     
    for( i = 0; i < env->number_of_parameters; i++ )
      for( j = i+1; j < env->number_of_parameters; j++ )
      {
        env->dependency_matrices[gomea_index][i][j] = randomRealUniform01(env); //env->dependency_matrices[gomea_index][i][i] + env->dependency_matrices[gomea_index][j][j] - env->dependency_matrices[gomea_index][i][j];
        env->dependency_matrices[gomea_index][j][i] = env->dependency_matrices[gomea_index][i][j];
      }
  } else {
    for( i = 0; i < env->number_of_parameters; i++ )
    {
      for( j = i+1; j < env->number_of_parameters; j++ )
      {
        indices                  = (int *) Malloc( 2*sizeof( int ) );
        indices[0]               = i;
        indices[1]               = j;
        cumulative_probabilities = estimateParametersForSingleBinaryMarginal( env, gomea_index, indices, 2, &factor_size );

        //dependency_matrices[gomea_index][j][i] = dependency_matrices[gomea_index][i][j];
        env->dependency_matrices[gomea_index][i][j] = 0.0;
        for( k = 0; k < factor_size; k++ )
        {
          if( k == 0 )
            p = cumulative_probabilities[k];
          else
            p = cumulative_probabilities[k]-cumulative_probabilities[k-1];
          if( p > 0 )
            env->dependency_matrices[gomea_index][i][j] += -p*log2(p);
        }

        env->dependency_matrices[gomea_index][j][i] = env->dependency_matrices[gomea_index][i][j];

        free( indices );
        free( cumulative_probabilities );
      }
      indices                  = (int *) Malloc( 1*sizeof( int ) );
      indices[0]               = i;
      cumulative_probabilities = estimateParametersForSingleBinaryMarginal( env, gomea_index, indices, 1, &factor_size );

      env->dependency_matrices[gomea_index][i][i] = 0.0;
      for( k = 0; k < factor_size; k++ )
      {
        if( k == 0 )
          p = cumulative_probabilities[k];
        else
          p = cumulative_probabilities[k]-cumulative_probabilities[k-1];
         if( p > 0 )
          env->dependency_matrices[gomea_index][i][i] += -p*log2(p);
      }

      free( indices );
      free( cumulative_probabilities );
    }

    /* Then transform into mutual information matrix MI(X,Y)=H(X)+H(Y)-H(X,Y) */
    for( i = 0; i < env->number_of_parameters; i++ )
      for( j = i+1; j < env->number_of_parameters; j++ )
      {
        env->dependency_matrices[gomea_index][i][j] = env->dependency_matrices[gomea_index][i][i] + env->dependency_matrices[gomea_index][j][j] - env->dependency_matrices[gomea_index][i][j];
        env->dependency_matrices[gomea_index][j][i] = env->dependency_matrices[gomea_index][i][j];
      }

  }
}

// void computeDependencyMatrixSpecificGOMEA( int gomea_index )
// {
//   int    i, j, k;
//   double p, entropy, average_distance;

//   for( i = 0; i < number_of_parameters; i++ )
//   {
//     for( j = i+1; j < number_of_parameters; j++ )
//     {
//       /* Compute entropy of probability that variable i has a smaller value than variable j, use inverted entropy */
//       p = 0.0;
//       for( k = 0; k < population_sizes[gomea_index]; k++ )
//       {
//         if( populations[gomea_index][k][i] < populations[gomea_index][k][j] )
//           p += 1.0;
//       }
//       p       /= (double) population_sizes[gomea_index];
//       entropy  = (p==0)||(p==1)?0:-(p*log2(p) + (1.0-p)*log2(1.0-p));

//       dependency_matrices[gomea_index][i][j] = 1.0-entropy;

//       /* Multiply by inverted average distance between variables */
//       average_distance = 0.0;
//       for( k = 0; k < population_sizes[gomea_index]; k++ )
//       {
//         average_distance += (populations[gomea_index][k][i] - populations[gomea_index][k][j])*(populations[gomea_index][k][i] - populations[gomea_index][k][j]);
//       }
//       average_distance /= (double) population_sizes[gomea_index];

//       dependency_matrices[gomea_index][i][j] *= 1.0-average_distance;

//       /* Create symmetric matrix */
//       dependency_matrices[gomea_index][j][i] = dependency_matrices[gomea_index][i][j];
//     }
//   }
// }

// /**
//  * Estimates the cumulative probability distribution of a
//  * single binary marginal.
//  */
// double *estimateParametersForSingleBinaryMarginalSpecificGOMEA( int gomea_index, int *indices, int number_of_indices, int *factor_size )
// {
//   int     i, j, index, power_of_two;
//   double *result;

//   *factor_size = (int) pow( 2, number_of_indices );
//   result       = (double *) Malloc( (*factor_size)*sizeof( double ) );

//   for( i = 0; i < (*factor_size); i++ )
//     result[i] = 0.0;

//   for( i = 0; i < population_sizes[gomea_index]; i++ )
//   {
//     index        = 0;
//     power_of_two = 1;
//     for( j = number_of_indices-1; j >= 0; j-- )
//     {
//       index += populations[gomea_index][i][indices[j]] ? power_of_two : 0;
//       power_of_two *= 2;
//     }

//     result[index] += 1.0;
//   }

//   for( i = 0; i < (*factor_size); i++ )
//     result[i] /= (double) population_sizes[gomea_index];

//   for( i = 1; i < (*factor_size); i++ )
//     result[i] += result[i-1];

//   result[(*factor_size)-1] = 1.0;

//   return( result );
// }

/**
 * Ensures that every FOS element is unique, i.e. there are no
 * duplicate linkage subsets.
 */
// void uniquifyFOSSpecificGOMEA( int gomea_index )
// {
//   short *FOSs_subset_is_duplicate;
//   int    i, j, k, q, *sorted, **FOSs_new, *FOSs_number_of_indices_new, FOSs_length_new;

//   /* First analyze which subsets are duplicates */
//   FOSs_subset_is_duplicate = (short *) Malloc( FOSs_length[gomea_index]*sizeof( short ) );
//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//     FOSs_subset_is_duplicate[i] = 0;

//   sorted = mergeSortIntegersDecreasing( FOSs_number_of_indices[gomea_index], FOSs_length[gomea_index] );

//   i = 0;
//   while( i < FOSs_length[gomea_index] )
//   {
//     /* Determine stretch of FOS elements that have the same length */
//     j = i+1;
//     while( (j < FOSs_length[gomea_index]) && (FOSs_number_of_indices[gomea_index][sorted[i]] == FOSs_number_of_indices[gomea_index][sorted[j]]) )
//       j++;

//     /* Check inside stretch for duplicates */
//     for( k = i; k < j-1; k++ )
//     {
//       if( FOSs_subset_is_duplicate[sorted[k]] )
//         continue;

//       for( q = k+1; q < j; q++ )
//       {
//         if( FOSs_subset_is_duplicate[sorted[q]] )
//           continue;

//         if( linkageSubsetsOfSameLengthAreDuplicates( FOSs[gomea_index][sorted[k]], FOSs[gomea_index][sorted[q]], FOSs_number_of_indices[gomea_index][sorted[k]] ) )
//           FOSs_subset_is_duplicate[sorted[q]] = 1;
//       }
//     }
//     i = j;
//   }

//   /* Then re-create the FOS without the duplicate sets */
//   FOSs_length_new = 0;
//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//   {
//     if( !FOSs_subset_is_duplicate[i] )
//       FOSs_length_new++;
//   }

//   FOSs_new                   = (int **) Malloc( FOSs_length_new*sizeof( int * ) );
//   FOSs_number_of_indices_new = (int *) Malloc( FOSs_length_new*sizeof( int ) );

//   j = 0;
//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//   {
//     if( !FOSs_subset_is_duplicate[i] )
//     {
//       FOSs_new[j] = (int *) Malloc( FOSs_number_of_indices[gomea_index][i]*sizeof( int ) );
//       for( k = 0; k < FOSs_number_of_indices[gomea_index][i]; k++ )
//         FOSs_new[j][k] = FOSs[gomea_index][i][k];

//       FOSs_number_of_indices_new[j] = FOSs_number_of_indices[gomea_index][i];

//       j++;
//     }
//   }

//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//     free( FOSs[gomea_index][i] );
//   free( FOSs[gomea_index] );
//   free( FOSs_number_of_indices[gomea_index] );
//   FOSs[gomea_index]                   = FOSs_new;
//   FOSs_number_of_indices[gomea_index] = FOSs_number_of_indices_new;
//   FOSs_length[gomea_index]            = FOSs_length_new;

//   free( sorted );
//   free( FOSs_subset_is_duplicate );
// }

// short linkageSubsetsOfSameLengthAreDuplicates( int *linkageSubset0, int *linkageSubset1, int length )
// {
//   short result, *linkageSubset0AsBitString, *linkageSubset1AsBitString;
//   int   i, *sorted0, *sorted1;

//   result = 0;
//   if( length == 1 )
//   {
//     if( linkageSubset0[0] == linkageSubset1[0] )
//       result = 1;
//   }
//   else if( length == number_of_parameters )
//   {
//     result = 1;
//   }
//   else if( length <= (number_of_parameters/2) )
//   {
//     sorted0 = mergeSortIntegersDecreasing( linkageSubset0, length );
//     sorted1 = mergeSortIntegersDecreasing( linkageSubset1, length );
//     result = 1;
//     for( i = 0; i < length; i++ )
//     {
//       if( linkageSubset0[sorted0[i]] != linkageSubset1[sorted1[i]] )
//       {
//         result = 0;
//         break;
//       }
//     }
//     free( sorted0 );
//     free( sorted1 );
//   }
//   else
//   {
//     linkageSubset0AsBitString = (short *) Malloc( number_of_parameters*sizeof( short ) );
//     linkageSubset1AsBitString = (short *) Malloc( number_of_parameters*sizeof( short ) );
//     for( i = 0; i < number_of_parameters; i++ )
//     {
//       linkageSubset0AsBitString[i] = 0;
//       linkageSubset1AsBitString[i] = 0;
//     }
//     for( i = 0; i < length; i++ )
//     {
//       linkageSubset0AsBitString[linkageSubset0[i]] = 1;
//       linkageSubset1AsBitString[linkageSubset1[i]] = 1;
//     }
//     result = 1;
//     for( i = 0; i < number_of_parameters; i++ )
//     {
//       if( linkageSubset0AsBitString[i] != linkageSubset1AsBitString[i] )
//       {
//         result = 0;
//         break;
//       }
//     }
//     free( linkageSubset0AsBitString );
//     free( linkageSubset1AsBitString );
//   }

//   return( result );
// }

void printFOSContentsSpecificGOMEA( struct Environment *env, int gomea_index )
{
  int i, j;

  for( i = 0; i < env->FOSs_length[gomea_index]; i++ )
  {
    printf( "# [" );
    for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
    {
      printf( "%d",env->FOSs[gomea_index][i][j] );
      if( j < env->FOSs_number_of_indices[gomea_index][i]-1 )
        printf( " " );
    }
    printf( "]\n" );
  }
  fflush( stdout );
}

/**
 * Computes the two-log of x.
 */
double math_log_two = log(2.0);
double log2( double x )
{
  return( log(x) / math_log_two );
}

/**
 * Generates new solutions.
 */
void generateAndEvaluateNewSolutionsToFillOffspringSpecificGOMEA( struct Environment *env, int gomea_index )
{
  char *solution;
  int    i, j;
  double obj, con;

  for( i = 0; i < env->population_sizes[gomea_index]; i++ )
  {
    if (env->halt_execution)
      break;
    solution = generateAndEvaluateNewSolutionBinarySpecificGOMEA( env, gomea_index, i%(env->population_sizes[gomea_index]), &obj, &con );

    for( j = 0; j < env->number_of_parameters; j++ )
      env->offsprings[gomea_index][i][j] = solution[j];

    env->objective_values_offsprings[gomea_index][i] = obj;
    env->constraint_values_offsprings[gomea_index][i] = con;

    free( solution );
  }
}

/**
 * Performs Genepool Optimal Mixing (for one solution in the population).
 * This is the _binary_ version, the original version is found below
 * gomea_index is the population index in the populations list.
 * parent_index is the solution we want to apply GOM to.
 * obj is a pointer to the objective value of this solution (so we can change it and it is also reflected in the value outside this function)
 * con is a pointer to the constraint value (but we can ignore this value)
 */
char *generateAndEvaluateNewSolutionBinarySpecificGOMEA( struct Environment *env, int gomea_index, int parent_index, double *obj, double *con )
{
  char    donor_parameters_are_the_same;
  char    *result, *backup, *backup_for_ilse;
  short   solution_has_changed;
  int     i, j, j_ilse_best, index, parameter_index_for_ilse[1];
  double  obj_backup, con_backup, obj_backup_for_ilse, con_backup_for_ilse, obj_ilse_best, con_ilse_best;

  solution_has_changed = 0;

  //Copy current solution bits to result
  result = (char *) Malloc( env->number_of_parameters*sizeof( char ) );
  for( i = 0; i < env->number_of_parameters; i++ )
    result[i] = env->populations[gomea_index][parent_index][i];

  //Copy current objective and constraint values to obj and con
  *obj = env->objective_values[gomea_index][parent_index];
  *con = env->constraint_values[gomea_index][parent_index];

  backup = (char *) Malloc( env->number_of_parameters*sizeof( char ) );
  for( i = 0; i < env->number_of_parameters; i++ )
    backup[i] = result[i];

  backup_for_ilse = NULL; /* Only needed to prevent compiler warnings. */
  if( env->use_ilse )
  {
    backup_for_ilse = (char *) Malloc( env->number_of_parameters*sizeof( char ) );
    for( i = 0; i < env->number_of_parameters; i++ )
      backup_for_ilse[i] = result[i];
  }

  obj_backup = *obj;
  con_backup = *con;

  /* Phase 1: optimal mixing with random donors */
  //Shuffle FOSs, to access randomly
  shuffleFOSSpecificGOMEA( env, gomea_index );
  if( env->use_ilse )
    shuffleFOSSubsetsSpecificGOMEA( env, gomea_index );

  //Iterate over all FOSs 
  for( i = 0; i < env->FOSs_length[gomea_index]; i++ )
  {
    if( env->FOSs_number_of_indices[gomea_index][i] == env->number_of_parameters )
      continue;
    
    if( env->halt_execution)
      break;

    //Choose random solution index as donor
    index = randomInt( env, env->population_sizes[gomea_index] );

    /* Convert index to binary representation and set factor variables. */
    //Donate bits from donor
    for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
      result[env->FOSs[gomea_index][i][j]] = env->populations[gomea_index][index][env->FOSs[gomea_index][i][j]];

    //Commented away random rescaling (applicable only for permutation)
    // // EDIT: RANDOM RESCALING
    // if( randomRealUniform01() < 0.1 )
    // {
    //   int    number_of_intervals;
    //   double min, max, range, random_left_bracket;

    //   number_of_intervals = number_of_parameters;

    //   min = result[FOSs[gomea_index][i][0]]; 
    //   max = result[FOSs[gomea_index][i][0]]; 
    //   for( j = 1; j < FOSs_number_of_indices[gomea_index][i]; j++ )
    //   {
    //     if( result[FOSs[gomea_index][i][j]] < min )
    //       min = result[FOSs[gomea_index][i][j]];
    //     if( result[FOSs[gomea_index][i][j]] > max )
    //       max = result[FOSs[gomea_index][i][j]];
    //   }
    //   range = max-min;

    //   random_left_bracket = ((double) randomInt(number_of_intervals))/((double)number_of_intervals);
    //   if( range > 0 )
    //   {
    //     for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
    //       result[FOSs[gomea_index][i][j]] = ((result[FOSs[gomea_index][i][j]]-min)/range)*(1.0/((double)number_of_intervals)) + random_left_bracket;
    //   }
    //   else
    //   {
    //     for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
    //       result[FOSs[gomea_index][i][j]] = random_left_bracket;
    //   }
    // }
    
    /* Test if the change is for the better */
    donor_parameters_are_the_same = 1;
    for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
    {
      if( backup[env->FOSs[gomea_index][i][j]] != result[env->FOSs[gomea_index][i][j]] )
      {
        donor_parameters_are_the_same = 0;
        break;
      }
    }
    //If there was a change in the parameters / bits
    if( !donor_parameters_are_the_same )
    {
      if( !env->use_ilse )
        problemEvaluation( env, env->problem_index, result, obj, con, env->FOSs_number_of_indices[gomea_index][i], env->FOSs[gomea_index][i], backup, obj_backup, con_backup, gomea_index );
      else
      {
        for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
          result[env->FOSs[gomea_index][i][j]] = backup[env->FOSs[gomea_index][i][j]];

        j_ilse_best         = 0;
        obj_ilse_best       = *obj;
        con_ilse_best       = *con;
        obj_backup_for_ilse = *obj;
        con_backup_for_ilse = *con;
        for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
        {
          if( result[env->FOSs[gomea_index][i][j]] != env->populations[gomea_index][index][env->FOSs[gomea_index][i][j]] )
          {
            result[env->FOSs[gomea_index][i][j]] = env->populations[gomea_index][index][env->FOSs[gomea_index][i][j]];
            parameter_index_for_ilse[0] = env->FOSs[gomea_index][i][j];
            problemEvaluation( env, env->problem_index, result, obj, con, 1, parameter_index_for_ilse, backup_for_ilse, obj_backup_for_ilse, con_backup_for_ilse, gomea_index );
          }
          if( (j == 0) || betterFitness( *obj, *con, obj_ilse_best, con_ilse_best ) || equalFitness( *obj, *con, obj_ilse_best, con_ilse_best ) )
          {
            j_ilse_best   = j;
            obj_ilse_best = *obj;
            con_ilse_best = *con;
          }
          backup_for_ilse[env->FOSs[gomea_index][i][j]] = env->populations[gomea_index][index][env->FOSs[gomea_index][i][j]];
          obj_backup_for_ilse = *obj;
          con_backup_for_ilse = *con;
        }
        for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
        {
          result[env->FOSs[gomea_index][i][j]]          = backup[env->FOSs[gomea_index][i][j]];
          backup_for_ilse[env->FOSs[gomea_index][i][j]] = backup[env->FOSs[gomea_index][i][j]];
        }
        for( j = 0; j <= j_ilse_best; j++ )
          result[env->FOSs[gomea_index][i][j]] = env->populations[gomea_index][index][env->FOSs[gomea_index][i][j]];
        *obj = obj_ilse_best;
        *con = con_ilse_best;
      }

      //If the new solution has equal or higher fitness, keep using it and mark it as changed
      if( betterFitness( *obj, *con, obj_backup, con_backup ) || equalFitness( *obj, *con, obj_backup, con_backup ) )
      {
        for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
          backup[env->FOSs[gomea_index][i][j]] = result[env->FOSs[gomea_index][i][j]];
        if( env->use_ilse )
        {
          for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
            backup_for_ilse[env->FOSs[gomea_index][i][j]] = result[env->FOSs[gomea_index][i][j]];
        }

        obj_backup = *obj;
        con_backup = *con;

        solution_has_changed = 1;
      }
      //Otherwise, reset solution to backup
      else
      {
        for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
          result[env->FOSs[gomea_index][i][j]] = backup[env->FOSs[gomea_index][i][j]];

        *obj = obj_backup;
        *con = con_backup;
      }
    }
  }

  /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
  // Peter's NIS threshold:
  // if( (!solution_has_changed) || (no_improvement_stretchs[gomea_index] > (10+10*(log(population_sizes[gomea_index])/log(10)))) ) 
  // Original NIS threshold (in paper)
  if( (!env->halt_execution) && ( (!solution_has_changed) || (env->no_improvement_stretchs[gomea_index] > (1+(log(env->population_sizes[gomea_index])/log(10))))) )
  {
    shuffleFOSSpecificGOMEA( env, gomea_index );
    if( env->use_ilse )
      shuffleFOSSubsetsSpecificGOMEA( env, gomea_index );

    solution_has_changed = 0;
    //Iterate over all FOSs
    for( i = 0; i < env->FOSs_length[gomea_index]; i++ )
    {
      /* Convert elite solution to binary representation and set factor variables. */
      //Donate bits from elitist donor
      for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
        result[env->FOSs[gomea_index][i][j]] = env->elitist_solution[env->FOSs[gomea_index][i][j]];

      /* Test if the change is for the better */
      donor_parameters_are_the_same = 1;
      for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
      {
        if( backup[env->FOSs[gomea_index][i][j]] != result[env->FOSs[gomea_index][i][j]] )
        {
          donor_parameters_are_the_same = 0;
          break;
        }
      }
      //If there was a change in the parameters / bits
      if( !donor_parameters_are_the_same )
      {
        if( !env->use_ilse )
          problemEvaluation( env, env->problem_index, result, obj, con, env->FOSs_number_of_indices[gomea_index][i], env->FOSs[gomea_index][i], backup, obj_backup, con_backup, gomea_index );
        else
        {
          for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
            result[env->FOSs[gomea_index][i][j]] = backup[env->FOSs[gomea_index][i][j]];

          j_ilse_best         = 0;
          obj_ilse_best       = *obj;
          con_ilse_best       = *con;
          obj_backup_for_ilse = *obj;
          con_backup_for_ilse = *con;
          for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
          {
            if( result[env->FOSs[gomea_index][i][j]] != env->elitist_solution[env->FOSs[gomea_index][i][j]] )
            {
              result[env->FOSs[gomea_index][i][j]] = env->elitist_solution[env->FOSs[gomea_index][i][j]];
              parameter_index_for_ilse[0] = env->FOSs[gomea_index][i][j];
              problemEvaluation( env, env->problem_index, result, obj, con, 1, parameter_index_for_ilse, backup_for_ilse, obj_backup_for_ilse, con_backup_for_ilse, gomea_index);
            }
            if( (j == 0) || betterFitness( *obj, *con, obj_ilse_best, con_ilse_best ) || equalFitness( *obj, *con, obj_ilse_best, con_ilse_best ) )
            {
              j_ilse_best   = j;
              obj_ilse_best = *obj;
              con_ilse_best = *con;
            }
            backup_for_ilse[env->FOSs[gomea_index][i][j]] = env->elitist_solution[env->FOSs[gomea_index][i][j]];
            obj_backup_for_ilse = *obj;
            con_backup_for_ilse = *con;
          }
          for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
          {
            result[env->FOSs[gomea_index][i][j]]          = backup[env->FOSs[gomea_index][i][j]];
            backup_for_ilse[env->FOSs[gomea_index][i][j]] = backup[env->FOSs[gomea_index][i][j]];
          }
          for( j = 0; j <= j_ilse_best; j++ )
            result[env->FOSs[gomea_index][i][j]] = env->elitist_solution[env->FOSs[gomea_index][i][j]];
          *obj = obj_ilse_best;
          *con = con_ilse_best;
        }

        //If the new solution has higher fitness, mark it as changed
        if( betterFitness( *obj, *con, obj_backup, con_backup ) )
        {
          for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
            backup[env->FOSs[gomea_index][i][j]] = result[env->FOSs[gomea_index][i][j]];
          if( env->use_ilse )
          {
            for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
              backup_for_ilse[env->FOSs[gomea_index][i][j]] = result[env->FOSs[gomea_index][i][j]];
          }

          obj_backup = *obj;
          con_backup = *con;

          solution_has_changed = 1;
        }
        //Otherwise, reset solution to backup
        else
        {
          for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
            result[env->FOSs[gomea_index][i][j]] = backup[env->FOSs[gomea_index][i][j]];

          *obj = obj_backup;
          *con = con_backup;
        }
      }
      //If the solution has been changed, stop Forced Improvement
      if( solution_has_changed )
        break;
    }

    //If the mixing with the elitist solution did not lead to an improvement, copy elitist solution into current solution.
    if( !solution_has_changed )
    {
      if( betterFitness( env->elitist_solution_objective_value, env->elitist_solution_constraint_value, *obj, *con ) )
        solution_has_changed = 1;

      for( i = 0; i < env->number_of_parameters; i++ )
        result[i] = env->elitist_solution[i];
      *obj = env->elitist_solution_objective_value;
      *con = env->elitist_solution_constraint_value;
    }
  }

  free( backup );
  if( env->use_ilse )
    free( backup_for_ilse );

  return( result );
}

// /**
//  * Performs Genepool Optimal Mixing (for one solution in the population).
//  */
// double *generateAndEvaluateNewSolutionSpecificGOMEA( int gomea_index, int parent_index, double *obj, double *con )
// {
//   char    donor_parameters_are_the_same;
//   double *result, *backup, *backup_for_ilse;
//   short   solution_has_changed;
//   int     i, j, j_ilse_best, index, parameter_index_for_ilse[1];
//   double  obj_backup, con_backup, obj_backup_for_ilse, con_backup_for_ilse, obj_ilse_best, con_ilse_best;

//   solution_has_changed = 0;

//   result = (double *) Malloc( number_of_parameters*sizeof( double ) );
//   for( i = 0; i < number_of_parameters; i++ )
//     result[i] = populations[gomea_index][parent_index][i];

//   *obj = objective_values[gomea_index][parent_index];
//   *con = constraint_values[gomea_index][parent_index];

//   backup = (double *) Malloc( number_of_parameters*sizeof( double ) );
//   for( i = 0; i < number_of_parameters; i++ )
//     backup[i] = result[i];

//   backup_for_ilse = NULL; /* Only needed to prevent compiler warnings. */
//   if( use_ilse )
//   {
//     backup_for_ilse = (double *) Malloc( number_of_parameters*sizeof( double ) );
//     for( i = 0; i < number_of_parameters; i++ )
//       backup_for_ilse[i] = result[i];
//   }

//   obj_backup = *obj;
//   con_backup = *con;

//   /* Phase 1: optimal mixing with random donors */
//   shuffleFOSSpecificGOMEA( gomea_index );
//   if( use_ilse )
//     shuffleFOSSubsetsSpecificGOMEA( gomea_index );

//   for( i = 0; i < FOSs_length[gomea_index]; i++ )
//   {
//     if( FOSs_number_of_indices[gomea_index][i] == number_of_parameters )
//       continue;

//     index = randomInt( population_sizes[gomea_index] );

//     /* Convert index to binary representation and set factor variables. */
//     for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//       result[FOSs[gomea_index][i][j]] = populations[gomea_index][index][FOSs[gomea_index][i][j]];

//     // EDIT: RANDOM RESCALING
//     if( randomRealUniform01() < 0.1 )
//     {
//       int    number_of_intervals;
//       double min, max, range, random_left_bracket;

//       number_of_intervals = number_of_parameters;

//       min = result[FOSs[gomea_index][i][0]]; 
//       max = result[FOSs[gomea_index][i][0]]; 
//       for( j = 1; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//       {
//         if( result[FOSs[gomea_index][i][j]] < min )
//           min = result[FOSs[gomea_index][i][j]];
//         if( result[FOSs[gomea_index][i][j]] > max )
//           max = result[FOSs[gomea_index][i][j]];
//       }
//       range = max-min;

//       random_left_bracket = ((double) randomInt(number_of_intervals))/((double)number_of_intervals);
//       if( range > 0 )
//       {
//         for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//           result[FOSs[gomea_index][i][j]] = ((result[FOSs[gomea_index][i][j]]-min)/range)*(1.0/((double)number_of_intervals)) + random_left_bracket;
//       }
//       else
//       {
//         for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//           result[FOSs[gomea_index][i][j]] = random_left_bracket;
//       }
//     }
    
//     /* Test if the change is for the better */
//     donor_parameters_are_the_same = 1;
//     for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//     {
//       if( backup[FOSs[gomea_index][i][j]] != result[FOSs[gomea_index][i][j]] )
//       {
//         donor_parameters_are_the_same = 0;
//         break;
//       }
//     }
//     if( !donor_parameters_are_the_same )
//     {
//       if( !use_ilse )
//         installedProblemEvaluation( problem_index, result, obj, con, FOSs_number_of_indices[gomea_index][i], FOSs[gomea_index][i], backup, obj_backup, con_backup );
//       else
//       {
//         for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//           result[FOSs[gomea_index][i][j]] = backup[FOSs[gomea_index][i][j]];

//         j_ilse_best         = 0;
//         obj_ilse_best       = *obj;
//         con_ilse_best       = *con;
//         obj_backup_for_ilse = *obj;
//         con_backup_for_ilse = *con;
//         for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//         {
//           if( result[FOSs[gomea_index][i][j]] != populations[gomea_index][index][FOSs[gomea_index][i][j]] )
//           {
//             result[FOSs[gomea_index][i][j]] = populations[gomea_index][index][FOSs[gomea_index][i][j]];
//             parameter_index_for_ilse[0] = FOSs[gomea_index][i][j];
//             installedProblemEvaluation( problem_index, result, obj, con, 1, parameter_index_for_ilse, backup_for_ilse, obj_backup_for_ilse, con_backup_for_ilse );
//           }
//           if( (j == 0) || betterFitness( *obj, *con, obj_ilse_best, con_ilse_best ) || equalFitness( *obj, *con, obj_ilse_best, con_ilse_best ) )
//           {
//             j_ilse_best   = j;
//             obj_ilse_best = *obj;
//             con_ilse_best = *con;
//           }
//           backup_for_ilse[FOSs[gomea_index][i][j]] = populations[gomea_index][index][FOSs[gomea_index][i][j]];
//           obj_backup_for_ilse = *obj;
//           con_backup_for_ilse = *con;
//         }
//         for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//         {
//           result[FOSs[gomea_index][i][j]]          = backup[FOSs[gomea_index][i][j]];
//           backup_for_ilse[FOSs[gomea_index][i][j]] = backup[FOSs[gomea_index][i][j]];
//         }
//         for( j = 0; j <= j_ilse_best; j++ )
//           result[FOSs[gomea_index][i][j]] = populations[gomea_index][index][FOSs[gomea_index][i][j]];
//         *obj = obj_ilse_best;
//         *con = con_ilse_best;
//       }

//       if( betterFitness( *obj, *con, obj_backup, con_backup ) || equalFitness( *obj, *con, obj_backup, con_backup ) )
//       {
//         for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//           backup[FOSs[gomea_index][i][j]] = result[FOSs[gomea_index][i][j]];
//         if( use_ilse )
//         {
//           for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//             backup_for_ilse[FOSs[gomea_index][i][j]] = result[FOSs[gomea_index][i][j]];
//         }

//         obj_backup = *obj;
//         con_backup = *con;

//         solution_has_changed = 1;
//       }
//       else
//       {
//         for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//           result[FOSs[gomea_index][i][j]] = backup[FOSs[gomea_index][i][j]];

//         *obj = obj_backup;
//         *con = con_backup;
//       }
//     }
//   }

//   /* Phase 2 (Forced Improvement): optimal mixing with elitist solution */
//   if( (!solution_has_changed) || (no_improvement_stretchs[gomea_index] > (10+10*(log(population_sizes[gomea_index])/log(10)))) )
//   {
//     shuffleFOSSpecificGOMEA( gomea_index );
//     if( use_ilse )
//       shuffleFOSSubsetsSpecificGOMEA( gomea_index );

//     solution_has_changed = 0;
//     for( i = 0; i < FOSs_length[gomea_index]; i++ )
//     {
//       /* Convert elite solution to binary representation and set factor variables. */
//       for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//         result[FOSs[gomea_index][i][j]] = elitist_solution[FOSs[gomea_index][i][j]];

//       /* Test if the change is for the better */
//       donor_parameters_are_the_same = 1;
//       for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//       {
//         if( backup[FOSs[gomea_index][i][j]] != result[FOSs[gomea_index][i][j]] )
//         {
//           donor_parameters_are_the_same = 0;
//           break;
//         }
//       }
//       if( !donor_parameters_are_the_same )
//       {
//         if( !use_ilse )
//           installedProblemEvaluation( problem_index, result, obj, con, FOSs_number_of_indices[gomea_index][i], FOSs[gomea_index][i], backup, obj_backup, con_backup );
//         else
//         {
//           for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//             result[FOSs[gomea_index][i][j]] = backup[FOSs[gomea_index][i][j]];

//           j_ilse_best         = 0;
//           obj_ilse_best       = *obj;
//           con_ilse_best       = *con;
//           obj_backup_for_ilse = *obj;
//           con_backup_for_ilse = *con;
//           for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//           {
//             if( result[FOSs[gomea_index][i][j]] != elitist_solution[FOSs[gomea_index][i][j]] )
//             {
//               result[FOSs[gomea_index][i][j]] = elitist_solution[FOSs[gomea_index][i][j]];
//               parameter_index_for_ilse[0] = FOSs[gomea_index][i][j];
//               installedProblemEvaluation( problem_index, result, obj, con, 1, parameter_index_for_ilse, backup_for_ilse, obj_backup_for_ilse, con_backup_for_ilse );
//             }
//             if( (j == 0) || betterFitness( *obj, *con, obj_ilse_best, con_ilse_best ) || equalFitness( *obj, *con, obj_ilse_best, con_ilse_best ) )
//             {
//               j_ilse_best   = j;
//               obj_ilse_best = *obj;
//               con_ilse_best = *con;
//             }
//             backup_for_ilse[FOSs[gomea_index][i][j]] = elitist_solution[FOSs[gomea_index][i][j]];
//             obj_backup_for_ilse = *obj;
//             con_backup_for_ilse = *con;
//           }
//           for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//           {
//             result[FOSs[gomea_index][i][j]]          = backup[FOSs[gomea_index][i][j]];
//             backup_for_ilse[FOSs[gomea_index][i][j]] = backup[FOSs[gomea_index][i][j]];
//           }
//           for( j = 0; j <= j_ilse_best; j++ )
//             result[FOSs[gomea_index][i][j]] = elitist_solution[FOSs[gomea_index][i][j]];
//           *obj = obj_ilse_best;
//           *con = con_ilse_best;
//         }

//         if( betterFitness( *obj, *con, obj_backup, con_backup ) )
//         {
//           for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//             backup[FOSs[gomea_index][i][j]] = result[FOSs[gomea_index][i][j]];
//           if( use_ilse )
//           {
//             for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//               backup_for_ilse[FOSs[gomea_index][i][j]] = result[FOSs[gomea_index][i][j]];
//           }

//           obj_backup = *obj;
//           con_backup = *con;

//           solution_has_changed = 1;
//         }
//         else
//         {
//           for( j = 0; j < FOSs_number_of_indices[gomea_index][i]; j++ )
//             result[FOSs[gomea_index][i][j]] = backup[FOSs[gomea_index][i][j]];

//           *obj = obj_backup;
//           *con = con_backup;
//         }
//       }
//       if( solution_has_changed )
//         break;
//     }

//     if( !solution_has_changed )
//     {
//       if( betterFitness( elitist_solution_objective_value, elitist_solution_constraint_value, *obj, *con ) )
//         solution_has_changed = 1;

//       for( i = 0; i < number_of_parameters; i++ )
//         result[i] = elitist_solution[i];
//       *obj = elitist_solution_objective_value;
//       *con = elitist_solution_constraint_value;
//     }
//   }

//   free( backup );
//   if( use_ilse )
//     free( backup_for_ilse );

//   return( result );
// }

/**
 * Shuffles the FOS (ordering), but not the contents
 * of the linkage subsets themselves.
 */
void shuffleFOSSpecificGOMEA( struct Environment *env, int gomea_index )
{
  int i, *order, **FOSs_new, *FOSs_number_of_indices_new;

  FOSs_new                   = (int **) Malloc( env->FOSs_length[gomea_index]*sizeof( int * ) );
  FOSs_number_of_indices_new = (int *) Malloc( env->FOSs_length[gomea_index]*sizeof( int ) );
  order                     = randomPermutation( env, env->FOSs_length[gomea_index] );
  for( i = 0; i < env->FOSs_length[gomea_index]; i++ )
  {
    FOSs_new[i]                   = env->FOSs[gomea_index][order[i]];
    FOSs_number_of_indices_new[i] = env->FOSs_number_of_indices[gomea_index][order[i]];
  }
  free( env->FOSs[gomea_index] );
  free( env->FOSs_number_of_indices[gomea_index] );
  env->FOSs[gomea_index]                   = FOSs_new;
  env->FOSs_number_of_indices[gomea_index] = FOSs_number_of_indices_new;

  free( order );
}

/**
 * Shuffles the linkage subsets (ordering) in the FOS, but not the FOS itself.
 */
void shuffleFOSSubsetsSpecificGOMEA( struct Environment *env, int gomea_index )
{
  int i, j, *order, *FOSs_subset_new;

  for( i = 0; i < env->FOSs_length[gomea_index]; i++ )
  {
    order = randomPermutation( env, env->FOSs_number_of_indices[gomea_index][i] );

    FOSs_subset_new = (int *) Malloc( env->FOSs_number_of_indices[gomea_index][i]*sizeof( int ) );
    for( j = 0; j < env->FOSs_number_of_indices[gomea_index][i]; j++ )
      FOSs_subset_new[j] = env->FOSs[gomea_index][i][order[j]];
    free( env->FOSs[gomea_index][i] );
    env->FOSs[gomea_index][i] = FOSs_subset_new;

    free( order );
  }
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=- Section Ezilaitini -=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Undoes GOMEA initializations.
 */
void ezilaitiniAllGOMEAs( struct Environment *env )
{
  int i;

  for( i = 0; i < env->number_of_GOMEAs; i++ )
    ezilaitiniSpecificGOMEA( env, i );

  free( env->FOSs_length );
  free( env->FOSs_number_of_indices );
  free( env->FOSs );
  free( env->dependency_matrices );
  free( env->populations );
  free( env->objective_values );
  free( env->constraint_values );
  free( env->average_objective_values );
  free( env->average_constraint_values );
  free( env->terminated );
  free( env->offsprings );
  free( env->objective_values_offsprings );
  free( env->constraint_values_offsprings );
  free( env->objective_values_best_of_generation );
  free( env->constraint_values_best_of_generation );
}

void ezilaitiniSpecificGOMEA( struct Environment *env, int gomea_index )
{
  int i;

  if( env->FOSs[gomea_index] != NULL )
  {
    for( i = 0; i < env->FOSs_length[gomea_index]; i++ )
      free( env->FOSs[gomea_index][i] );
    free( env->FOSs[gomea_index] );
    free( env->FOSs_number_of_indices[gomea_index] );
  }

  for( i = 0; i < env->number_of_parameters; i++ )
    free( env->dependency_matrices[gomea_index][i] );
  free( env->dependency_matrices[gomea_index] );

  ezilaitiniSpecificGOMEAMemoryForPopulationAndOffspring( env, gomea_index );
}

/**
 * Initializes the memory for a single GOMEA thread, for the population only.
 */
void ezilaitiniSpecificGOMEAMemoryForPopulationAndOffspring( struct Environment *env, int gomea_index )
{
  int i;

  for( i = 0; i < env->population_sizes[gomea_index]; i++ )
    free( env->offsprings[gomea_index][i] );

  for( i = 0; i < env->population_sizes[gomea_index]; i++ )
    free( env->populations[gomea_index][i] );

  free( env->populations[gomea_index] );
  free( env->objective_values[gomea_index] );
  free( env->constraint_values[gomea_index] );
  free( env->offsprings[gomea_index] );
  free( env->objective_values_offsprings[gomea_index] );
  free( env->constraint_values_offsprings[gomea_index] );
}

void ezilaitiniValueAndSetOfSolutionsToReach(struct Environment *env)
{
  int i;

  for( i = 0; i < env->number_of_solutions_in_sostr; i++ )
    free( env->sostr[i] );
  free( env->sostr );
}

// void ezilaitiniProblem( int index )
// {
//   int j;
//   switch( index )
//   {
//     case  0: break;
//     default:
//       for (j = 0; j < njobs; j++)
//         free(ptime[j]);
//   }
// }
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Time -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
long getMilliSecondsRunning(struct Environment *env)
{
  return( getMilliSecondsRunningSinceTimeStamp( env->timestamp_start ) );
}

long getMilliSecondsRunningAfterInit(struct Environment *env)
{
  return( getMilliSecondsRunningSinceTimeStamp( env->timestamp_start_after_init ) );
}

long getMilliSecondsRunningSinceTimeStamp( long timestamp )
{
  long timestamp_now, difference;

  timestamp_now = getCurrentTimeStampInMilliSeconds();

  difference = timestamp_now-timestamp;

  return( difference );
}

long getCurrentTimeStampInMilliSeconds()
{
  struct timeval tv;
  long   result;

  gettimeofday( &tv, NULL );
  result = (tv.tv_sec * 1000) + (tv.tv_usec / 1000);

  return( result );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Run -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
/**
 * Initializes the VOSOSTR, the random number generator and the problem and runs the GOMEA.
 */
void run( struct Environment *env )
{
  env->timestamp_start = getCurrentTimeStampInMilliSeconds();

  // initializeValueAndSetOfSolutionsToReach( env );

  initializeRandomNumberGenerator( env );

  // initializeProblem( problem_index );

  if( env->print_verbose_overview )
    printVerboseOverview( env );

  env->timestamp_start_after_init = getCurrentTimeStampInMilliSeconds();

  multiPopGOMEA( env );

  // writeRunningTime( env, (char *) "total_running_time.dat" );

  //ezilaitiniProblem( problem_index );

  ezilaitiniValueAndSetOfSolutionsToReach( env );
}

void multiPopGOMEA(struct Environment *env)
{
  env->maximum_number_of_GOMEAs                  = 25;
  env->number_of_subgenerations_per_GOMEA_factor = 4;
  env->base_population_size                      = 1;

  env->number_of_GOMEAs               = 0;
  env->number_of_generations          = 0;
  env->number_of_evaluations          = 0;
  env->number_of_bit_flip_evaluations = 0;
  env->minimum_GOMEA_index            = 0;

  env->population_sizes                   = (int *) Malloc( env->maximum_number_of_GOMEAs*sizeof( int ) );
  env->number_of_subgenerations_per_GOMEA = (int *) Malloc( env->maximum_number_of_GOMEAs*sizeof( int ) );
  env->no_improvement_stretchs            = (int *) Malloc( env->maximum_number_of_GOMEAs*sizeof( int ) );
  env->elitist_solution                   = (char *) Malloc( env->number_of_parameters*sizeof( char ) );
  env->ltga_result.elitist_solution               = (char *) Malloc( env->number_of_parameters*sizeof( char ) );

  while( !checkTermination( env ) )
  {
    if( env->number_of_GOMEAs < env->maximum_number_of_GOMEAs )
      initializeNewGOMEA( env );

    // if( env->write_generational_statistics )
    //   writeGenerationalStatistics( env );

    // if( env->write_generational_solutions )
    //   writeGenerationalSolutions( env, 0 );

    generationalStepAllGOMEAs( env );

    env->number_of_generations++;
  }

  // if( env->write_generational_statistics )
  //   writeGenerationalStatistics( env );

  // if( env->write_generational_solutions )
  //   writeGenerationalSolutions( env, 1 );

  ezilaitiniAllGOMEAs( env );

  free( env->elitist_solution );
  free( env->no_improvement_stretchs );
  free( env->number_of_subgenerations_per_GOMEA );
  free( env->population_sizes );
}
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/





/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=- Section Main -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * The main function:
 * - interpret parameters on the command line
 * - run the algorithm with the interpreted parameters
 */
// int main( int argc, char **argv )
// {
//   interpretCommandLine( argc, argv );

//   run();

//   return( 0 );
// }
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

struct Environment getNewEnvironment() {
  struct Environment env;

  env.terminated = 0;

  env.elitist_solution = 0;
  env.sostr = 0;
  env.populations = 0;
  env.offsprings = 0;

  env.halt_execution = 0;
  env.write_generational_statistics = 0;
  env.write_generational_solutions = 0;
  env.print_verbose_overview = 0;
  env.print_FOSs_contents = 0;
  env.use_ilse = 0;
  env.use_random_linkage = 0;
  env.use_fitness_variance_tolerance = 0;
  env.use_premature_stopping = 0;
  env.vosostr_hit_status = 0;
  env.vtr_exists = 0;
  env.sostr_exists = 0;

  env.problem_index = 0;
  env.FOSs_structure_index = 0;
  env.number_of_parameters = 0;
  env.number_of_solutions_in_sostr = 0;
  env.number_of_generations = 0;
  env.population_sizes = 0;
  env.base_population_size = 0;
  env.number_of_subgenerations_per_GOMEA = 0;
  env.number_of_subgenerations_per_GOMEA_factor = 0;
  env.no_improvement_stretchs = 0;
  env.number_of_GOMEAs = 0;
  env.maximum_number_of_GOMEAs = 0;
  env.FOSs = 0;
  env.FOSs_number_of_indices = 0;
  env.FOSs_length = 0;
  env.minimum_GOMEA_index = 0;

  env.maximum_number_of_evaluations = 0;
  env.maximum_number_of_milliseconds = 0;
  env.timestamp_start = 0;
  env.timestamp_start_after_init = 0;
  env.number_of_evaluations = 0;
  env.elitist_solution_number_of_evaluations = 0;
  env.elitist_solution_hitting_time = 0;
  env.vosostr_number_of_evaluations = 0;
  env.vosostr_hitting_time = 0;

  env.number_of_bit_flip_evaluations = 0;
  env.elitist_solution_number_of_bit_flip_evaluations = 0;
  env.vosostr_number_of_bit_flip_evaluations = 0;

  env.elitist_solution_objective_value = 0;
  env.elitist_solution_constraint_value = 0;
  env.vtr = 0;
  env.objective_values = 0;
  env.constraint_values = 0;
  env.objective_values_offsprings = 0;
  env.constraint_values_offsprings = 0;
  env.objective_values_best_of_generation = 0;
  env.constraint_values_best_of_generation = 0;
  env.average_objective_values = 0;
  env.average_constraint_values = 0;
  env.dependency_matrices = 0;

  env.random_seed = 0;
  env.random_seed_changing = 0;


  env.input_parameters.M = 0;
  env.input_parameters.k = 0;
  env.input_parameters.o = 0;
  env.input_parameters.b = 0;
  env.cliques = 0;
  env.codomain_values = 0;
  env.global_optima_score = 0;

  env.ltga_result.elitist_solution = 0;
  env.ltga_result.elitist_fitness = 0;
  env.ltga_result.elitist_evaluations = 0;
  env.ltga_result.elitist_time = 0;
  env.ltga_result.global_optima_found = 0;
  env.ltga_result.pop_size_glob_opt = 0;

  env.fitness_variance_tolerance = 0;

  return env;
}

int setParameters(struct Environment *env, struct LTGAParameters *run_parameters, uint **cliques_memory, double **codomain_memory ) {
  env->write_generational_statistics    = run_parameters->write_generational_statistics; //Voor al deze variabelen geldt dat 0 false en > 0 of 1 true is.
  env->write_generational_solutions     = run_parameters->write_generational_solutions;
  env->print_verbose_overview           = run_parameters->print_verbose_overview;
  env->print_FOSs_contents              = run_parameters->print_FOSs_contents;
  env->use_ilse                         = run_parameters->use_ilse;
  env->use_random_linkage               = run_parameters->use_random_linkage;
  env->vtr_exists                       = run_parameters->use_value_to_reach;
  env->use_fitness_variance_tolerance   = run_parameters->use_fitness_variance_tolerance;
  env->use_premature_stopping           = run_parameters->use_premature_stopping;
  
  env->problem_index = 1; //we don't use this
  env->FOSs_structure_index = 1; //set fos to LT

  env->maximum_number_of_milliseconds = run_parameters->maximum_number_of_milliseconds;
  env->maximum_number_of_evaluations = run_parameters->maximum_number_of_evaluations;

  env->input_parameters.M = run_parameters->m;
  env->input_parameters.k = run_parameters->k;
  env->input_parameters.o = run_parameters->o;
  env->input_parameters.b = run_parameters->b;

  env->global_optima_score = run_parameters->global_optima_score;
  env->fitness_variance_tolerance = run_parameters->fitness_variance_tolerance;

  if (env->input_parameters.M == 0 || env->input_parameters.k <= env->input_parameters.o) {
    return -1;
  }
  env->number_of_parameters = (env->input_parameters.M - 1) * (env->input_parameters.k - env->input_parameters.o) + env->input_parameters.k;

  env->cliques = cliques_memory;
  env->codomain_values = codomain_memory;
  return 0;
}

struct LTGAResultC run_with_parameters(struct LTGAParameters run_parameters, uint **cliques_memory, double **codomain_memory )
{
  //interpretCommandLine( argc, argv );

  struct LTGAResultC ltga_result;
  struct Environment env = getNewEnvironment();

  int parameters_success = setParameters(&env, &run_parameters, cliques_memory, codomain_memory);

  if (parameters_success == -1 ) {
    ltga_result = env.ltga_result;
    return ltga_result;
  }

  run( &env );

  ltga_result = env.ltga_result;
  return ltga_result;
}

uint * *const get_memory_pointer_cliques(uint m) {
  uint * *const cliques_memory = (uint * *const) Malloc( m * sizeof( uint * const ));

  return cliques_memory;
}

double * *const get_memory_pointer_codomain(uint m) {
  double * *const codomain_memory = (double * *const) Malloc( m * sizeof( double * const ));

  return codomain_memory;
}

void free_memory_pointer_cliques(uint * *const cliques_memory) {
  free(cliques_memory);
}

void free_memory_pointer_codomain(double * *const codomain_memory) {
  free(codomain_memory);
}

void free_memory_pointer_elitist_solution(char * elitist_solution) {
  free(elitist_solution);
}

void print_clique_memory(uint **cliques_memory, uint m, uint k) {
  for(uint i = 0; i < m; i++) {
    printf("Clique %u: ", i);
    for (uint j = 0; j < k; j++) {
      printf("x%u, ", cliques_memory[i][j]);
    }
    printf("\n");
  }
}

void print_codomain_memory(double **codomain_memory, uint m, uint k) {
  for(uint i = 0; i < m; i++) {
    printf("Clique %u: \n", i);
    for (int j = 0; j < (1 << k); j++) {
      printf("%lf \n", codomain_memory[i][j]);
    }
  }
}
/*-=-=-=-=-=-=-=-=-=-=-=- Section Important Changes =-=-=-=-=-=-=-=-=-=-=-=-*/
/**
 * Removed Selection before learning FOS (LT), which was done in the binary version of LTGA.
 * Use Mutual Information to calculate Dependency Matrices
 * Use binary format to store solutions
 * Get clique tree and codomain values from file/function arguments
 * Calculate fitness for binary solutions from this clique tree and codomain values
 */
/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/
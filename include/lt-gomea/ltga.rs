#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

use libc;
use structopt::StructOpt;

use super::algorithm::{Algorithm, Measurement};
use problem_generator::problem::clique_tree::{CliqueTree, SolutionFit};

///Result of LTGA call, it's important to stress that when we don't use_value_to_reach,
/// then the returned global_optima_found will always be false and pop_size_glob_opt 0
pub struct LTGAResult {
    elitist_solution: Vec<u32>,
    elitist_fitness: f64,
    elitist_evaluations: u32,
    elitist_time: u64,
    global_optima_found: bool,
    pop_size_glob_opt: u32,
}

impl LTGAResult {
    pub fn from_LTGAResultC(ltga_result_c: LTGAResultC, problem_size: usize) -> Self {
        let mut solution = Vec::with_capacity(problem_size);
        unsafe {
            for i in 0..problem_size as usize {
                solution.push(*ltga_result_c.elitist_solution.offset(i as isize) as u32);
            }
        }
        LTGAResult {
            elitist_solution: solution,
            elitist_fitness: ltga_result_c.elitist_fitness,
            elitist_evaluations: ltga_result_c.elitist_evaluations as u32,
            elitist_time: ltga_result_c.elitist_time as u64,
            global_optima_found: ltga_result_c.global_optima_found != 0,
            pop_size_glob_opt: ltga_result_c.pop_size_glob_opt as u32,
        }
    }
}

#[derive(Debug, Clone)]
pub struct LTGA {
    problem_size: usize,
    ltga_parameters: LTGAParameters,
}

impl Algorithm for LTGA {
    fn init() {}
    fn run(&mut self, clique_tree: &CliqueTree, max_number_evaluations: u32) -> Measurement {
        let ltga_result = self.run_ltga(clique_tree, max_number_evaluations);

        let elitist_solutionfit = SolutionFit {
            solution: ltga_result.elitist_solution,
            fitness: ltga_result.elitist_fitness,
        };

        //If we pass the global optima value and tell it to stop when it is done, check if we agree that a global optimum has been found
        if self.ltga_parameters.use_value_to_reach != 0 {
            if ltga_result.global_optima_found
                != clique_tree.is_global_optimum(&elitist_solutionfit)
            {
                println!(
                    "LTGA result optimum: {} , clique tree optimum: {}",
                    ltga_result.elitist_fitness, clique_tree.glob_optima_score
                );
                let equal_solutions =
                    elitist_solutionfit.solution == clique_tree.glob_optima_strings[0];
                println!(
                    "LTGA result optimum and clique tree optimum solutions equal: {}",
                    equal_solutions
                );
                println!(
                    "LTGA optimum to be found: {}, clique tree optimum to be found: {}",
                    self.ltga_parameters.global_optima_score, clique_tree.glob_optima_score
                );
            }
            assert_eq!(
                ltga_result.global_optima_found,
                clique_tree.is_global_optimum(&elitist_solutionfit)
            );
        }

        Measurement {
            fitness: ltga_result.elitist_fitness,
            number_evaluations: ltga_result.elitist_evaluations,
            global_optima_found: ltga_result.global_optima_found,
            pop_size_glob_opt: ltga_result.pop_size_glob_opt,
        }
    }
}

impl LTGA {
    pub fn new(problem_size: usize, ltga_parameters: LTGAParameters) -> Self {
        LTGA {
            problem_size,
            ltga_parameters,
        }
    }
    pub fn new_default() -> Self {
        let ltga_parameters = LTGAParameters {
            write_generational_statistics: 0,
            write_generational_solutions: 0,
            print_verbose_overview: 0,
            print_FOSs_contents: 0,
            use_ilse: 0,
            use_random_linkage: 0,
            use_value_to_reach: 1,
            use_fitness_variance_tolerance: 1,
            use_premature_stopping: 0,
            maximum_number_of_evaluations: 0,
            maximum_number_of_milliseconds: 0,
            m: 0,
            k: 0,
            o: 0,
            b: 0,
            global_optima_score: 0.0,
            fitness_variance_tolerance: 0.00001,
        };
        LTGA::new(0, ltga_parameters)
    }

    pub fn get_ltga_instance_from_options(ltga_options: &Option<LTGAOptions>) -> LTGA {
        let ltga_options_clone = ltga_options.clone();
        match ltga_options_clone {
            None => LTGA::new_default(),
            Some(ltga_options) => match ltga_options {
                LTGAOptions::LTGAOptions {
                    use_ilse,
                    use_random_linkage,
                    use_value_to_reach,
                    use_premature_stopping,
                    fitness_variance_tolerance,
                    maximum_number_of_evaluations,
                    maximum_number_of_milliseconds,
                } => {
                    let ilse = if use_ilse { 1 } else { 0 };
                    let random_linkage = if use_random_linkage { 1 } else { 0 };
                    let value_to_reach = if use_value_to_reach { 1 } else { 0 };
                    let premature_stopping = if use_premature_stopping { 1 } else { 0 };

                    let fitness_variance_flag = if fitness_variance_tolerance.is_some() {
                        1
                    } else {
                        0
                    };
                    let fitness_variance_value = fitness_variance_tolerance.unwrap_or(0.00001);
                    let maximum_number_of_evaluations_value =
                        maximum_number_of_evaluations.unwrap_or_default();
                    let maximum_number_of_milliseconds_value =
                        maximum_number_of_milliseconds.unwrap_or_default();

                    let ltga_parameters = LTGAParameters::new(
                        0,
                        0,
                        0,
                        0,
                        ilse,
                        random_linkage,
                        value_to_reach,
                        fitness_variance_flag,
                        premature_stopping,
                        maximum_number_of_evaluations_value,
                        maximum_number_of_milliseconds_value as u64,
                        0,
                        0,
                        0,
                        0,
                        0.0,
                        fitness_variance_value,
                    );
                    LTGA::new(0, ltga_parameters)
                }
            },
        }
    }

    pub fn run_ltga(
        &mut self,
        clique_tree: &CliqueTree,
        maximum_number_evaluations: u32,
    ) -> LTGAResult {
        let m = clique_tree.input_parameters.m;
        let k = clique_tree.input_parameters.k;
        let o = clique_tree.input_parameters.o;
        let b = clique_tree.input_parameters.b;

        let problem_size = (clique_tree.input_parameters.m - 1)
            * (clique_tree.input_parameters.k - clique_tree.input_parameters.o)
            + clique_tree.input_parameters.k;
        self.problem_size = problem_size as usize;
        self.ltga_parameters.maximum_number_of_evaluations = maximum_number_evaluations as i64;
        self.ltga_parameters.maximum_number_of_milliseconds = maximum_number_evaluations as i64;
        self.ltga_parameters.m = m;
        self.ltga_parameters.k = k;
        self.ltga_parameters.o = o;
        self.ltga_parameters.b = b;
        self.ltga_parameters.global_optima_score = clique_tree.glob_optima_score;

        unsafe {
            let cliques_memory: *mut *const u32 = get_memory_pointer_cliques(m) as *mut *const u32;
            let codomain_memory: *mut *const f64 =
                get_memory_pointer_codomain(m) as *mut *const f64;

            fill_clique_memory(&clique_tree.cliques, cliques_memory, m);
            fill_codomain_memory(&clique_tree.codomain_values, codomain_memory, m);

            // print_clique_memory(cliques_memory, m, k);
            // print_codomain_memory(codomain_memory, m, k);
            // for i in 0..m {
            //     print!("Clique {}: ", i);
            //     for j in 0..k {
            //         let number = *(*cliques_memory.offset(i as isize)).offset(j as isize);
            //         print!("X{}, ", number );
            //     }
            //     println!("");
            // }

            // for i in 0..m {
            //     println!("Clique {}: ", i);
            //     for j in 0..(1 << k) {
            //         let codomain_value = *(*codomain_memory.offset(i as isize)).offset(j as isize);
            //         println!("v: {}, ", codomain_value );
            //     }
            // }

            //Call LTGA here!
            let ltga_result_c = run_with_parameters(
                self.ltga_parameters,
                cliques_memory as *mut *mut u32,
                codomain_memory as *mut *mut f64,
            );

            // for i in 0..problem_size as usize {
            //     ltga_result_solution.push( *ltga_result.elitist_solution.offset(i as isize) as u32);
            // }
            let ltga_result =
                LTGAResult::from_LTGAResultC(ltga_result_c, self.problem_size as usize);

            // assert!(clique_tree.glob_optima_strings.contains(&ltga_result.elitist_solution));
            // assert!(is_equal_fitness(clique_tree.glob_optima_score, ltga_result.elitist_fitness));

            free_memory_pointer_elitist_solution(ltga_result_c.elitist_solution);
            free_memory_pointer_cliques(cliques_memory as *mut *mut u32);
            free_memory_pointer_codomain(codomain_memory as *mut *mut f64);
            ltga_result
        }
    }
}

#[derive(StructOpt, Debug, Clone)]
pub enum LTGAOptions {
    #[structopt(name = "ltga-options")]
    LTGAOptions {
        ///Flag for use of ILSE
        #[structopt(short = "i")]
        use_ilse: bool,
        ///Flag for use of random linkage
        #[structopt(short = "r")]
        use_random_linkage: bool,
        ///Flag for use of value to reach (global optima score)
        #[structopt(short = "v")]
        use_value_to_reach: bool,
        ///Flag for use of premature stopping (smaller populations when lower average fitness)
        #[structopt(short = "p")]
        use_premature_stopping: bool,
        ///Parameter for the use and value of fitness variance tolerance
        #[structopt(short = "f")]
        fitness_variance_tolerance: Option<f64>,
        ///Parameter for the maximum number of evaluations an algorithm is allowed to run
        #[structopt(short = "e")]
        maximum_number_of_evaluations: Option<u32>,
        ///Parameter for the maximum number of milliseconds an algorithm is allowed to run
        #[structopt(short = "t")]
        maximum_number_of_milliseconds: Option<u32>,
    },
}

impl LTGAParameters {
    pub fn new(
        write_generational_statistics: i16,
        write_generational_solutions: i16,
        print_verbose_overview: i16,
        print_FOSs_contents: i16,
        use_ilse: i16,
        use_random_linkage: i16,
        use_value_to_reach: i16,
        use_fitness_variance_tolerance: i16,
        use_premature_stopping: i16,
        maximum_number_of_evaluations: u32,
        maximum_number_of_milliseconds: u64,
        m: u32,
        k: u32,
        o: u32,
        b: u32,
        global_optima_score: f64,
        fitness_variance_tolerance: f64,
    ) -> Self {
        LTGAParameters {
            write_generational_statistics,
            write_generational_solutions,
            print_verbose_overview,
            print_FOSs_contents,
            use_ilse,
            use_random_linkage,
            use_value_to_reach,
            use_fitness_variance_tolerance,
            use_premature_stopping,
            maximum_number_of_evaluations: maximum_number_of_evaluations as i64,
            maximum_number_of_milliseconds: maximum_number_of_milliseconds as i64,
            m,
            k,
            o,
            b,
            global_optima_score,
            fitness_variance_tolerance,
        }
    }
}

fn fill_clique_memory(cliques: &Vec<Vec<u32>>, cliques_memory: *mut *const libc::c_uint, m: u32) {
    unsafe {
        for i in 0..m as usize {
            *cliques_memory.offset(i as isize) = cliques[i].as_ptr();
        }
    }
}

fn fill_codomain_memory(
    codomain: &Vec<Vec<f64>>,
    codomain_memory: *mut *const libc::c_double,
    m: u32,
) {
    unsafe {
        for i in 0..m as usize {
            *codomain_memory.offset(i as isize) = codomain[i].as_ptr();
        }
    }
}

#include "../../src/Template/Experiments/IOHprofiler_experimenter.hpp"
#include "../EvolutionaryAlgorithm.hpp"


int budget;
int budget_scale = 5;
int budget_power = 2;

int DEFAULT_MU = 1;
int DEFAULT_LAMBDA = 1;
int DEFAULT_ML = 1;
IOHprofiler_random random_generator(1);

void _run_experiment() {
  std::string configName = "./RLSconfiguration.ini";
  IOHprofiler_experimenter<int> experimenter(configName,randomLocalSearch);
  experimenter._set_independent_runs(11);
  experimenter._run();
}

int main(){
  _run_experiment();
  return 0;
}
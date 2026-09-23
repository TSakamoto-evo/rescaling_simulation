#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <filesystem>
#include "parameter.hpp"
#include "population.hpp"

int main(){
  // haploid population size in the original scale
  int base_pop_size = 2000000;
  // rescaling factor
  int q_val = 125;

  Parameter para;
  para.pop_size = base_pop_size / q_val;

  para.r_rate = 1e-8 * q_val;
  para.u_rate = 3e-9 * q_val;
  para.length = 1e+6;

  // 0: interval between the reco-mut event is fixed
  // 1: interval between the reco-mut event is adaptively determined based on existing fitness variation
  para.simu_mode = 1;
  para.renew = para.pop_size / 100;
  para.var_target = 0.001;

  // proportion of neutral, deleterious, and beneficial mutations
  para.pn = 0.5;
  para.pb = 0.0;
  para.pd = 0.5;

  // 0: selection strength is fixed, 1: selection strength is drawn from gamma
  para.mode_b = 1;
  para.b_mean = 250.0 / 2.0 / para.pop_size;

  // 0: selection strength is fixed, 1: selection strength is drawn from gamma
  para.mode_d = 0;
  para.d_s = 100.0 / 2.0 / para.pop_size;

  std::string mut_index = "mut_index_list.txt";
  std::string fixed_list = "fixed_mutation_list.txt";

  std::ofstream log("log.txt");

  if(std::filesystem::exists("finish_simu.txt")){
    return(0);
  }

  para.check_parameter();
  Population pop(para);
  bool initialize = pop.initialize_from_files("ind.txt", "pop.txt", "rand.txt", "mut.txt", 
    mut_index, fixed_list);

  std::ofstream ofs1(mut_index, std::ios::app);
  std::ofstream ofs2(fixed_list, std::ios::app);

  if(initialize == 0){
    pop.record_full_state("ind.txt", "pop.txt", "rand.txt", "mut.txt", ofs1, ofs2);
    log << "No initial files" << std::endl;
  }else{
    log << "Use initial files" << std::endl;
  }

  std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

  // burn-in
  while(pop.return_gen() < 7 * para.pop_size){
    pop.one_step(ofs1, ofs2);

    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    auto time = now - start;
    double minute = std::chrono::duration_cast<std::chrono::minutes>(time).count();

    if(minute > 60){
      pop.record_full_state("ind.txt", "pop.txt", "rand.txt", "mut.txt", ofs1, ofs2);
      start = std::chrono::system_clock::now();
    }
  }

  pop.output_statistics(100, 100000, static_cast<int>(para.length / 250), 0, "");

  std::ofstream ofs3("finish_simu.txt");
  ofs3.close();
}
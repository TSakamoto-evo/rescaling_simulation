#ifndef POPULATION
#define POPULATION

#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <iostream>
#include <string>
#include <filesystem>
#include <locale>
#include <sstream>
#include <omp.h>
#include "parameter.hpp"
#include "individual.hpp"
#include "fenwick.hpp"

class Population{
private:
  Parameter para;
  double gen;
  int next_rec_gen;
  int rec_interval;
  int since_last_rec;

  // pop state
  std::vector<Individual> pop;
  std::vector<Individual> pop_next;

  // fenwick trees
  Fenwick birth_tree;
  Fenwick death_tree;

  double base_birth_rate;
  double base_death_rate;

  // rand
  std::mt19937 mt;

  std::discrete_distribution<> det_mut_class;
  std::uniform_real_distribution<> uni;
  std::uniform_int_distribution<> uni_pop;

  std::gamma_distribution<> b_dist;
  std::gamma_distribution<> d_dist;

  std::vector<int> indices_hap;

  // neutral variants
  std::vector<int> neu_allele_counts;
  std::vector<int> neu_origin_time;
  std::vector<size_t> unused_neu_ids;

  // selective variants
  std::vector<int> sele_allele_counts;
  std::vector<int> sele_origin_time;
  std::vector<double> sele_all_s;
  std::vector<size_t> unused_sele_ids;

  std::unordered_set<size_t> neu_fixed_list;
  std::unordered_set<size_t> sele_fixed_list;

  // used in recombination & mutation
  std::vector<std::vector<int>> label;
  std::vector<std::vector<double>> recpos;

  std::vector<std::vector<size_t>> mut_ids;
  std::vector<std::vector<double>> mut_pos;
  std::vector<std::vector<double>> mut_s;

  std::vector<int> new_label1, new_label2;
  std::vector<double> new_recpos1, new_recpos2;

  double to_double_hex(const std::string& s) const;

public:
  Population(const Parameter input_para);
  
  void selection();
  void reco_mut();

  void one_step(std::ofstream& mut_index_file, std::ofstream& fixed_file);
  double return_s();
  void fixed_site(std::ofstream& fixed_file);
  void renew_trees();
  double return_gen() const{ return(gen); };
  
  double calculate_var_fit(double& total_birth, double& total_death);

  void remove_fixed_site();
  void output_statistics(const int sample_size, const int ld_site, const int bin, 
    const bool sample_mode, const std::string tag);

  void record_full_state(const std::string ind_file, 
    const std::string pop_file, const std::string rand_file, 
    const std::string mut_file,
    std::ofstream& mut_index_file, std::ofstream& fixed_file);
  bool initialize_from_files(const std::string ind_file, 
    const std::string pop_file, const std::string rand_file, 
    const std::string mut_file, const std::string mut_index_file, 
    const std::string fixed_file);

  void shrink_vector();
};

#endif
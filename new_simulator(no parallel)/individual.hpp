#ifndef INDIVIDUAL
#define INDIVIDUAL

#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <limits>

class Individual{
private:
  std::vector<double> neu_pos;
  std::vector<size_t> neu_index;
  std::vector<double> sele_pos;
  std::vector<size_t> sele_index;

  double log_fitness;

public:
  Individual(){};
  Individual(const std::vector<double>& input_neu_pos, 
    const std::vector<size_t>& input_neu_index, 
    const std::vector<double>& input_sele_pos, 
    const std::vector<size_t>& input_sele_index, 
    const double input_log_fitness);
  
  void calculate_fitness(const std::vector<double>& sele_all_s);
  double ret_log_fitness() const{ return(log_fitness); };

  void ret_haplotype(const double start, const double end,
    std::vector<double>& ret_neu_pos, std::vector<size_t>& ret_neu_index, 
    std::vector<double>& ret_sele_pos, std::vector<size_t>& ret_sele_index, 
    const std::unordered_set<size_t>& neu_fixed_list, 
    const std::unordered_set<size_t>& sele_fixed_list, 
    double& sel_s, const std::vector<double>& sele_all_s);
  
  void ret_genome(std::vector<double>& ret_neu_pos, std::vector<size_t>& ret_neu_index, 
    std::vector<double>& ret_sele_pos, std::vector<size_t>& ret_sele_index) const;

  void remove_fixed_site(const std::unordered_set<size_t>& neu_fixed_list, 
    const std::unordered_set<size_t>& sele_fixed_list, 
    const std::vector<double>& sele_all_s);
  void count_neu_ac(std::unordered_map<size_t, int>& count, 
    std::unordered_map<size_t, double>& pos) const;
  void count_sele_ac(std::unordered_map<size_t, int>& count, 
    std::unordered_map<size_t, double>& pos) const;
  bool find_neu_locus(const size_t id) const;
  bool find_sele_locus(const size_t id) const;

  std::vector<double>& ref_neu_pos(){ return(neu_pos); };
  std::vector<size_t>& ref_neu_idx(){ return(neu_index); };
  std::vector<double>& ref_sele_pos(){ return(sele_pos); };
  std::vector<size_t>& ref_sele_idx(){ return(sele_index); };
  double& ref_log_fitness(){ return(log_fitness); };

  void shrink_to_fit();
  size_t ret_length_s() const{ return(sele_pos.size()); };

  void copy_all(std::vector<double>& ret_neu_pos, std::vector<size_t>& ret_neu_index, 
    std::vector<double>& ret_sele_pos, std::vector<size_t>& ret_sele_index, 
    double& ret_log_fitness, 
    std::vector<int>& neu_allele_counts, std::vector<int>& sele_allele_counts) const;
  void change_allele_counts(std::vector<int>& neu_allele_counts, 
    std::vector<size_t>& unused_neu_ids,
    std::vector<int>& sele_allele_counts, std::vector<size_t>& unused_sele_ids) const;
  void swap_all(std::vector<double>& ret_neu_pos, std::vector<size_t>& ret_neu_index, 
    std::vector<double>& ret_sele_pos, std::vector<size_t>& ret_sele_index, 
    double& ret_log_fitness);
};

#endif
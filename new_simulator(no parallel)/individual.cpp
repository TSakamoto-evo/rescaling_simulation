#include "individual.hpp"

Individual::Individual(const std::vector<double>& input_neu_pos, 
  const std::vector<size_t>& input_neu_index, 
  const std::vector<double>& input_sele_pos, 
  const std::vector<size_t>& input_sele_index,
  const double input_log_fitness){

  neu_pos = input_neu_pos;
  neu_index = input_neu_index;
  sele_pos = input_sele_pos;
  sele_index = input_sele_index;
  log_fitness = input_log_fitness;
}

void Individual::calculate_fitness(const std::vector<double>& sele_all_s){
  log_fitness = 0.0;
  for(const auto& id: sele_index){
    log_fitness += sele_all_s.at(id);
  }
}

void Individual::ret_haplotype(const double start, const double end,
  std::vector<double>& ret_neu_pos, std::vector<size_t>& ret_neu_index, 
  std::vector<double>& ret_sele_pos, std::vector<size_t>& ret_sele_index, 
  const std::unordered_set<size_t>& neu_fixed_list, 
  const std::unordered_set<size_t>& sele_fixed_list, 
  double& sel_s, const std::vector<double>& sele_all_s){
  
  // neutral
  {
    const auto& Hnpos = neu_pos;
    const auto& Hnid = neu_index;

    auto it_begin = std::lower_bound(Hnpos.begin(), Hnpos.end(), start);
    auto it_end = std::lower_bound(it_begin, Hnpos.end(), end);

    size_t i_begin = std::distance(Hnpos.begin(), it_begin);
    size_t i_end = std::distance(Hnpos.begin(), it_end);

    size_t anc_size = ret_neu_pos.size();
    ret_neu_pos.resize(anc_size + (i_end - i_begin));
    ret_neu_index.resize(anc_size + (i_end - i_begin));

    size_t skip = 0;
    for(size_t i = i_begin; i < i_end; i++){
      if(neu_fixed_list.count(Hnid.at(i)) != 0){
        skip++;
      }else{
        size_t j = anc_size + (i - i_begin) - skip;
        ret_neu_pos.at(j) = Hnpos.at(i);
        ret_neu_index.at(j) = Hnid.at(i);
      }
    }

    if(skip > 0){
      ret_neu_pos.resize(anc_size + (i_end - i_begin) - skip);
      ret_neu_index.resize(anc_size + (i_end - i_begin) - skip);
    }
  }

  // selection
  {
    const auto& Hspos = sele_pos;
    const auto& Hsid = sele_index;

    auto it_begin = std::lower_bound(Hspos.begin(), Hspos.end(), start);
    auto it_end = std::lower_bound(it_begin, Hspos.end(), end);

    size_t i_begin = std::distance(Hspos.begin(), it_begin);
    size_t i_end = std::distance(Hspos.begin(), it_end);

    size_t anc_size = ret_sele_pos.size();
    ret_sele_pos.resize(anc_size + (i_end - i_begin));
    ret_sele_index.resize(anc_size + (i_end - i_begin));

    size_t skip = 0;

    for(size_t i = i_begin; i < i_end; i++){
      if(sele_fixed_list.count(Hsid.at(i)) != 0){
        skip++;
      }else{
        size_t j = anc_size + (i - i_begin) - skip;
        ret_sele_pos.at(j) = Hspos.at(i);
        ret_sele_index.at(j) = Hsid.at(i);

        sel_s += sele_all_s.at(Hsid.at(i));
      }
    }

    if(skip > 0){
      ret_sele_pos.resize(anc_size + (i_end - i_begin) - skip);
      ret_sele_index.resize(anc_size + (i_end - i_begin) - skip);
    }
  }
}

void Individual::ret_genome(std::vector<double>& ret_neu_pos, std::vector<size_t>& ret_neu_index, 
  std::vector<double>& ret_sele_pos, std::vector<size_t>& ret_sele_index) const{

  ret_neu_pos = neu_pos;
  ret_neu_index = neu_index;
  ret_sele_pos = sele_pos;
  ret_sele_index = sele_index;
}

void Individual::remove_fixed_site(const std::unordered_set<size_t>& neu_fixed_list, 
  const std::unordered_set<size_t>& sele_fixed_list, const std::vector<double>& sele_all_s){

  // neutral
  {
    size_t skip = 0;
    size_t ini_size = neu_index.size();

    for(size_t i = 0; i < ini_size; i++){
      if(neu_fixed_list.count(neu_index.at(i)) != 0){
        skip++;
      }else if(skip > 0){
        neu_pos.at(i - skip) = neu_pos.at(i);
        neu_index.at(i - skip) = neu_index.at(i);
      }
    }

    if(skip > 0){
      neu_pos.resize(ini_size - skip);
      neu_index.resize(ini_size - skip);
    }
  }

  // selective
  {
    log_fitness = 0.0;

    size_t skip = 0;
    size_t ini_size = sele_index.size();

    for(size_t i = 0; i < ini_size; i++){
      if(sele_fixed_list.count(sele_index.at(i)) != 0){
        skip++;
      }else{
        log_fitness += sele_all_s.at(sele_index.at(i));
        if(skip > 0){
          sele_pos.at(i - skip) = sele_pos.at(i);
          sele_index.at(i - skip) = sele_index.at(i);
        }
      }
    }

    if(skip > 0){
      sele_pos.resize(ini_size - skip);
      sele_index.resize(ini_size - skip);
    }
  }
}

void Individual::count_neu_ac(std::unordered_map<size_t, int>& count, 
  std::unordered_map<size_t, double>& pos) const{
  
  for(size_t i = 0; i < neu_pos.size(); i++){
    if(count.count(neu_index.at(i))){
      count.at(neu_index.at(i))++;
    }else{
      count.emplace(neu_index.at(i), 1);
      pos.emplace(neu_index.at(i), neu_pos.at(i));
    }
  }
}

void Individual::count_sele_ac(std::unordered_map<size_t, int>& count, 
  std::unordered_map<size_t, double>& pos) const{
  
  for(size_t i = 0; i < sele_pos.size(); i++){
    if(count.count(sele_index.at(i))){
      count.at(sele_index.at(i))++;
    }else{
      count.emplace(sele_index.at(i), 1);
      pos.emplace(sele_index.at(i), sele_pos.at(i));
    }
  }
}

bool Individual::find_neu_locus(const size_t id) const{
  if(std::find(neu_index.begin(), neu_index.end(), id) != neu_index.end()){
    return(1);
  }else{
    return(0);
  }
}

bool Individual::find_sele_locus(const size_t id) const{
  if(std::find(sele_index.begin(), sele_index.end(), id) != sele_index.end()){
    return(1);
  }else{
    return(0);
  }
}

void Individual::shrink_to_fit(){
  neu_pos.shrink_to_fit();
  neu_index.shrink_to_fit();
  sele_pos.shrink_to_fit();
  sele_index.shrink_to_fit();
}

void Individual::copy_all(std::vector<double>& ret_neu_pos, std::vector<size_t>& ret_neu_index, 
  std::vector<double>& ret_sele_pos, std::vector<size_t>& ret_sele_index, 
  double& ret_log_fitness, 
  std::vector<int>& neu_allele_counts, std::vector<int>& sele_allele_counts) const{

  // neutral
  {
    ret_neu_pos = neu_pos;
    ret_neu_index = neu_index;

    for(size_t i = 0; i < neu_index.size(); i++) {
      neu_allele_counts.at(neu_index.at(i))++;
    }
  }

  // selective
  {
    ret_sele_pos = sele_pos;
    ret_sele_index = sele_index;

    for(size_t i = 0; i < sele_index.size(); i++) {
      sele_allele_counts.at(sele_index.at(i))++;
    }
    ret_log_fitness = log_fitness;
  }
}

void Individual::change_allele_counts(std::vector<int>& neu_allele_counts, 
  std::vector<size_t>& unused_neu_ids,
  std::vector<int>& sele_allele_counts, std::vector<size_t>& unused_sele_ids) const{
  
  for(size_t i = 0; i < neu_index.size(); i++){
    size_t id = neu_index.at(i);
    neu_allele_counts.at(id)--;

    if(neu_allele_counts.at(id) == 0){
      unused_neu_ids.push_back(id);
    }
  }

  for(size_t i = 0; i < sele_index.size(); i++){
    size_t id = sele_index.at(i);
    sele_allele_counts.at(id)--;

    if(sele_allele_counts.at(id) == 0){
      unused_sele_ids.push_back(id);
    }
  }
}

void Individual::swap_all(std::vector<double>& ret_neu_pos, std::vector<size_t>& ret_neu_index, 
  std::vector<double>& ret_sele_pos, std::vector<size_t>& ret_sele_index, double& ret_log_fitness){

  std::swap(neu_pos, ret_neu_pos);
  std::swap(neu_index, ret_neu_index);
  std::swap(sele_pos, ret_sele_pos);
  std::swap(sele_index, ret_sele_index);
  std::swap(log_fitness, ret_log_fitness);
}
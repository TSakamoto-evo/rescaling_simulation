#include "population.hpp"

Population::Population(const Parameter input_para){
  para = input_para;
  gen = 0.0;
  next_rec_gen = 1;
  since_last_rec = 0;

  if(para.simu_mode == 0){
    rec_interval = para.renew;
  }else{
    rec_interval = para.pop_size / 2;
  }

  pop.clear();
  pop_next.clear();

  neu_allele_counts.clear();
  neu_origin_time.clear();
  unused_neu_ids.clear();

  sele_allele_counts.clear();
  sele_origin_time.clear();
  sele_all_s.clear();
  unused_sele_ids.clear();

  neu_fixed_list.clear();
  sele_fixed_list.clear();

  birth_tree.initialize(para.pop_size);
  death_tree.initialize(para.pop_size);
  base_birth_rate = 0.0;
  base_death_rate = 0.0;

  std::vector<double> tmp_pos;
  std::vector<size_t> tmp_index;

  for(int i = 0; i < para.pop_size; i++){
    pop.emplace_back(tmp_pos, tmp_index, tmp_pos, tmp_index, 0.0);

    birth_tree.add(i, 1.0);
    death_tree.add(i, 1.0);
  }
  pop_next = pop;

  std::random_device seed;
  mt.seed(seed());

  std::vector<double> prop = {para.pn, para.pb, para.pd};
  det_mut_class = std::discrete_distribution<>(prop.begin(), prop.end());
  uni = std::uniform_real_distribution<>(0.0, 1.0);
  uni_pop = std::uniform_int_distribution<>(0, para.pop_size - 1);

  if(para.mode_b == 1){
    b_dist = std::gamma_distribution<>(para.b_k, para.b_mean / para.b_k);
  }

  if(para.mode_d == 1){
    d_dist = std::gamma_distribution<>(para.d_k, para.d_mean / para.d_k);
  }

  indices_hap = std::vector<int>(para.pop_size);
  std::iota(indices_hap.begin(), indices_hap.end(), 0);

  label.resize(para.pop_size);
  recpos.resize(para.pop_size);

  mut_ids.resize(para.pop_size);
  mut_pos.resize(para.pop_size);
  mut_s.resize(para.pop_size);
}

void Population::selection(){
  double rand_birth = birth_tree.total() * uni(mt);
  int birth_ind = birth_tree.lower_bound(rand_birth);

  double rand_death = death_tree.total() * uni(mt);
  int death_ind = death_tree.lower_bound(rand_death);

  if(birth_ind != death_ind){
    auto& indiv = pop_next.at(death_ind);

    std::vector<double>& neu_pos = indiv.ref_neu_pos();
    std::vector<size_t>& neu_index = indiv.ref_neu_idx();

    std::vector<double>& sele_pos = indiv.ref_sele_pos();
    std::vector<size_t>& sele_index = indiv.ref_sele_idx();

    double& log_fitness = indiv.ref_log_fitness();
    double old_log_fitness = pop.at(death_ind).ret_log_fitness();

    pop.at(birth_ind).copy_all(neu_pos, neu_index, sele_pos, sele_index, log_fitness, 
      neu_allele_counts, sele_allele_counts);
    
    double add_birth_tree = std::exp(log_fitness - base_birth_rate) - std::exp(old_log_fitness - base_birth_rate);
    birth_tree.add(death_ind, add_birth_tree);

    double add_death_tree = std::exp(-(log_fitness - base_death_rate)) - std::exp(-(old_log_fitness - base_death_rate));
    death_tree.add(death_ind, add_death_tree);

    pop.at(death_ind).change_allele_counts(neu_allele_counts, unused_neu_ids, sele_allele_counts, unused_sele_ids);
    std::swap(pop.at(death_ind), pop_next.at(death_ind));
  }

  since_last_rec++;
}

void Population::reco_mut(){
  for(int i = 0; i < para.pop_size; i++){
    label.at(i).clear();
    recpos.at(i).clear();
    label.at(i).push_back(i);
    recpos.at(i).push_back(para.length);

    mut_ids.at(i).clear();
    mut_pos.at(i).clear();
    mut_s.at(i).clear();
  }

  std::poisson_distribution<> det_reco(para.pop_size * para.r_rate * para.length * since_last_rec / para.pop_size);
  int reco_event = det_reco(mt);

  for(int i = 0; i < reco_event; i++){
    int hap1 = uni_pop(mt);
    int hap2 = uni_pop(mt);

    while(hap1 == hap2){
      hap1 = uni_pop(mt);
      hap2 = uni_pop(mt);
    }

    {
      double break_point = uni(mt) * para.length;

      new_label1.clear();
      new_label2.clear();
      new_recpos1.clear();
      new_recpos2.clear();

      int index1 = 0;
      while(recpos.at(hap1).at(index1) < break_point){
        new_label1.push_back(label.at(hap1).at(index1));
        new_recpos1.push_back(recpos.at(hap1).at(index1));

        index1++;
      }
      new_label1.push_back(label.at(hap1).at(index1));
      new_recpos1.push_back(break_point);

      int index2 = 0;
      while(recpos.at(hap2).at(index2) < break_point){
        new_label2.push_back(label.at(hap2).at(index2));
        new_recpos2.push_back(recpos.at(hap2).at(index2));

        index2++;
      }
      new_label2.push_back(label.at(hap2).at(index2));
      new_recpos2.push_back(break_point);

      for(int j = index1; j < static_cast<int>(label.at(hap1).size()); j++){
        new_label2.push_back(label.at(hap1).at(j));
        new_recpos2.push_back(recpos.at(hap1).at(j));
      }

      for(int j = index2; j < static_cast<int>(label.at(hap2).size()); j++){
        new_label1.push_back(label.at(hap2).at(j));
        new_recpos1.push_back(recpos.at(hap2).at(j));
      }

      if(uni(mt) < 0.5){
        label.at(hap1) = new_label1;
        label.at(hap2) = new_label2;

        recpos.at(hap1) = new_recpos1;
        recpos.at(hap2) = new_recpos2;
      }else{
        label.at(hap1) = new_label2;
        label.at(hap2) = new_label1;

        recpos.at(hap1) = new_recpos2;
        recpos.at(hap2) = new_recpos1;
      }
    }
  }

  std::poisson_distribution<> det_mut_num(para.pop_size * para.u_rate * para.length * since_last_rec * 2.0 / para.pop_size);
  int mut_num = det_mut_num(mt);

  for(int i = 0; i < mut_num; i++){
    int hap = uni_pop(mt);
    double pos = uni(mt) * para.length;
    double s = return_s();

    size_t id;

    if(s == 0.0){
      if(unused_neu_ids.size() == 0){
        id = neu_allele_counts.size();
        neu_allele_counts.push_back(1);
        neu_origin_time.push_back(next_rec_gen - 1);
      }else{
        id = unused_neu_ids.back();
        unused_neu_ids.pop_back();
        neu_allele_counts.at(id) = 1;
        neu_origin_time.at(id) = next_rec_gen - 1;
      }
    }else{
      if(unused_sele_ids.size() == 0){
        id = sele_allele_counts.size();
        sele_allele_counts.push_back(1);
        sele_origin_time.push_back(next_rec_gen - 1);
        sele_all_s.push_back(s);
      }else{
        id = unused_sele_ids.back();
        unused_sele_ids.pop_back();
        sele_allele_counts.at(id) = 1;
        sele_origin_time.at(id) = next_rec_gen - 1;
        sele_all_s.at(id) = s;
      }
    }

    mut_ids.at(hap).push_back(id);
    mut_pos.at(hap).push_back(pos);
    mut_s.at(hap).push_back(s);
  }

  #pragma omp parallel for schedule(dynamic)
  for(int i = 0; i < para.pop_size; i++){
    auto& indiv = pop_next.at(i);

    std::vector<double>& neu_pos = indiv.ref_neu_pos();
    neu_pos.clear();
    std::vector<size_t>& neu_index = indiv.ref_neu_idx();
    neu_index.clear();
    std::vector<double>& sele_pos = indiv.ref_sele_pos();
    sele_pos.clear();
    std::vector<size_t>& sele_index = indiv.ref_sele_idx();
    sele_index.clear();
    double& log_fitness = indiv.ref_log_fitness();
    log_fitness = 0.0;

    {
      double start = 0.0;

      if(mut_ids.at(i).size() == 0 && recpos.at(i).size() == 1 && 
        neu_fixed_list.size() == 0 && sele_fixed_list.size() == 0){

        pop.at(i).swap_all(neu_pos, neu_index, sele_pos, sele_index, log_fitness);
      }else{
        size_t n = mut_pos.at(i).size();
        std::vector<size_t> ord(n);
        std::iota(ord.begin(), ord.end(), 0);

        std::sort(ord.begin(), ord.end(), [&](size_t i1, size_t i2) {
          return mut_pos.at(i).at(i1) < mut_pos.at(i).at(i2);
        });

        std::vector<size_t> new_ids(n);
        std::vector<double> new_pos(n);
        std::vector<double> new_s(n);

        for (size_t j = 0; j < n; j++) {
          new_ids.at(j) = mut_ids.at(i).at(ord.at(j));
          new_pos.at(j) = mut_pos.at(i).at(ord.at(j));
          new_s.at(j)   = mut_s.at(i).at(ord.at(j));
        }

        size_t mut_pos_idx = 0;
        for(size_t j = 0; j < recpos.at(i).size(); j++){
          int parent = label.at(i).at(j);
          
          while(mut_pos_idx < new_ids.size() && new_pos.at(mut_pos_idx) < recpos.at(i).at(j)){
            pop.at(parent).ret_haplotype(start, new_pos.at(mut_pos_idx), neu_pos, neu_index, sele_pos, sele_index, 
              neu_fixed_list, sele_fixed_list, log_fitness, sele_all_s);
            
            if(new_s.at(mut_pos_idx) == 0.0){
              neu_pos.push_back(new_pos.at(mut_pos_idx));
              neu_index.push_back(new_ids.at(mut_pos_idx));
            }else{
              sele_pos.push_back(new_pos.at(mut_pos_idx));
              sele_index.push_back(new_ids.at(mut_pos_idx));
              log_fitness += new_s.at(mut_pos_idx);
            }

            start = new_pos.at(mut_pos_idx);
            mut_pos_idx++;
          }

          pop.at(parent).ret_haplotype(start, recpos.at(i).at(j), neu_pos, neu_index, sele_pos, sele_index, 
            neu_fixed_list, sele_fixed_list, log_fitness, sele_all_s);
          start = recpos.at(i).at(j);
        }     
      }
    }
  }

  for(int i = 0; i < para.pop_size; i++){
    if(mut_ids.at(i).size() == 0 && recpos.at(i).size() == 1 && 
        neu_fixed_list.size() == 0 && sele_fixed_list.size() == 0){
      continue;
    }else{
      double old_log_fitness = pop.at(i).ret_log_fitness();
      double log_fitness = pop_next.at(i).ret_log_fitness();

      double add_birth = std::exp(log_fitness - base_birth_rate) - std::exp(old_log_fitness - base_birth_rate);
      double add_death = std::exp(-(log_fitness - base_death_rate)) - std::exp(-(old_log_fitness - base_death_rate));

      birth_tree.add(i, add_birth);
      death_tree.add(i, add_death);
    }
  }

  std::swap(pop, pop_next);

  for(const auto& i: neu_fixed_list){
    unused_neu_ids.push_back(i);
    neu_allele_counts.at(i) = 0;
  }
  neu_fixed_list.clear();

  for(const auto& i: sele_fixed_list){
    unused_sele_ids.push_back(i);
    sele_allele_counts.at(i) = 0;
  }
  sele_fixed_list.clear();

  since_last_rec = 0;
}

void Population::fixed_site(std::ofstream& fixed_file){
  std::vector<double>& neu_pos = pop.at(0).ref_neu_pos();
  std::vector<size_t>& neu_index = pop.at(0).ref_neu_idx();
  std::vector<double>& sele_pos = pop.at(0).ref_sele_pos();
  std::vector<size_t>& sele_index = pop.at(0).ref_sele_idx();

  neu_fixed_list.clear();

  for(size_t id = 0; id < neu_allele_counts.size(); id++){
    if(neu_allele_counts.at(id) == para.pop_size){
      bool found = 0;
      for(size_t j = 0; j < neu_pos.size(); j++){
        if(neu_index.at(j) == id){
          double position = neu_pos.at(j);
          int origin = neu_origin_time.at(id);

          fixed_file << next_rec_gen << "\t" << origin << "\t" << position << "\t" << 0.0 << "\n";
          found = 1;
          break;
        }
      }
      if(!found){
        std::cerr << "index does not exist" << std::endl;
        std::exit(1);
      }

      neu_fixed_list.emplace(id);
    }
  }

  sele_fixed_list.clear();

  for(size_t id = 0; id < sele_allele_counts.size(); id++){
    if(sele_allele_counts.at(id) == para.pop_size){
      bool found = 0; 
      for(size_t j = 0; j < sele_pos.size(); j++){
        if(sele_index.at(j) == id){
          double position = sele_pos.at(j);
          int origin = sele_origin_time.at(id);
          double s = sele_all_s.at(id);

          fixed_file << next_rec_gen << "\t" << origin << "\t" << position << "\t" << s << "\n";
          found = 1;
          break;
        }
      }
      if(!found){
        std::cerr << "index does not exist" << std::endl;
        std::exit(1);
      }

      sele_fixed_list.emplace(id);
    }
  }
}

void Population::renew_trees(){
  std::vector<double> log_fitness(para.pop_size);
  base_birth_rate = pop.at(0).ret_log_fitness();
  base_death_rate = pop.at(0).ret_log_fitness();

  for(int i = 0; i < para.pop_size; i++){
    log_fitness.at(i) = pop.at(i).ret_log_fitness();

    if(base_birth_rate < log_fitness.at(i)){
      base_birth_rate = log_fitness.at(i);
    }
    if(base_death_rate > log_fitness.at(i)){
      base_death_rate = log_fitness.at(i);
    }
  }

  birth_tree.initialize(para.pop_size);
  death_tree.initialize(para.pop_size);

  for(int i = 0; i < para.pop_size; i++){
    birth_tree.add(i, std::exp(log_fitness.at(i) - base_birth_rate));
    death_tree.add(i, std::exp(-(log_fitness.at(i) - base_death_rate)));
  }
}

void Population::one_step(std::ofstream& mut_index_file, std::ofstream& fixed_file){
  selection();

  if(since_last_rec >= rec_interval || since_last_rec >= para.pop_size / 2){
    reco_mut();
  }

  gen += 2.0 / para.pop_size;

  if(gen >= next_rec_gen){
    double total_birth, total_death;

    fixed_site(fixed_file);
    double var_fit = calculate_var_fit(total_birth, total_death);

    if(para.simu_mode == 1){
      if(var_fit > 0.0){
        double tau = para.var_target / var_fit * 2.0 * para.pop_size;

        if(tau < 1.0){
          rec_interval = 1;
        }else if(tau > para.pop_size / 2){
          rec_interval = para.pop_size / 2;
        }else{
          rec_interval = static_cast<int>(std::floor(tau));
        }
      }else{
        rec_interval = para.pop_size / 2;
      }
    }

    double birth_tree_total = birth_tree.total();
    double death_tree_total = death_tree.total();

    if(std::abs(total_birth - birth_tree_total) > 1e-6 || std::abs(total_death - death_tree_total) > 1e-6 ||
      birth_tree_total < para.pop_size / 100.0 || birth_tree_total > para.pop_size * 100.0 || 
      death_tree_total < para.pop_size / 100.0 || death_tree_total > para.pop_size * 100.0 || 
      next_rec_gen % para.pop_size == 0){

      renew_trees();
    }

    mut_index_file << next_rec_gen << "\t" << var_fit << "\t" << rec_interval << std::endl;
    next_rec_gen++;
  }
}

double Population::calculate_var_fit(double& total_birth, double& total_death){
  std::vector<double> log_fitness(para.pop_size);

  double birth_mean = 0.0;
  double death_mean = 0.0;

  for(int i = 0; i < para.pop_size; i++){
    log_fitness.at(i) = pop.at(i).ret_log_fitness();

    double birth_fit = std::exp(log_fitness.at(i) - base_birth_rate);
    double death_fit = std::exp(-(log_fitness.at(i) - base_death_rate));

    birth_mean += birth_fit;
    death_mean += death_fit;
  }

  total_birth = birth_mean;
  total_death = death_mean;

  birth_mean /= para.pop_size;
  death_mean /= para.pop_size;

  double sum = 0.0;
  double sum_sq = 0.0;

  for(int i = 0; i < para.pop_size; i++){
    double birth_fit = std::exp(log_fitness.at(i) - base_birth_rate) / birth_mean;
    double death_fit = std::exp(-(log_fitness.at(i) - base_death_rate)) / death_mean;

    double total_fit = birth_fit - death_fit;

    sum += total_fit;
    sum_sq += total_fit * total_fit;
  }

  sum /= para.pop_size;
  sum_sq /= para.pop_size;

  return(sum_sq - sum * sum);
}


double Population::return_s(){
  int mut_class = det_mut_class(mt);

  double s_coef = 0.0;
  if(mut_class == 1){
    if(para.mode_b == 0){
      s_coef = para.b_s;
    }else{
      s_coef = b_dist(mt);
    }
  }else if(mut_class == 2){
    if(para.mode_d == 0){
      s_coef = -para.d_s;
    }else{
      s_coef = -d_dist(mt);
    }
  }

  return(s_coef);
}

void Population::remove_fixed_site(){
  for(int i = 0; i < para.pop_size; i++){
    pop.at(i).remove_fixed_site(neu_fixed_list, sele_fixed_list, sele_all_s);
  }

  for(const auto& i: neu_fixed_list){
    unused_neu_ids.push_back(i);
    neu_allele_counts.at(i) = 0;
  }
  neu_fixed_list.clear();

  for(const auto& i: sele_fixed_list){
    unused_sele_ids.push_back(i);
    sele_allele_counts.at(i) = 0;
  }
  sele_fixed_list.clear();
}

void Population::output_statistics(const int sample_size, const int ld_site, const int bin, 
  const bool sample_mode, const std::string tag){

  int pop_size = para.pop_size;
  int s_size = sample_size;

  if(sample_mode == 0){
    if(s_size <= 0 || s_size > pop_size){
      std::cerr << "All individuals are used to calculate statistics" << std::endl;
      s_size = pop_size;
    }
  }

  std::vector<int> sampled_index(s_size);

  if(sample_mode == 0){
    std::shuffle(indices_hap.begin(), indices_hap.end(), mt);
    for(int i = 0; i < s_size; i++){
      sampled_index.at(i) = indices_hap.at(i);
    }
  }else{
    for(int i = 0; i < s_size; i++){
      sampled_index.at(i) = uni_pop(mt);
    }
  }

  // SFS
  std::vector<int> ac_n(s_size - 1);
  std::vector<int> ac_b(s_size - 1);
  std::vector<int> ac_d(s_size - 1);

  std::unordered_map<size_t, int> neu_sample_ac;
  std::unordered_map<size_t, double> neu_sample_pos;

  // neutral
  {
    for(int i = 0; i < s_size; i++){
      int ind = sampled_index.at(i);
      pop.at(ind).count_neu_ac(neu_sample_ac, neu_sample_pos);
    }
    std::vector<size_t> fixed;

    for(const auto& i: neu_sample_ac){
      int ac = i.second;

      if(ac < s_size){
        ac_n.at(ac - 1)++;
      }else if(ac == s_size){
        fixed.push_back(i.first);
      }else{
        std::cerr << "error" << std::endl;
        std::exit(1);
      }
    }

    for(const auto& i: fixed){
      neu_sample_ac.erase(i);
      neu_sample_pos.erase(i);
    }
  }

  std::unordered_map<size_t, int> sele_sample_ac;
  std::unordered_map<size_t, double> sele_sample_pos;

  // selective
  {
    for(int i = 0; i < s_size; i++){
      int ind = sampled_index.at(i);
      pop.at(ind).count_sele_ac(sele_sample_ac, sele_sample_pos);
    }
    std::vector<size_t> fixed;

    for(const auto& i: sele_sample_ac){
      int ac = i.second;
      double s = sele_all_s.at(i.first);

      if(ac < s_size){
        if(s > 0.0){
          ac_b.at(ac - 1)++;
        }else{
          ac_d.at(ac - 1)++;
        }
      }else if(ac == s_size){
        fixed.push_back(i.first);
      }else{
        std::cerr << "error" << std::endl;
        std::exit(1);
      }
    }
    
    for(const auto& i: fixed){
      sele_sample_ac.erase(i);
      sele_sample_pos.erase(i);
    }
  }
  

  std::ofstream output_ac("sample_sfs" + tag + ".txt");
  output_ac << "ac\tneutral\tbeneficial\tdeleterious\n";
  for(int i = 0; i < s_size - 1; i++){
    output_ac << i + 1 << "\t" << ac_n.at(i) << "\t" << ac_b.at(i) << "\t" << ac_d.at(i) << "\n";
  }
  output_ac.close();

  // LD
  size_t n_site = ld_site;
  if(static_cast<int>(neu_sample_ac.size() + sele_sample_ac.size()) < ld_site || ld_site < 0){
    n_site = neu_sample_ac.size() + sele_sample_ac.size();
  }

  std::vector<long long int> site_index(neu_sample_ac.size() + sele_sample_ac.size());
  int tmp = 0;
  for(const auto& i: neu_sample_ac){
    site_index.at(tmp) = i.first;
    tmp++;
  }
  for(const auto& i: sele_sample_ac){
    site_index.at(tmp) = -static_cast<long long int>(i.first) - 1;
    tmp++;
  }
  std::shuffle(site_index.begin(), site_index.end(), mt);

  std::vector<std::vector<bool>> genetic_matrix(s_size, std::vector<bool>(n_site));
  std::vector<double> site_pos;

  for(size_t i = 0; i < n_site; i++){
    if(site_index.at(i) >= 0){
      // neutral site
      size_t id = site_index.at(i);
      site_pos.push_back(neu_sample_pos.at(id));

      for(int j = 0; j < s_size; j++){
        int ind = sampled_index.at(j);

        bool exist = pop.at(ind).find_neu_locus(id);
        genetic_matrix.at(j).at(i) = exist;
      }
    }else{
      size_t id = -(site_index.at(i) + 1);
      site_pos.push_back(sele_sample_pos.at(id));

      for(int j = 0; j < s_size; j++){
        int ind = sampled_index.at(j);

        bool exist = pop.at(ind).find_sele_locus(id);
        genetic_matrix.at(j).at(i) = exist;
      }
    }
  }

  std::vector<double> r_sq_bin(bin);
  std::vector<int> bin_count(bin);

  for(size_t i = 0; i < n_site; i++){
    for(size_t j = i + 1; j < n_site; j++){
      int count_11 = 0;
      int count_10 = 0;
      int count_01 = 0;

      double dist = std::abs(site_pos.at(i) - site_pos.at(j));

      for(int k = 0; k < s_size; k++){
        if(genetic_matrix.at(k).at(i) == 1 && genetic_matrix.at(k).at(j) == 1){
          count_11++;
        }else if(genetic_matrix.at(k).at(i) == 1 && genetic_matrix.at(k).at(j) == 0){
          count_10++;
        }else if(genetic_matrix.at(k).at(i) == 0 && genetic_matrix.at(k).at(j) == 1){
          count_01++;
        }
      }

      double fa = 1.0 * (count_11 + count_10) / s_size;
      double fb = 1.0 * (count_11 + count_01) / s_size;
      double fab = 1.0 * count_11 / s_size;

      double d = fab - fa * fb;
      double r_sq = d * d / fa / (1.0 - fa) / fb / (1.0 - fb);

      int bin_index = std::floor(bin * dist / para.length);
      if(bin_index >= bin){
        bin_index = bin - 1;
      }

      r_sq_bin.at(bin_index) += r_sq;
      bin_count.at(bin_index)++;
    }
  }

  std::ofstream output_ld("output_ld" + tag + ".txt");
  output_ld << "bin\tcount\tr_sq\n";
  for(int i = 0; i < bin; i++){
    if(bin_count.at(i) > 0){
      output_ld << i << "\t" << bin_count.at(i) << "\t" << r_sq_bin.at(i) / bin_count.at(i) << "\n";
    }
  }
  output_ld.close();
}

void Population::record_full_state(const std::string ind_file, 
  const std::string pop_file, const std::string rand_file, 
  const std::string mut_file,
  std::ofstream& mut_index_file, std::ofstream& fixed_file){

  mut_index_file << std::flush;
  fixed_file << std::flush;
  
  // individual class
  std::vector<double> neu_pos, sele_pos;
  std::vector<size_t> neu_index, sele_index;

  int pop_size = para.pop_size;
  std::ofstream ofs1("tmp_" + ind_file);
  ofs1 << std::hexfloat;

  std::unordered_map<size_t, double> link_neu_id_pos;
  std::unordered_map<size_t, double> link_sele_id_pos; 

  for(int i = 0; i < pop_size; i++){
    pop.at(i).ret_genome(neu_pos, neu_index, sele_pos, sele_index);
    // neutral
    ofs1 << i << "\t" << 0;
    for(size_t j = 0; j < neu_index.size(); j++){
      ofs1 << "\t" << neu_index.at(j);

      if(link_neu_id_pos.count(neu_index.at(j)) == 0){
        link_neu_id_pos.emplace(neu_index.at(j), neu_pos.at(j));
      }
    }
    ofs1 << "\n";

    // selective
    ofs1 << i << "\t" << 1;
    for(size_t j = 0; j < sele_index.size(); j++){
      ofs1 << "\t" << sele_index.at(j);

      if(link_sele_id_pos.count(sele_index.at(j)) == 0){
        link_sele_id_pos.emplace(sele_index.at(j), sele_pos.at(j));
      }
    }
    ofs1 << "\n";
    ofs1 << i << "\t" << 2 << "\t" << pop.at(i).ret_log_fitness() << "\n";
  }
  ofs1.close();

  std::ofstream ofs4("tmp_" + mut_file);
  ofs4 << std::hexfloat;

  ofs4 << 0;
  for(size_t i = 0; i < unused_neu_ids.size(); i++){
    ofs4 << "\t" << unused_neu_ids.at(i);
  }
  ofs4 << "\n";

  ofs4 << 1;
  for(size_t i = 0; i < unused_sele_ids.size(); i++){
    ofs4 << "\t" << unused_sele_ids.at(i);
  }
  ofs4 << "\n";

  for(size_t i = 0; i < neu_allele_counts.size(); i++){
    if(link_neu_id_pos.count(i) > 0){
      ofs4 << 2 << "\t" << i << "\t" << neu_allele_counts.at(i) << "\t" << neu_origin_time.at(i) << 
        "\t" << link_neu_id_pos.at(i) << "\n";
    }else{
      ofs4 << 2 << "\t" << i << "\t" << neu_allele_counts.at(i) << "\t" << neu_origin_time.at(i) << 
        "\t" << -1.0 << "\n";
    }
  }

  for(size_t i = 0; i < sele_allele_counts.size(); i++){
    if(link_sele_id_pos.count(i) > 0){
      ofs4 << 3 << "\t" << i << "\t" << sele_allele_counts.at(i) << "\t" << sele_origin_time.at(i) << 
        "\t" << link_sele_id_pos.at(i) << "\t" << sele_all_s.at(i) << "\n";
    }else{
      ofs4 << 3 << "\t" << i << "\t" << sele_allele_counts.at(i) << "\t" << sele_origin_time.at(i) << 
        "\t" << -1.0 << "\t" << sele_all_s.at(i) << "\n";
    }
  }
  ofs4.close();

  // population class
  std::ofstream ofs2("tmp_" + pop_file);
  ofs2 << std::hexfloat;
  ofs2 << 0 << "\t" << gen << "\t" << next_rec_gen << "\t" << base_birth_rate << "\t" << base_death_rate << "\t" << 
    rec_interval << "\t" << since_last_rec << "\n";

  ofs2 << 1;
  for(const auto& i: neu_fixed_list){
    ofs2 << "\t" << i;
  }
  ofs2 << "\n";

  ofs2 << 2;
  for(const auto& i: sele_fixed_list){
    ofs2 << "\t" << i;
  }
  ofs2 << "\n";

  ofs2 << 3;
  for(const auto& i: indices_hap){
    ofs2 << "\t" << i;
  }
  ofs2 << "\n";

  int n;
  std::vector<double> bit;

  birth_tree.return_all(n, bit);
  ofs2 << 4 << "\t" << 0 << "\t" << n << "\n";
  ofs2 << 4 << "\t" << 1;
  for(const auto& i: bit){
    ofs2 << "\t" << i;
  }
  ofs2 << "\n";

  death_tree.return_all(n, bit);
  ofs2 << 5 << "\t" << 0 << "\t" << n << "\n";
  ofs2 << 5 << "\t" << 1;
  for(const auto& i: bit){
    ofs2 << "\t" << i;
  }
  ofs2 << "\n";

  ofs2.close();

  // rand
  std::ofstream ofs3("tmp_" + rand_file);
  ofs3 << mt;
  ofs3.close();

  std::filesystem::rename("tmp_" + ind_file, ind_file);
  std::filesystem::rename("tmp_" + pop_file, pop_file);
  std::filesystem::rename("tmp_" + rand_file, rand_file);
  std::filesystem::rename("tmp_" + mut_file, mut_file);
}

double Population::to_double_hex(const std::string& s) const{
  return std::stod(s);
}

bool Population::initialize_from_files(const std::string ind_file, 
  const std::string pop_file, const std::string rand_file, 
  const std::string mut_file, const std::string mut_index_file, 
  const std::string fixed_file){

  if(std::filesystem::exists(ind_file) && std::filesystem::exists(pop_file) && 
    std::filesystem::exists(rand_file) && std::filesystem::exists(mut_file) &&
    std::filesystem::exists(mut_index_file) && std::filesystem::exists(fixed_file)){

    std::unordered_map<size_t, double> ret_neu_pos;
    std::unordered_map<size_t, double> ret_sele_pos;

    // mutation file
    {
      neu_allele_counts.clear();
      neu_origin_time.clear();
      unused_neu_ids.clear();

      sele_allele_counts.clear();
      sele_origin_time.clear();
      sele_all_s.clear();
      unused_sele_ids.clear();

      std::ifstream ifs(mut_file);
      if(!ifs){
        std::cerr << "Fail to open the mut_file!" << std::endl;
        std::exit(1);
      }

      std::string line;

      while (getline(ifs, line)){
        std::istringstream iss(line);
        std::string tmp_list;
        std::vector<std::string> list;

        while(getline(iss, tmp_list, '\t')){
          list.push_back(tmp_list);
        }

        if(list.at(0) == "0"){
          for(size_t i = 1; i < list.size(); i++){
            unused_neu_ids.push_back(std::stoull(list.at(i)));
          }
        }else if(list.at(0) == "1"){
          for(size_t i = 1; i < list.size(); i++){
            unused_sele_ids.push_back(std::stoull(list.at(i)));
          }
        }else if(list.at(0) == "2"){
          size_t id = std::stoull(list.at(1));
          int ac = std::stoi(list.at(2));
          int origin = std::stoi(list.at(3));
          double pos = to_double_hex(list.at(4));

          neu_allele_counts.push_back(ac);
          neu_origin_time.push_back(origin);
          ret_neu_pos.emplace(id, pos);
        }else if(list.at(0) == "3"){
          size_t id = std::stoull(list.at(1));
          int ac = std::stoi(list.at(2));
          int origin = std::stoi(list.at(3));
          double pos = to_double_hex(list.at(4));
          double s = to_double_hex(list.at(5));

          sele_allele_counts.push_back(ac);
          sele_origin_time.push_back(origin);
          ret_sele_pos.emplace(id, pos);
          sele_all_s.push_back(s);
        }else{
          std::cerr << "error in the mut_file format" << std::endl;
          std::exit(1);
        }
      }
    }

    // individual file
    {
      pop.clear();

      std::vector<double> neu_pos, sele_pos;
      std::vector<size_t> neu_index, sele_index;

      std::ifstream ifs(ind_file);
      if(!ifs){
        std::cerr << "Fail to open the ind_file!" << std::endl;
        std::exit(1);
      }

      std::string line;

      while (getline(ifs, line)){
        std::istringstream iss(line);
        std::string tmp_list;
        std::vector<std::string> list;

        while(getline(iss, tmp_list, '\t')){
          list.push_back(tmp_list);
        }

        if(list.at(1) == "0"){
          neu_index.clear();
          neu_pos.clear();

          for(int i = 2; i < static_cast<int>(list.size()); i++){
            size_t add_index = std::stoull(list.at(i));

            if(ret_neu_pos.count(add_index) == 0 || ret_neu_pos.at(add_index) < 0.0){
              std::cerr << "no index: neu " << add_index << std::endl;
            }

            neu_index.push_back(add_index);
            neu_pos.push_back(ret_neu_pos.at(add_index));
          }
        }else if(list.at(1) == "1"){
          sele_index.clear();
          sele_pos.clear();

          for(int i = 2; i < static_cast<int>(list.size()); i++){
            size_t add_index = std::stoull(list.at(i));

            if(ret_sele_pos.count(add_index) == 0 || ret_sele_pos.at(add_index) < 0.0){
              std::cerr << "no index: sele " << add_index << std::endl;
            }

            sele_index.push_back(add_index);
            sele_pos.push_back(ret_sele_pos.at(add_index));
          }
        }else if(list.at(1) == "2"){
          double log_fitness = to_double_hex(list.at(2));
          pop.emplace_back(neu_pos, neu_index, sele_pos, sele_index, log_fitness);
        }else{
          std::cerr << "error in the ind_file format" << std::endl;
          std::exit(1);
        }
      }

      ifs.close();
    }

    // population file
    {
      neu_fixed_list.clear();
      sele_fixed_list.clear();
      indices_hap.clear();

      int birth_n = -1;
      int death_n = -1;
      std::vector<double> birth_bit, death_bit;

      std::ifstream ifs(pop_file);
      if(!ifs){
        std::cerr << "Fail to open the pop_file!" << std::endl;
        std::exit(1);
      }

      std::string line;

      while (getline(ifs, line)){
        std::istringstream iss(line);
        std::string tmp_list;
        std::vector<std::string> list;

        while(getline(iss, tmp_list, '\t')){
          list.push_back(tmp_list);
        }

        if(list.at(0) == "0"){
          gen = to_double_hex(list.at(1));
          next_rec_gen = std::stoi(list.at(2));
          base_birth_rate = to_double_hex(list.at(3));
          base_death_rate = to_double_hex(list.at(4));
          rec_interval = std::stoi(list.at(5));
          since_last_rec = std::stoi(list.at(6));
        }else if(list.at(0) == "1"){
          for(int i = 1; i < static_cast<int>(list.size()); i++){
            neu_fixed_list.emplace(std::stoull(list.at(i)));
          }
        }else if(list.at(0) == "2"){
          for(int i = 1; i < static_cast<int>(list.size()); i++){
            sele_fixed_list.emplace(std::stoull(list.at(i)));
          }
        }else if(list.at(0) == "3"){
          for(int i = 1; i < static_cast<int>(list.size()); i++){
            indices_hap.push_back(std::stoi(list.at(i)));
          }
        }else if(list.at(0) == "4"){
          if(list.at(1) == "0"){
            birth_n = std::stoi(list.at(2));
          }else if(list.at(1) == "1"){
            for(int i = 2; i < static_cast<int>(list.size()); i++){
              double add = to_double_hex(list.at(i));
              birth_bit.push_back(add);
            }
          }else{
            std::cerr << "error in pop_file format" << std::endl;
            std::exit(1);
          }
        }else if(list.at(0) == "5"){
          if(list.at(1) == "0"){
            death_n = std::stoi(list.at(2));
          }else if(list.at(1) == "1"){
            for(int i = 2; i < static_cast<int>(list.size()); i++){
              double add = to_double_hex(list.at(i));
              death_bit.push_back(add);
            }
          }else{
            std::cerr << "error in pop_file format" << std::endl;
            std::exit(1);
          }         
        }else{
          std::cerr << "error in pop_file format" << std::endl;
          std::exit(1);
        }
      }

      ifs.close();

      if(birth_n < 0 || death_n < 0){
        std::cerr << "error in pop_file format" << std::endl;
        std::exit(1);
      }

      birth_tree.input_all(birth_n, birth_bit);
      death_tree.input_all(death_n, death_bit);
    }

    // rand
    {
      std::ifstream ifs(rand_file);
      if(!ifs){
        std::cerr << "Fail to open the rand_file!" << std::endl;
        std::exit(1);
      }

      ifs >> mt;
      ifs.close();
    }

    // renew mut_index file
    {
      std::ofstream ofs("tmp_" + mut_index_file);
      std::ifstream ifs(mut_index_file);
      if(!ifs){
        std::cerr << "Fail to open the mut_index_file!" << std::endl;
        std::exit(1);
      }

      std::string line;

      while (getline(ifs, line)){
        std::istringstream iss(line);
        std::string tmp_list;
        std::vector<std::string> list;

        while(getline(iss, tmp_list, '\t')){
          list.push_back(tmp_list);
        }

        if(std::stod(list.at(0)) < next_rec_gen){
          ofs << line << "\n";
        }
      }

      ifs.close();
      ofs.close();

      std::filesystem::rename("tmp_" + mut_index_file, mut_index_file);
    }

    // renew fixed file
    {
      std::ofstream ofs("tmp_" + fixed_file);
      std::ifstream ifs(fixed_file);
      if(!ifs){
        std::cerr << "Fail to open the mut_index_file!" << std::endl;
        std::exit(1);
      }

      std::string line;

      while (getline(ifs, line)){
        std::istringstream iss(line);
        std::string tmp_list;
        std::vector<std::string> list;

        while(getline(iss, tmp_list, '\t')){
          list.push_back(tmp_list);
        }

        if(std::stod(list.at(0)) < next_rec_gen){
          ofs << line << "\n";
        }
      }

      ifs.close();
      ofs.close();

      std::filesystem::rename("tmp_" + fixed_file, fixed_file);
    }

    pop_next.resize(pop.size());

    return(1);
  }else{
    return(0);
  }
}

void Population::shrink_vector(){
  for(auto& i: pop){
    i.shrink_to_fit();
  }
  for(auto& j: pop_next){
    j.shrink_to_fit();
  }
}
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>


int main(){
  int site_number = 0;
  int bin_size = 250;
  int max_snps = 100000;

  std::vector<std::vector<bool>> genotype;
  std::vector<int> pos;

  std::vector<int> ac_n;
  std::vector<int> ac_b;
  std::vector<int> ac_d;

  std::random_device seed;
  std::mt19937 mt(seed());
  std::uniform_real_distribution<> uni(0.0, 1.0);

  std::ifstream ifs("sample.vcf");
  if(!ifs){
    std::cerr << "Error: fail to open the vcf file!" << std::endl;
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

    if(list.at(0)[0] != '#'){
      std::istringstream iss2(list.at(7));
      std::vector<std::string> info;
      while(getline(iss2, tmp_list, ';')){
        info.push_back(tmp_list);
      }

      if(ac_n.size() == 0){
        ac_n = std::vector<int>(2 * (list.size() - 9) + 1);
        ac_b = std::vector<int>(2 * (list.size() - 9) + 1);
        ac_d = std::vector<int>(2 * (list.size() - 9) + 1);
      }

      int mut_type = -1;
      for(const auto& i: info){
        if(i == "MT=1"){
          mut_type = 1;
        }else if(i == "MT=2"){
          mut_type = 2;
        }else if(i == "MT=3"){
          mut_type = 3;
        }
      }

      if(mut_type == -1){
        std::cerr << "undefined mut type" << std::endl;
      }


      std::vector<bool> this_site(2 * (list.size() - 9));

      int ac = 0;
      for(size_t i = 9; i < list.size(); i++){
        if(list.at(i) == "0|0"){
          this_site.at((i - 9) * 2) = 0;
          this_site.at((i - 9) * 2 + 1) = 0;
        }else if(list.at(i) == "0|1"){
          this_site.at((i - 9) * 2) = 0;
          this_site.at((i - 9) * 2 + 1) = 1;
          ac += 1;
        }else if(list.at(i) == "1|0"){
          this_site.at((i - 9) * 2) = 1;
          this_site.at((i - 9) * 2 + 1) = 0;
          ac += 1;
        }else if(list.at(i) == "1|1"){
          this_site.at((i - 9) * 2) = 1;
          this_site.at((i - 9) * 2 + 1) = 1;
          ac += 2;
        }else{
          std::cerr << "unexpected genotype" << std::endl;
        }
      }

      if(ac > 0 && ac < 2 * (static_cast<int>(list.size()) - 9)){
        site_number++;

        if(mut_type == 1){
          ac_n.at(ac)++;
        }else if(mut_type == 2){
          ac_b.at(ac)++;
        }else if(mut_type == 3){
          ac_d.at(ac)++;
        }

        if(static_cast<int>(genotype.size()) < max_snps){
          genotype.push_back(this_site);
          pos.push_back(std::stoi(list.at(1)));
        }else{
          int rep = (site_number * uni(mt));
          if(rep == site_number){
            rep--;
          }

          if(rep < max_snps){
            genotype.at(rep) = this_site;
            pos.at(rep) = std::stoi(list.at(1));
          }
        }
      }
    }
  }

  std::ofstream output_ac("sample_sfs.txt");
  output_ac << "ac\tneutral\tbeneficial\tdeleterious\n";
  for(int i = 1; i < static_cast<int>(ac_n.size()) - 1; i++){
    output_ac << i << "\t" << ac_n.at(i) << "\t" << ac_b.at(i) << "\t" << ac_d.at(i) << "\n";
  }
  output_ac.close();

  std::vector<double> r_sq_bin;
  std::vector<int> bin_count;

  for(size_t i = 0; i < genotype.size(); i++){
    for(size_t j = i + 1; j < genotype.size(); j++){
      int count_11 = 0;
      int count_10 = 0;
      int count_01 = 0;

      double dist = std::abs(pos.at(i) - pos.at(j));
      int s_size = static_cast<int>(genotype.at(0).size());

      for(int k = 0; k < s_size; k++){
        if(genotype.at(i).at(k) == 1 && genotype.at(j).at(k) == 1){
          count_11++;
        }else if(genotype.at(i).at(k) == 1 && genotype.at(j).at(k) == 0){
          count_10++;
        }else if(genotype.at(i).at(k) == 0 && genotype.at(j).at(k) == 1){
          count_01++;
        }
      }

      double fa = 1.0 * (count_11 + count_10) / s_size;
      double fb = 1.0 * (count_11 + count_01) / s_size;
      double fab = 1.0 * count_11 / s_size;

      double d = fab - fa * fb;
      double r_sq = d * d / fa / (1.0 - fa) / fb / (1.0 - fb);

      int bin_index = std::floor(dist / bin_size);
      while(bin_index >= static_cast<int>(r_sq_bin.size())){
        r_sq_bin.push_back(0.0);
        bin_count.push_back(0);
      }

      r_sq_bin.at(bin_index) += r_sq;
      bin_count.at(bin_index)++;
    }
  }

  std::ofstream output_ld("output_ld.txt");
  output_ld << "bin\tcount\tr_sq\n";
  for(size_t i = 0; i < r_sq_bin.size(); i++){
    if(bin_count.at(i) > 0){
      output_ld << i << "\t" << bin_count.at(i) << "\t" << r_sq_bin.at(i) / bin_count.at(i) << "\n";
    }
  }
  output_ld.close();

  std::ofstream ofs("finish_simu.txt");
  ofs.close();
}
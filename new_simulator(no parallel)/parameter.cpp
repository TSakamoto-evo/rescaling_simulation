#include "parameter.hpp"

void Parameter::check_parameter(){
    if(r_rate < 0.0){
      std::cerr << "Error: Negative recombination rate" << std::endl;
      std::exit(1);
    }

    if(u_rate < 0.0){
      std::cerr << "Error: Negative mutation rate" << std::endl;
      std::exit(1);
    }

    if(length <= 0.0){
      std::cerr << "Error: Non-positive genome length" << std::endl;
      std::exit(1);
    }

    double sum = pn + pb + pd;
    if(sum <= 0.0 || pn < 0.0 || pb < 0.0 || pd < 0.0){
      std::cerr << "Error: Invalid mutational proportion" << std::endl;
      std::exit(1);
    }

    pn /= sum;
    pb /= sum;
    pd /= sum;
    
    if((mode_b == 0 && b_s < 0.0) || (mode_d == 0 && d_s < 0.0)){
      std::cerr << "Error: Negative coefficient" << std::endl;
      std::exit(1);
    }

    if((mode_b == 1 && (b_k <= 0.0 || b_mean <= 0.0)) || (mode_d == 1 && (d_k <= 0.0 || d_mean <= 0.0))){
      std::cerr << "Error: Non-positive parameters for gamma distribution" << std::endl;
      std::exit(1);
    }
  }

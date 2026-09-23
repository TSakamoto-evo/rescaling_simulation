#ifndef PARAMETER
#define PARAMETER

#include <vector>
#include <iostream>

class Parameter{
public:
  int pop_size;

  double r_rate = 0.0;
  double u_rate = 0.0;
  double length = 1.0;

  // 0: fixed, 1: flex
  bool simu_mode = 0;
  int renew = 1;
  double var_target = 0.01;

  double pn;
  double pb;
  double pd;

  // 0: fixed 1: gamma
  bool mode_b = 0;
  double b_s = 0.0;
  double b_k = 1.0;
  double b_mean = 0.01;

  // 0: fixed 1: gamma
  bool mode_d = 0;
  double d_s = 0.0;
  double d_k = 1.0;
  double d_mean = 0.01;

  void check_parameter();
};

#endif
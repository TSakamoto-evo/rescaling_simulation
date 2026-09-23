#ifndef FENWICK
#define FENWICK

#include <vector>
#include <algorithm>

class Fenwick{
private: 
  int n;
  std::vector<double> bit;

public:
  void initialize(int input_n){
    n = input_n;
    
    if (static_cast<int>(bit.size()) < n + 1) {
      bit.resize(input_n + 1);
    }
    std::fill(bit.begin(), bit.end(), 0.0);
  };

  void return_all(int& ret_n, std::vector<double>& ret_bit){
    ret_n = n;
    ret_bit = bit;
  }

  void input_all(const int input_n, const std::vector<double>& input_bit){
    n = input_n;
    bit = input_bit;
  }

  void add(int idx, double delta){
    for(int i = idx + 1; i <= n; i += i & -i) {
      bit[i] += delta;
    }
  };

  double sum(int i) const{
    double s = 0;
    for (int x = i + 1; x > 0; x -= x & -x){
      s += bit[x];
    }
    return s;
  }

  double total() const{
    return sum(n - 1);
  }

  int lower_bound(double x) const{
    int idx = 0;
    double acc = 0.0;

    int bit_mask = 1;
    while ((bit_mask << 1) <= n){
      bit_mask <<= 1;
    }

    for (int k = bit_mask; k > 0; k >>= 1) {
      int next = idx + k;
      if (next <= n && acc + bit[next] < x) {
        acc += bit[next];
        idx = next;
      }
    }
    return idx;
  }
};

#endif
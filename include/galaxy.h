#ifndef GALAXY_H
#define GALAXY_H

#include <cmath>
#include <vector>

#include "units.h"

class Galaxy {
public:
  Galaxy(const size_t& size);
  
  virtual ~Galaxy();
  
  void init(const double& D0, const double& vA, const double& sigma, const double& H);
  
  size_t get_size() const {
    return z.size();
  }

  double get_z(const size_t& i) const {
    return z.at(i);
  }

  double get_v(const size_t& i) const {
    return v.at(i);
  }

  double get_v() const {
    size_t i = get_plane();
    return get_v(i);
  }
  
  double get_Q(const size_t& i) const {
    return Q.at(i);
  }

  double get_D(const size_t& i) const {
    return D.at(i);
  }

  double get_D() const {
    size_t i = get_plane();
    return get_D(i);
  }
  
  double get_dz() const {
    return dz;
  }

  size_t get_plane() const {
    return i_plane;
  }
  
    
private:
  size_t i_plane;
  double dz;
  std::vector<double> z; 
  std::vector<double> Q;
  std::vector<double> D;
  std::vector<double> v;
};

#endif

#include <algorithm>
#include <cmath>

#include <chrono>
#include <thread>         // std::this_thread::sleep_for

#include <iomanip>
#include <iostream>
#include <vector>
#include "tridiag.h"

// TIME UNITS
static const double sec = 1.;
static const double year = 3.154e+7 * sec;
static const double kyr = 1e3 * year;
static const double Myr = 1e6 * year;
static const double Gyr = 1e9 * year;

// LENGTH UNITS
static const double meter = 1.;
static const double cm = 1e-2 * meter;
static const double km = 1e3 * meter;
static const double parsec = 3.1e18 * cm;
static const double kpc = 1e3 * parsec;

// COMBINED UNITS
static const double cm2 = cm * cm;

#define pow2(A) ((A)*(A))

class ProgressBar {
public:
  ProgressBar() {
  }
  virtual ~ProgressBar() {
    std::cout << std::endl;
  }
  void print(float progress) const {
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cout << "=";
      else if (i == pos) std::cout << ">";
      else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
  }
private:
  int barWidth = 70;
};

class Galaxy {
public:
  Galaxy(const size_t& size) : z(std::vector<double> (size)), Q(std::vector<double> (size)), D(std::vector<double> (size)), v(std::vector<double> (size)) {
  }
  
  void init(const double& D0, const double& vA, const double& sigma, const double& H) {
    dz = 2 * H / (double)(z.size() - 1);
    for (int i = 0; i < z.size(); ++i) {
      const double z_ = -H + (double)i * dz;
      z.at(i) = z_;
      Q.at(i) = std::exp(-(z_ * z_) / 2. / (sigma * sigma));
      D.at(i) = D0;
      v.at(i) = (z_ < 0.) ? -vA : vA;
    }
  }
  
  size_t get_size() const {
    return z.size();
  }

  double get_Q(const size_t& i) const {
    return Q.at(i);
  }

  double get_D(const size_t& i) const {
    return D.at(i);
  }

  double get_dz() const {
    return dz;
  }
  
  virtual ~Galaxy() {
    z.clear(); Q.clear(); D.clear(); v.clear();
  }
  
private:
  double dz;
  std::vector<double> z; 
  std::vector<double> Q;
  std::vector<double> D;
  std::vector<double> v;
};

class DiffusionAdvectionSolver1D {
public:
  DiffusionAdvectionSolver1D(const double& dt, Galaxy *galaxy) {
    this->dt = dt;
    this->galaxy = galaxy; 
    size = galaxy->get_size();
  }
  
  virtual ~DiffusionAdvectionSolver1D() {
  }
  
  virtual void do_step(std::vector<double>& N) = 0;
  
  void run(const double& max_time, const ProgressBar* bar) {
    std::vector<double> N = std::vector<double> (size);
    for(double t = double(); t < max_time; t += dt) {
      //std::cout << std::scientific << t / kyr << "\n";
      do_step(N);
      bar->print(t / max_time);
    }
  }
  
protected:
  size_t size;
  double dt;
  Galaxy *galaxy;
};

class CrankNicholson : public DiffusionAdvectionSolver1D {
public:
  CrankNicholson(const double& dt, Galaxy *galaxy) : DiffusionAdvectionSolver1D(dt, galaxy) {
    rhs.resize(size - 2);
    N_next.resize(size - 2);
    central_diag.resize(size - 2);
    upper_diag.resize(size - 3);
    lower_diag.resize(size - 3);
  }
  
  virtual void do_step(std::vector<double>& N) {
    const double dt_half = 0.5 * dt;
    for (int i = 1; i < size - 1; ++i) {
      double Dz = galaxy->get_D(i);
      double Dz_up = (i < size - 1) ? galaxy->get_D(i + 1) : Dz;
      double Dz_do = (i > 0) ? galaxy->get_D(i - 1) : Dz;

      double Dz_dz2 = Dz / pow2(galaxy->get_dz());
      double dDz_4_dz2 = (Dz_up - Dz_do) / 4 / pow2(galaxy->get_dz());

      double U = Dz_dz2 + dDz_4_dz2;
      double C = 2. * Dz_dz2;
      double L = Dz_dz2 - dDz_4_dz2;
      
      central_diag.at(i - 1) = 1 - dt_half * C;
      if (i != 1) {
	lower_diag.at(i - 2) = -dt_half * L;
      }
      if (i != size - 2) {
      	upper_diag.at(i - 1) = -dt_half * U;
      }
      rhs.at(i - 1) = N.at(i) * (1 + dt_half * C);
      rhs.at(i - 1) += (i != 0) ? dt_half * L * N.at(i - 1) : 0;
      rhs.at(i - 1) += (i != size - 2) ? dt_half * U * N.at(i + 1) : 0;
      rhs.at(i - 1) += dt * galaxy->get_Q(i);
    }
    
    gsl_linalg_solve_tridiag(central_diag, upper_diag, lower_diag, rhs, N_next);
    
    for (int i = 1; i < size - 1; ++i)
      N.at(i) = std::max(N_next.at(i - 1), 0.);
  }
  
protected:
  std::vector<double> rhs;
  std::vector<double> N_next;
  std::vector<double> central_diag;
  std::vector<double> upper_diag;
  std::vector<double> lower_diag;
};

int main( int argc , char **argv ) {
  const double dt = 0.1 * kyr;
  const double D0 = 3e28 * cm2 / sec;
  const double vA = 10 * km / sec;
  const double sigma = 100 * parsec;
  const double H = 4 * kpc;
  const int size = pow(2, 4) + 1;
  
  Galaxy galaxy(size);
  galaxy.init(D0, vA, sigma, H);

  CrankNicholson CN(dt, &galaxy);

  ProgressBar bar;
  CN.run(100 * Myr, &bar);
}

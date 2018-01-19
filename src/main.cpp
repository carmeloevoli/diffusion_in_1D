#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

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

template <typename T, int size>
class Grid {
  T memblock [size];
public:
  Grid(const T& value) {
    for (int i = 0; i < size; i++)
      memblock[i] = value;
  }
  void set (int x, T value){
    memblock[x] = value;
  }
  T get (int x) {
    return memblock[x];
  }
  void dump(const T& dt) {
    std::cout << dt << "\n";
    for (int i = 0; i < size; i++)
      std::cout << memblock[i] << "\n";
  }
};

template<typename T, int size>
class Galaxy {
public:
  Galaxy(const T& D0, const T& vA, const T& sigma, const T& H) {
    init_grids(D0, vA, sigma, H);
  }
  virtual ~Galaxy() {
  }
private:
  T z [size];
  T Q [size];
  T D [size];
  T v [size];
  
  void init_grids(const T& D0, const T& vA, const T& sigma, const T& H) {
    const T dz = 2 * H / (T)(size - 1);
    for (size_t i = 0; i < size; ++i) {
      const T z_ = -H + (T)i * dz;
      z[i] = z_;
      Q[i] = std::exp(-(z_ * z_) / 2. / (sigma * sigma));
      D[i] = D0;
      v[i] = (z_ < 0.) ? -vA : vA;
    }
  }
};

template<typename T, int size>
class Stepper {
public:
  virtual ~Stepper() {}
  virtual void do_step(Grid<T, size>& N, const T& dt) = 0;
};

template<typename T, int size>
class CrankNicholson : public Stepper<T,size> {
public:
  CrankNicholson(Galaxy<T,size> *galaxy_) : galaxy(galaxy_) {
  }
  void do_step(Grid<T, size>& N, const T& dt) {
    T dt_half = 0.5 * dt;
    for (int i = 1; i < size - 1; ++i) {
      T L, C, U;
      
      central_diagonal[i] = 1 - dt_half * C;
      if (i != 0) {
	lower_diagonal[i - 1] = -dt_half * L;
      }
      if (i != size - 2) {
	upper_diagonal[i] = -dt_half * U;
      }
      
      rhs[i] = N.get(i) * (2 - central_diagonal[i]);
      rhs[i] += (i != 0) ? dt_half * L * N.get(i - 1) : 0;
      rhs[i] += (i != size - 2) ? dt_half * U * N.get(i + 1) : 0;
      //rhs[i] += dt * Q.get(ik, iz) / (double) number_of_operators;
    }
    
    //gsl_linalg_solve_tridiag(central_diagonal, upper_diagonal, lower_diagonal, rhs, W_up);
    
    for (int i = 1; i < size - 1; ++i)
      N.set(i) = std::max(N_up[i - 1], 0);
  }
private:
  Galaxy<T, size> *galaxy;
  T rhs [size - 2];
  T N_up[size - 2];
  T central_diagonal [size - 2];
  T upper_diagonal [size - 3];
  T lower_diagonal [size - 3];
};

int main( int argc , char **argv ) {
  const double dt = 0.1 * kyr;
  const double D0 = 3e28 * cm2 / sec;
  const double vA = 10 * km / sec;
  const double sigma = 100 * parsec;
  const double H = 4 * kpc;
  
  Grid<double, 16> N(0);
  Galaxy<double, 16> galaxy(D0, vA, sigma, H);
  
  CrankNicholson<double, 16> CN (&galaxy);
  
  for(double t = 0; t < Myr; t += dt) {
    CN.do_step(N, dt);
    N.dump(dt);
  }
}

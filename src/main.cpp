#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "bar.h"
#include "galaxy.h"
#include "units.h"
#include "tridiag.h"

std::string generate_output_filename(const double& t, const size_t& size) {
	std::stringstream sstream;
	sstream << "output/" << "N";
	sstream << "_size_" << size;
	sstream << "_t_" << t;
	sstream << ".txt";
	std::string out = sstream.str();
	return out;
}

class DiffusionAdvectionSolver1D {
public:
  DiffusionAdvectionSolver1D(const double& dt, Galaxy *galaxy) {
    this->dt = dt;
    this->galaxy = galaxy; 
    size = galaxy->get_size();
    double courant_advection = galaxy->get_v() * dt / galaxy->get_dz();
    double courant_diffusion = 2. * galaxy->get_D() * dt / pow2(galaxy->get_dz());
    std::cout << " ... time step in Myr is " << dt / Myr << "\n";
    std::cout << " ... grid step in kpc is " << galaxy->get_dz() / kpc << "\n";
    std::cout << " ... Courant number for advection is " << courant_advection << "\n";
    std::cout << " ... Courant number for diffusion is " << courant_diffusion << "\n";
   }
  
  virtual ~DiffusionAdvectionSolver1D() {
  }
  
  virtual void do_step(std::vector<double>& N) = 0;

  void dump(const std::vector<double>& N, const double& t) {
    std::string filename = generate_output_filename(t / Myr, size);
    std::ofstream outfile(filename.c_str());
    outfile << "# N [au]\n";
    outfile << std::scientific << std::setprecision(4);
    for (size_t i = 0; i < N.size(); ++i) {
      outfile << galaxy->get_z(i) / kpc << " " << N.at(i) << "\n";
    }
    outfile.close();
  }
  
  void run(const double& max_time, const double& dump_time, const ProgressBar* bar) {
    std::vector<double> N = std::vector<double> (size);
    int counter = 0;
    for(double t = double(); t < max_time; t += dt, ++counter) {
      do_step(N);
      bar->print(t / max_time);
      if (counter % (int)(dump_time / dt) == 0)
	dump(N, t);
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
      double dv_dz = galaxy->get_v(i) / galaxy->get_dz();
      double Dz = galaxy->get_D(i);
      double Dz_up = (i < size - 1) ? galaxy->get_D(i + 1) : Dz;
      double Dz_do = (i > 0) ? galaxy->get_D(i - 1) : Dz;
      double Dz_dz2 = Dz / pow2(galaxy->get_dz());
      double dDz_4_dz2 = (Dz_up - Dz_do) / 4 / pow2(galaxy->get_dz());

      double U = 0., C = 0., L = 0.;

      U += Dz_dz2 + dDz_4_dz2;
      C += 2. * Dz_dz2;
      L += Dz_dz2 - dDz_4_dz2;
      
      if (i > galaxy->get_plane()) {
	L += galaxy->get_v(i - 1) / galaxy->get_dz();
	C += galaxy->get_v(i) / galaxy->get_dz();
	U += 0;
      }
      else if (i < galaxy->get_plane()) {
	L += 0;
	C += galaxy->get_v(i) / galaxy->get_dz();
	U += galaxy->get_v(i + 1) / galaxy->get_dz();
      }
      else  {
	L += -galaxy->get_v(i - 1) / 2. / galaxy->get_dz();
	C += 0;
	U += -galaxy->get_v(i + 1) / 2. / galaxy->get_dz();
      }
      
      central_diag.at(i - 1) = 1 + dt_half * C;
      if (i != 1) {
	lower_diag.at(i - 2) = -dt_half * L;
      }
      if (i != size - 2) {
      	upper_diag.at(i - 1) = -dt_half * U;
      }
      rhs.at(i - 1) = N.at(i) * (1 - dt_half * C);
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
  const double dt = 1 * kyr;
  const double D0 = 3e28 * cm2 / sec;
  const double vA = 30 * km / sec;
  const double sigma = 50 * parsec;
  const double H = 4 * kpc;
  const int size = pow(2, 8) + 1;
  
  Galaxy galaxy(size);
  galaxy.init(D0, vA, sigma, H);

  CrankNicholson CN(dt, &galaxy);

  ProgressBar bar;
  CN.run(550 * Myr, 10 * Myr, &bar);
}

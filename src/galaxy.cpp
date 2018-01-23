#include "galaxy.h"

Galaxy::Galaxy(const size_t& size) : z(std::vector<double> (size)), Q(std::vector<double> (size)), D(std::vector<double> (size)), v(std::vector<double> (size)) {}

Galaxy::~Galaxy() {
  z.clear(); Q.clear(); D.clear(); v.clear();
}

void Galaxy::init(const double& D0, const double& vA, const double& sigma, const double& H) {
  dz = 2 * H / (double)(z.size() - 1);
  i_plane = (z.size() - 1) / 2;
  for (int i = 0; i < z.size(); ++i) {
    const double z_ = -H + (double)i * dz;
    z.at(i) = z_;
    Q.at(i) = 1e-30 / std::sqrt(2.0 * M_PI * pow2(sigma)) * std::exp(-pow2(z_) / 2. / pow2(sigma));
    D.at(i) = D0;
    v.at(i) = (z_ < 0.) ? vA : vA;
  }
}

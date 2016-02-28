/*
 * velocity.cpp
 *
 *  Created on: Feb 28, 2016
 *      Author: rice
 */

#include "velocity.h"
#include "logger.h"

Velocity::Velocity(const std::vector<float> &vel, int _nx, int _nz) :
  dat(vel), nx(_nx), nz(_nz) {
  if (dat.size() != size_t(nx * nz)) {
    ERROR() << format("%d elements in velocity, but nx*nz=%d") % dat.size() % (nx * nz);
    exit(0);
  }
}

Velocity::Velocity(int _nx, int _nz) : dat(_nx *_nz), nx(_nx), nz(_nz) {
}

/*< expand domain of 'org' to 'exp': source(org)-->destination(exp) >*/
void expand(Velocity &exp, const Velocity &org, int nb) {
  int nx = org.nx;
  int nz = org.nz;
  int nxpad = nx + 2 * nb;
  int nzpad = nz + nb;
  const std::vector<float> &a = org.dat;
  std::vector<float> &b = exp.dat;

  /// internal
  for (int ix = 0; ix < nx; ix++) {
    for (int iz = 0; iz < nz; iz++) {
      b[(nb + ix) * nzpad + iz] = a[ix * nx + iz];
    }
  }

  /// boundary
  for (int ix = 0; ix < nxpad; ix++) {
    for (int iz = 0; iz < nb; iz++) {
      b[ix * nzpad + (nzpad - iz - 1)] = b[ix * nzpad + (nzpad - nb - 1)];  /* bottom*/
    }
  }

  for (int ix = 0; ix < nb; ix++) {
    for (int iz = 0; iz < nzpad; iz++) {
      b[ix * nzpad + iz] = b[nb * nzpad + iz];/* left */
      b[(nxpad - ix - 1) * nzpad + iz] = b[(nxpad - nb - 1) * nzpad + iz]; /* right */
    }
  }
}

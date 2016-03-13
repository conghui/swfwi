/*
 * mpi-fwi-params.h
 *
 *  Created on: Feb 26, 2016
 *      Author: rice
 */

#ifndef SRC_MPI_FWI2D_MPI_FWI_PARAMS_H_
#define SRC_MPI_FWI2D_MPI_FWI_PARAMS_H_

extern "C"
{
#include <rsf.h>
}

class MpiFwiParams {
public:
  static MpiFwiParams &instance();

private:
  MpiFwiParams();
  MpiFwiParams(const MpiFwiParams &);
  void operator=(const MpiFwiParams &);
  ~MpiFwiParams();
  void getInputParams();
  void putOutputParams();
  void check();
  void calVars();

private:
  static MpiFwiParams *ins;

public:
  sf_file vinit;        /* initial velocity model, unit=m/s */
  sf_file shots;        /* recorded shots from exact velocity model */
  sf_file vupdates;     /* updated velocity in iterations */
  sf_file grads;        /* gradient in iterations */
  sf_file illums;       /* source illumination in iterations */
  sf_file objs;         /* values of objective function in iterations */
  bool verb;            // verbosity
  int nz;
  int nx;
  float dz;
  float dx;
  bool precon;
  int niter;
  int rbell;
  int nt;
  int ng;
  int ns;
  float dt;
  float amp;
  float fm;
  int sxbeg;
  int szbeg;
  int gxbeg;
  int gzbeg;
  int jsx;
  int jsz;
  int jgx;
  int jgz;
  int csd;
  int nb;

public: // calculated
  int nk; // # of shots for each process
  int numProc; // total number of MPI process
  int rank;    // current process
  const char *obsDataFileName;
};

#endif /* SRC_MPI_FWI2D_MPI_FWI_PARAMS_H_ */

/*
 * mpi-global-params.h
 *
 *  Created on: Feb 24, 2016
 *      Author: rice
 */

#ifndef SRC_MPIFWI_MPI_GLOBAL_PARAMS_H_
#define SRC_MPIFWI_MPI_GLOBAL_PARAMS_H_

extern "C"
{
#include <rsf.h>
}

class MpiGlobalParams {
public:
  static MpiGlobalParams &instance();

private:
  MpiGlobalParams();
  MpiGlobalParams(const MpiGlobalParams &);
  void operator=(const MpiGlobalParams &);
  ~MpiGlobalParams();
  void getInputParams();
  void putOutputParams();
  void check();
  void calVars();

private:
  static MpiGlobalParams *ins;

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

public: // calculated
  int nk; // # of shots for each process
  int numProc; // total number of MPI process
  int rank;    // current process
};
#endif /* SRC_MPIFWI_MPI_GLOBAL_PARAMS_H_ */

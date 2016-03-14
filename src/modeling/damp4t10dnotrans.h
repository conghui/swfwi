/*
 * damp4t10dnotrans.h
 *
 *  Created on: Mar 14, 2016
 *      Author: rice
 */

#ifndef SRC_MODELING_DAMP4T10DNOTRANS_H_
#define SRC_MODELING_DAMP4T10DNOTRANS_H_


extern "C" {
#include <rsf.h>
}
#include <boost/function.hpp>
#include "velocity.h"
#include "shot-position.h"

class Damp4t10dNotrans {
public:
  Damp4t10dNotrans(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float dt, float dx, float fm, int nb, int nt);

  Velocity expandDomain(const Velocity &vel);


  void stepForward(float *p0, float *p1) const;
  void stepBackward(float *p0, float *p1) const;
  void bindVelocity(const Velocity &_vel);
  void addSource(float *p, const float *source, const ShotPosition &pos) const;
  void addSource(float *p, const float *source, int is) const;
  void addEncodedSource(float *p, const float *encsrc) const;
  void recordSeis(float *seis_it, const float *p) const;
  void maskGradient(float *grad) const;
  void refillBoundary(float *vel) const;
  void sfWriteVel(sf_file file) const;

  void removeDirectArrival(float* data) const;
  void subEncodedSource(float *p, const float *source) const;
  void refillVelStencilBndry();

  std::vector<float> initBndryVector(int nt) const;
  void writeBndry(float* _bndr, const float* p, int it) const;
  void readBndry(const float* _bndr, float* p, int it) const;

public:
  const Velocity &getVelocity() const;
  Velocity &getVelocity();
  const ShotPosition &getAllSrcPos() const;
  const ShotPosition &getAllGeoPos() const;
  int getns() const;
  int getng() const;
  float getdt() const;
  float getdx() const;
  int getnt() const;
  int getnx() const;
  int getnz() const;

private:
  void manipSource(float *p, const float *source, const ShotPosition &pos, boost::function2<float, float, float> op) const;
  void recordSeis(float *seis_it, const float *p, const ShotPosition &geoPos) const;
  void subSource(float *p, const float *source, const ShotPosition &pos) const;
  void removeDirectArrival(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float* data, int nt, float t_width) const;

private:
  const static int EXFDBNDRYLEN = 0;
  const static int FDLEN = 5;

private:
  const Velocity *vel;
  const ShotPosition *allSrcPos;
  const ShotPosition *allGeoPos;
  float dt;
  float dx;
  float fm;
  int nb;
  int nt;
  mutable int bndrSize;

private:
  std::vector<float> bndr;

};
#endif /* SRC_MODELING_DAMP4T10DNOTRANS_H_ */

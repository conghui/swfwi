/*
 * damp4t10d.h
 *
 *  Created on: Feb 29, 2016
 *      Author: rice
 */

#ifndef SRC_FM2D_DAMP4T10D_H_
#define SRC_FM2D_DAMP4T10D_H_

extern "C" {
#include <rsf.h>
}
#include <boost/function.hpp>
#include "velocity.h"
#include "shot-position.h"

class Damp4t10d {
public:
  Damp4t10d(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float dt, float dx, float fm, int nb, int nt);

  Velocity expandDomain(const Velocity &vel);


  void stepForward(float *p0, float *p1) const;
  void stepBackward(float *p0, float *p1) const;
  void bindVelocity(const Velocity &_vel);
  void addSource(float *p, const float *source, const ShotPosition &pos) const;
  void addSource(float *p, const float *source, int is) const;
  void subSource(float *p, const float *source, const ShotPosition &pos) const;
  void addEncodedSource(float *p, const float *encsrc) const;
  void recordSeis(float *seis_it, const float *p) const;
  void maskGradient(float *grad) const;
  void scaleGradient(float *grad) const;
  void refillBoundary(float *vel) const;
  void sfWriteVel(const std::vector<float> &exvel, sf_file file) const;

  void fwiRemoveDirectArrival(float* data, int shot_id) const;
  void removeDirectArrival(float* data) const;
  void subEncodedSource(float *p, const float *source) const;
  void refillVelStencilBndry();

  std::vector<float> initBndryVector(int nt) const;
  void writeBndry(float* _bndr, const float* p, int it) const;
  void readBndry(const float* _bndr, float* p, int it) const;

  void FwiForwardModeling(const std::vector<float> &encsrc, std::vector<float> &dcal, int shot_id) const;
  void EssForwardModeling(const std::vector<float> &encsrc, std::vector<float> &dcal) const;

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
  void removeDirectArrival(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float* data, int nt, float t_width) const;

private:
  const static int EXFDBNDRYLEN = 6;

private:
  const Velocity *vel;
  const ShotPosition *allSrcPos;
  const ShotPosition *allGeoPos;
  float dt;
  float dx;
  float fm;
  int bx0, bxn;
  int bz0, bzn;
  int nt;
  mutable int bndrSize;
  mutable int bndrWidth;

private:
  std::vector<float> bndr;

};

#endif /* SRC_FM2D_DAMP4T10D_H_ */

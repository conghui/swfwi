/*
 * zjh4t10dspongevelnotrans.h
 *
 *  Created on: Mar 10, 2016
 *      Author: rice
 */

#ifndef SRC_MODELING_ZJH4T10DSPONGEVELNOTRANS_H_
#define SRC_MODELING_ZJH4T10DSPONGEVELNOTRANS_H_

extern "C" {
#include <rsf.h>
}
#include <boost/function.hpp>
#include "velocity.h"
#include "shot-position.h"

class Zjh4t10dSpongeVelNoTrans {
public:
  Zjh4t10dSpongeVelNoTrans(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float dt, float dx, float fm, int nb, int nt);

  Velocity expandDomain(const Velocity &vel);


  void stepForward(float *p0, float *p1) const;
  void stepBackward(float *p0, float *p1) const;
  void stepBackwardLapIllum(float *lap, float *illum, float *p0, float *p1) const;
  void bindVelocity(const Velocity &_vel);
  void addSource(float *p, const float *source, const ShotPosition &pos) const;
  void addSource(float *p, const float *source, int is) const;
  void addEncodedSource(float *p, const float *encsrc) const;
  void recordSeis(float *seis_it, const float *p) const;
  const Velocity &getVelocity() const;
  void maskGradient(float *grad) const;
  void refillBoundary(float *vel) const;
  void sfWrite(const std::vector<float> &model, sf_file file) const;

  void removeDirectArrival(float* data) const;
  void subEncodedSource(float *p, const float *source) const;

  std::vector<float> initBndryVector(int nt) const;
  void writeBndry(float *bndr, const float *p, int it) const;
  void readBndry(const float *bndr, float *p, int it) const;
public:
  int getns() const;
  int getng() const;
  float getdt() const;
  float getdx() const;
  int getnt() const;

private:
  void manipSource(float *p, const float *source, const ShotPosition &pos, boost::function2<float, float, float> op) const;
  void recordSeis(float *seis_it, const float *p, const ShotPosition &geoPos) const;
  void subSource(float *p, const float *source, const ShotPosition &pos) const;
  void removeDirectArrival(const ShotPosition &allSrcPos, const ShotPosition &allGeoPos, float* data, int nt, float t_width) const;

private:
  const static int FDLEN = 5;
  const static int EXFDBNDRYLEN = 0;

private:
  const Velocity *vel;
  const ShotPosition *allSrcPos;
  const ShotPosition *allGeoPos;
  float dt;
  float dx;
  float fm;
  int nb;
  int nt;

private:
  std::vector<float> bndr;
  mutable int bndrSize;

};
#endif /* SRC_MODELING_ZJH4T10DSPONGEVELNOTRANS_H_ */



/*
 * enquistabc2d.h
 *
 *  Created on: Mar 1, 2016
 *      Author: rice
 */

#ifndef SRC_COMMON_ENQUISTABC2D_H_
#define SRC_COMMON_ENQUISTABC2D_H_

#include "velocity.h"
#include "shot-position.h"
#include <boost/function.hpp>

class EnquistAbc2d {
public:
  EnquistAbc2d(float dt, float dx, float dz);
  void stepForward(const float *p0, const float *p1, float *p2) const;
  void stepBackward(float *illum, float *lap, const float *p0, const float *p1, float *p2) const;
  void addSource(float *p, const float *source, const ShotPosition &pos) const;
  void subSource(float *p, const float *source, const ShotPosition &pos) const;
  void bindVelocity(const Velocity &vel);
  void recordSeis(float *seis_it, const float *p, const ShotPosition &geoPos) const;
  void writeBndry(float *bndr, const float *p, int it) const;
  void readBndry(const float *bndr, float *p, int it) const;
  std::vector<float> initBndryVector(int nt) const;

private:
  void manipSource(float *p, const float *source, const ShotPosition &pos, boost::function2<float, float, float> op) const;
private:
  const Velocity *vel;
  float dt;
  float dx;
  float dz;
  float dtx;
  float dtz;
};

#endif /* SRC_COMMON_ENQUISTABC2D_H_ */

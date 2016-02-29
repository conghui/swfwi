/*
 * sf-velocity-reader.h
 *
 *  Created on: Feb 26, 2016
 *      Author: rice
 */

#ifndef SRC_SF_VELOCITY_READER_H_
#define SRC_SF_VELOCITY_READER_H_

extern "C" {
#include <rsf.h>
}

#include "velocity.h"

class SfVelocityReader {
public:
  static Velocity read(sf_file file, int nx, int nz);

public:
  SfVelocityReader(sf_file &f);
  void readAndBcast(float *vv, size_t count, int rank);
  void read(float *vv, size_t count);
  virtual ~SfVelocityReader();

private:
  sf_file &file;
};


#endif /* SRC_SF_VELOCITY_READER_H_ */

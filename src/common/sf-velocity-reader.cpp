/*
 * sf-velocity-reader.cpp
 *
 *  Created on: Feb 26, 2016
 *      Author: rice
 */

#include <mpi.h>
#include "sf-velocity-reader.h"

#include "logger.h"

SfVelocityReader::SfVelocityReader(sf_file& f) : file(f) {
}

void SfVelocityReader::readAndBcast(float* vv, size_t count,
    int rank)
{
  if (rank == 0) {
    INFO() << format("rank %d is reading velocity") % rank;
    sf_floatread(vv, count, file);
  }

  // broadcast the velocity
  if (rank == 0) {
    INFO() << format("rank %d is broadcasting velocity") % rank;
  }

  MPI_Bcast(vv, count, MPI_FLOAT, 0, MPI_COMM_WORLD);
}

void SfVelocityReader::read(float* vv, size_t count) {
  INFO() << "reading velocity";
  sf_floatread(vv, count, file);
}

Velocity SfVelocityReader::read(sf_file file, int nx, int nz) {
  Velocity v(nx, nz);
  sf_seek(file, 0, SEEK_SET);
  sf_floatread(&v.dat[0], nx * nz, file);

  return v;
}

SfVelocityReader::~SfVelocityReader() {
}


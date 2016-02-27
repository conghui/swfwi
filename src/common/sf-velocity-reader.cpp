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

SfVelocityReader::~SfVelocityReader() {
}


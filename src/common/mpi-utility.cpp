/*
 * mpi-utility.cpp
 *
 *  Created on: Feb 26, 2016
 *      Author: rice
 */

#include <cstdlib>
#include "mpi-utility.h"

void MpiInplaceReduce(void *buf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, buf, count, datatype, op, root, comm);
  } else {
    MPI_Reduce(buf, NULL, count, datatype, op, root, comm);
  }
}

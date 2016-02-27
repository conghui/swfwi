/*
 * mpi-utility.h
 *
 *  Created on: Feb 26, 2016
 *      Author: rice
 */

#ifndef SRC_COMMON_MPI_UTILITY_H_
#define SRC_COMMON_MPI_UTILITY_H_

#include <mpi.h>

void MpiInplaceReduce(void *buf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);

#endif /* SRC_COMMON_MPI_UTILITY_H_ */

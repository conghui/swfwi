/*
 * shotdata-reader.cpp
 *
 *  Created on: Feb 26, 2016
 *      Author: rice
 */

#include <cstdio>
#include <mpi.h>
#include "shotdata-reader.h"
#include "logger.h"
#include "common.h"

void ShotDataReader::read(const char *datapath, float* dobs, int nshots, int nt, int ng) {
  int nproc;
  int rank;

  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_File fh;
  int file_open_error = MPI_File_open(MPI_COMM_WORLD, datapath,
                      MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

  if (file_open_error != MPI_SUCCESS) {

    char error_string[BUFSIZ];
    int length_of_error_string, error_class;

    MPI_Error_class(file_open_error, &error_class);
    MPI_Error_string(error_class, error_string, &length_of_error_string);
    ERROR() << format("%d: %s") % rank % error_string;

    MPI_Error_string(file_open_error, error_string, &length_of_error_string);
    ERROR() << format("%d: %s") % rank % error_string;

    MPI_Abort(MPI_COMM_WORLD, file_open_error);
  }

  MPI_File_seek(fh, rank * nt * ng * sizeof(float), MPI_SEEK_SET);
  for (int is = rank, k = 0; is < nshots; is += nproc, k++) {
    std::vector<float> trans(nt * ng);
    MPI_Status status;

    MPI_File_read(fh, &trans[0], trans.size(), MPI_FLOAT, &status);
    matrix_transpose(&trans[0], &dobs[k * trans.size()], nt, ng);

    MPI_File_seek(fh, nt * ng * (nproc - 1) * sizeof(float), MPI_SEEK_CUR);
  }

  MPI_File_close(&fh);
}



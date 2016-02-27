/*
 * shotdata-reader.h
 *
 *  Created on: Feb 26, 2016
 *      Author: rice
 */

#ifndef SRC_COMMON_SHOTDATA_READER_H_
#define SRC_COMMON_SHOTDATA_READER_H_

#include <mpi.h>

class ShotDataReader {
public:
  static void read(const char *datapath, float *dobs, int nshots, int nt, int ng);
};

#endif /* SRC_COMMON_SHOTDATA_READER_H_ */

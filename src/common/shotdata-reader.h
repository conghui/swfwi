/*
 * shotdata-reader.h
 *
 *  Created on: Feb 26, 2016
 *      Author: rice
 */

#ifndef SRC_COMMON_SHOTDATA_READER_H_
#define SRC_COMMON_SHOTDATA_READER_H_

extern "C" {
#include <rsf.h>
}
#include <mpi.h>
#include <vector>

class ShotDataReader {
public:
  static void parallelRead(const char *datapath, float *dobs, int nshots, int nt, int ng);
  static void serialRead(sf_file file, float *dobs, int nshots, int nt, int ng);
  static void readAndEncode(sf_file file, const std::vector<int> &codes, float *dobs, int nshots, int nt, int ng);
};

#endif /* SRC_COMMON_SHOTDATA_READER_H_ */

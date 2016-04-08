/*
 * Environment.cpp
 *
 *  Created on: Apr 3, 2016
 *      Author: rice
 */

#include "environment.h"
#include "logger.h"
#include <cstdlib>
#include <unistd.h>
#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>


static bool isdir(const std::string &dirpath) {
  struct stat sb;
  if (stat(dirpath.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
    return true;
  }

  return false;
}

void Environment::setDatapath() {
  std::string datapath = getcwd() + "/rsf/";

  if (!isdir(datapath)) {
    if( mkdir(datapath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0) {
      ERROR() << "mkdir faild";
      exit(0);
    }

  }

  setenv("DATAPATH", datapath.c_str(), 1);
}

std::string Environment::getcwd() {
  char cwd[1024];
  if (::getcwd(cwd, sizeof(cwd)) != NULL)
    std::fprintf(stdout, "Current working dir: %s\n", cwd);
  else
    std::perror("getcwd() error");

  return cwd;
}

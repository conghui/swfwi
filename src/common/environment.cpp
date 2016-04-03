/*
 * Environment.cpp
 *
 *  Created on: Apr 3, 2016
 *      Author: rice
 */

#include "environment.h"
#include <cstdlib>
#include <unistd.h>
#include <cstdio>


void Environment::setDatapath() {
  std::string datapath = getcwd() + "/rsf/";
  std::string mkdircmd = std::string("mkdir -p ") + datapath;
  std::system(mkdircmd.c_str());
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

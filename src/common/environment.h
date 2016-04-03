/*
 * Environment.h
 *
 *  Created on: Apr 3, 2016
 *      Author: rice
 */

#ifndef SRC_COMMON_ENVIRONMENT_H_
#define SRC_COMMON_ENVIRONMENT_H_

#include <string>

class Environment {
public:
  static void setDatapath();

private:
  static std::string getcwd();
};

#endif /* SRC_COMMON_ENVIRONMENT_H_ */

/*
 * random-code.h
 *
 *  Created on: Feb 27, 2016
 *      Author: rice
 */

#ifndef SRC_FWI_RANDOM_CODE_H_
#define SRC_FWI_RANDOM_CODE_H_

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <vector>

// This is a typedef for a random number generator.
typedef boost::minstd_rand base_generator_type;
typedef boost::uniform_int<> distribution_type;
typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;

class RandomCodes {
public:
  explicit RandomCodes(int seed);

public:
  std::vector<int> genPlus1Minus1(int nshots);
  int nextRand();

private:
  base_generator_type generator;
  gen_type codes_gen;
  boost::generator_iterator<gen_type> codes;
};
#endif /* SRC_FWI_RANDOM_CODE_H_ */

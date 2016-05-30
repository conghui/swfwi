#include "random-code.h"
RandomCodes::RandomCodes(int seed) :
  generator(seed), codes_gen(generator, distribution_type(0, 1)), codes(&codes_gen)
{

}

/**
 * return -1 or +1
 */
int RandomCodes::nextRand() {
  return *codes++ * 2 - 1;
}

std::vector<int> RandomCodes::genPlus1Minus1(int nshots) {
  std::vector<int> codes(nshots);
  for (int i = 0; i < nshots; i++) {
    codes[i] = nextRand();
  }

  return codes;
}


/*
int main()
{
  int nshots = 20;
  int seed = 1;
  RandomCodes r(seed);
  std::vector<int> essfwi_codes = r.genPlus1Minus1(nshots);

  for (unsigned i = 0; i < essfwi_codes.size(); i++) {
    std::cout << essfwi_codes[i] << " ";
  }
  std::cout << std::endl;

  RandomCodes r2(seed);
  std::vector<int> essfwi_codes_2 = r2.genPlus1Minus1(nshots);

  for (unsigned i = 0; i < essfwi_codes_2.size(); i++) {
    std::cout << essfwi_codes_2[i] << " ";
  }
  std::cout << std::endl;

  //////////////////////// another instance use different seeds /////////////////
  seed = 2;
  RandomCodes rr(seed);
  std::vector<int> enkf_codes = rr.genPlus1Minus1(nshots);

  for (unsigned i = 0; i < enkf_codes.size(); i++) {
    std::cout << enkf_codes[i] << " ";
  }
  std::cout << std::endl;

  return 0;
}
*/

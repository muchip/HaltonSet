#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "HaltonSet.hpp"

#define ZETA 3.

double testfunc(int dim, const double* x) {
  double retval = 0.6;

  for (int k = 1; k <= dim; ++k) retval += 0.2 * std::pow(k, -ZETA) * x[k - 1];
  return 1. / retval;
}

double testfunc2(int dim, const double* x) {
  double retval = 1.;
  double cst = 1. / 6  + 1.;
  for (int k = 0; k < dim; ++k) retval *= x[k] * (x[k] - 1) + cst;
  return retval;
}

double pairwisePlus(double* pts, int length) {
  if (length < 3) {
    double sum = 0;
    for (int i = 0; i < length; ++i) sum += pts[i];
    return sum;
  } else {
    int m = length / 2;
    return pairwisePlus(pts, m) + pairwisePlus(pts + m, length - m);
  }
}

/**   \brief to come
*
*            to write everything into a binary file and read it by
*            Matlab, use
*            std::ofstream myfile;
*            myfile.open("HaltonSeq.dat", std::ios::out | std::ios::binary);
*            const Eigen::VectorXd &pHV = HS.get_HaltonVector();
*            myfile.write((const char *)&(pHV(0)), dim * sizeof(double));
*            can be read to matlab according to
*            A = fread(fopen('HaltonSeq.dat','rb'),[dim length],'double');
*/
int main(int argc, char* argv[]) {
  struct timeval start;
  struct timeval stop;
  double dtime = 0;
  int maxq = 12;
  int dim = 100000;
  double result = 0;
  HaltonSet HS(dim, 0);
  std::fstream myfile;
  myfile.open("HaltonSeq.dat", std::ios::out | std::ios::binary);
  for (int q = maxq; q <= maxq; ++q) {
    gettimeofday(&start, NULL);
    HS.reset();
    result = 0;
    for (int i = 0; i < (1 << q); ++i) {
      auto pHV = HS.get_HaltonVector();
      myfile.write((const char *)pHV.data(), dim * sizeof(double));
      //result += testfunc2(pHV.size(), pHV.data());
      HS.next();
    }
    gettimeofday(&stop, NULL);
    //result /= (1 << q);
    //results(q) = pairwisePlus(intVals.data(), intVals.size());
    //if (q < maxq) 
     //std::cout << "error: " << std::abs(1-result) << " numPts: " << (1 << q) << "\n";
    myfile.close();
    dtime = stop.tv_sec - start.tv_sec + 1e-6 * (stop.tv_usec - start.tv_usec);
    std::cout << "time taken:" << dtime << std::endl;
  }

  (void)argc;
  (void)argv;

  return 0;
}

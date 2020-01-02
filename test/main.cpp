#include <sys/time.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "../HaltonSet.hpp"
#include "primes.hpp"

#define ZETA 4.

double testfunc(int dim, const double* x) {
  double retval = 0.6;
  for (int k = 1; k <= dim; ++k)
    retval += 0.2 * std::pow(k, -ZETA) * (2 * x[k - 1] - 1);
  return 1. / retval;
}

double testfunc2(int dim, const double* x) {
  double retval = 1.;
  double cst = 1. / 6 + 1.;
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
  //////////////////////////////////////////////////////////////////////////////
  // test prime numbers
  //////////////////////////////////////////////////////////////////////////////
  {
    gettimeofday(&start, NULL);
    HaltonSet<10000, 0> HS;
    auto p_Hvec = HS.get_primes();
    gettimeofday(&stop, NULL);
    for (auto i = 0; i < 10000; ++i)
      assert(p_Hvec[i] == primes[i] && "primes appear to be false");
    std::cout << "Prime test passed. ";
    dtime = stop.tv_sec - start.tv_sec + 1e-6 * (stop.tv_usec - start.tv_usec);
    std::cout << "Time taken: " << dtime << " sec." << std::endl;
  }
  //////////////////////////////////////////////////////////////////////////////
  // test convergence
  //////////////////////////////////////////////////////////////////////////////
  {
    constexpr unsigned int dim = 1000;
    constexpr unsigned int maxq = 20;
    double result = 0;
    double reference = 0;
    HaltonSet<dim, 100> HS;
    // compute reference
    gettimeofday(&start, NULL);
    HS.reset();
    for (auto i = 0; i < (1 << (maxq + 1)); ++i) {
      auto pHV = HS.get_HaltonVector();
      reference += testfunc(pHV.size(), pHV.data());
      HS.next();
    }
    reference /= (1 << (maxq + 1));

    for (auto q = 0; q <= maxq; ++q) {
      HS.reset();
      result = 0;
      for (auto i = 0; i < (1 << q); ++i) {
        auto pHV = HS.get_HaltonVector();
        result += testfunc(pHV.size(), pHV.data());
        HS.next();
      }
      result /= (1 << q);
      auto error = std::abs(reference - result) / std::abs(reference);
      std::cout << std::setprecision(10) << "error: " << error
                << "\t numPts: 2^" << q << "\t dim: " << dim << std::endl;
    }
    gettimeofday(&stop, NULL);
    dtime = stop.tv_sec - start.tv_sec + 1e-6 * (stop.tv_usec - start.tv_usec);
    std::cout << "Time taken: " << dtime << " sec." << std::endl;
  }
#if 0
  //////////////////////////////////////////////////////////////////////////////
  // screen output
  //////////////////////////////////////////////////////////////////////////////
  {
    HaltonSet<6, 0> HS;
    for (auto i = 0; i < 20; ++i) {
      auto bla = HS.get_HaltonVector();
      for (auto j = 0; j < 6; ++j) std::cout << bla[j] << " ";
      std::cout << std::endl;
      HS.next();
    }
  }
#endif
  //////////////////////////////////////////////////////////////////////////////
  // write Halton set sample to file
  //////////////////////////////////////////////////////////////////////////////
  {
    constexpr unsigned int dim = 10000;
    constexpr unsigned int npts = 1e5;
    std::fstream myfile;
    HaltonSet<dim, 0> HS;
    std::cout << "writing " << npts
              << " points of the Halton sequence in dimension " << dim
              << " to file" << std::endl;
    HS.reset();
    myfile.open("HaltonSeq.dat", std::ios::out | std::ios::binary);
    for (auto i = 0; i < npts; ++i) {
      auto pHV = HS.get_HaltonVector();
      myfile.write((const char*)pHV.data(), dim * sizeof(double));
      HS.next();
    }
    myfile.close();
  }
  (void)argc;
  (void)argv;
  return 0;
}

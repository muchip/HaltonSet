#ifndef HALTONSET_CLASS_H_
#define HALTONSET_CLASS_H_

#include <array>
#include <limits>
#include <vector>

template <typename T, unsigned int M, unsigned int S, unsigned int B>
class HaltonSetBase {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //  constructors
  //////////////////////////////////////////////////////////////////////////////
  HaltonSetBase() {
    init_primes();
    init_HaltonVector();
  }
  //////////////////////////////////////////////////////////////////////////////
  //  core routine
  //  -updates bAdic representation of the current index
  //  -evaluates radical inverse of each representation by Horner's method
  //////////////////////////////////////////////////////////////////////////////
  void next() {
    long double radInverse = 0;
    long double bInv = 0;
    unsigned digit = 0;

    for (auto i = 0; i < M; ++i) {
      // update index by adding 1 to each of the _M bAdic representations
      for (digit = 0; digit < B; ++digit) {
        // check if increasing index will generate a carry, if not just
        // add it, otherwise handle
        if (bAdic_[i][digit] + 1 < primes_[i]) {
          ++bAdic_[i][digit];
          break;
        } else
          bAdic_[i][digit] = 0;
      }
      assert(digit < B && "overflow in bAdic representation");
      // update length of bAdic representation
      if (maxDigit_[i] < digit) maxDigit_[i] = digit;
      // evaluate the radical inverse by the Horner scheme
      bInv = (long double)1 / (long double)primes_[i];
      radInverse = bAdic_[i][maxDigit_[i]];
      for (int j = maxDigit_[i] - 1; j >= 0; --j)
        radInverse = bInv * radInverse + bAdic_[i][j];
      radInverse *= bInv;
      // store the current point to the interfacing vector
      HaltonVector_[i] = (T)radInverse;
    }
    return;
  }
  // method to reset the Halton vector to its initial state
  void reset(void) {
    HaltonVector_.resize(M);
    maxDigit_.resize(M);
    bAdic_.resize(M);
    // reset Halton vector
    std::fill(std::begin(HaltonVector_), std::end(HaltonVector_), 0);
    // reset maxDigit
    std::fill(std::begin(maxDigit_), std::end(maxDigit_), 0);
    // reset bAdic
    for (auto i = 0; i < M; ++i)
      std::fill(std::begin(bAdic_[i]), std::end(bAdic_[i]), 0);
    return;
  }
  //////////////////////////////////////////////////////////////////////////////
  //  getter
  //////////////////////////////////////////////////////////////////////////////
  const std::vector<T> &get_HaltonVector(void) const { return HaltonVector_; }
  const std::vector<unsigned int> &get_primes(void) const { return primes_; }
  //////////////////////////////////////////////////////////////////////////////
  //  private members
  //////////////////////////////////////////////////////////////////////////////
 private:
  // set up array of prime numbers
  void init_primes() {
    primes_.resize(M);
    unsigned int number_of_primes = 0;
    // treat 2 and 3 explicitly to shortcout afterwards
    unsigned int current_integer = 3;
    if (M > 0) {
      primes_[0] = 2;
      ++number_of_primes;
    }
    // regard that this is only called if M > 1 in this case 3 is the next prime
    // all primes following 2 have to be odd!
    while (number_of_primes < M) {
      // compute next prime number
      while (!isPrime(current_integer, number_of_primes)) current_integer += 2;
      // store prime new prime
      primes_[number_of_primes] = current_integer;
      ++number_of_primes;
      current_integer += 2;
    }
    return;
  }
  // compare with respect to the already found numbers if a number is prime
  bool isPrime(unsigned int test_number, int number_of_primes) {
    unsigned int upper = std::sqrt(test_number);
    for (auto i = 0; i < number_of_primes && primes_[i] <= upper; ++i)
      if (!(test_number % primes_[i])) return false;
    return true;
  }
  /** \brief initializes the HaltonVector with respect to the warm up defined
   *          in _skip
   */
  void init_HaltonVector(void) {
    reset();
    for (auto i = 0; i < S; ++i) next();
    return;
  }
  // member variables (using std::vector since std::array is stack allocated
  // and strongly limited in storage)
  std::vector<unsigned int> primes_;
  std::vector<std::array<unsigned int, B>> bAdic_;
  std::vector<unsigned int> maxDigit_;
  std::vector<T> HaltonVector_;
};

/**
 *    \brief Specialization with double Halton points and
 *           maximum base length 40bits, i.e. for b=2
 *           maximum index is 1099511627776-1 \approx 10^12
 **/
template <unsigned int M, unsigned int S>
using HaltonSet = HaltonSetBase<double, M, S, 40u>;

#endif

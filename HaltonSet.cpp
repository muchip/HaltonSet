#include "HaltonSet.hpp"

HaltonSet::HaltonSet(int M) : _skip(100), _M(M) {
  HaltonSet::init_primes();
  HaltonSet::init_HaltonVector();
}

HaltonSet::HaltonSet(int M, int skip) : _skip(skip), _M(M) {
  HaltonSet::init_primes();
  HaltonSet::init_HaltonVector();
}

const std::vector<double> &HaltonSet::get_HaltonVector(void) const {
  return _HaltonVector;
}

bool HaltonSet::isPrime(std::vector<int> &primes, int test) {
  int upper = std::sqrt(test);
  for (int i = 0; i <= upper; ++i)
    if (!(test % primes[i]))
      return false;

  return true;
}

void HaltonSet::init_primes(void) {
  int testVal = 0;
  _primes.clear();
  _primes.push_back(2);
  _primes.push_back(3);
  _primes.push_back(5);

  while (int(_primes.size()) < _M) {
    testVal = _primes.back() + 2;
    while (!isPrime(_primes, testVal))
      testVal += 2;
    _primes.push_back(testVal);
  }
}

void HaltonSet::init_HaltonVector(void) {
  // choose 40bits for _bAdic,
  // i.e. for b=2 maximum index is 1099511627776-1 \approx 10^12
  _bAdic.resize(_M);
  for (auto i = _bAdic.begin(); i != _bAdic.end(); ++i) {
    (*i).resize(40);
    std::fill((*i).begin(), (*i).end(), 0);
  }
  _HaltonVector.resize(_M);
  _maxDigit.resize(_M);

  std::fill(_HaltonVector.begin(), _HaltonVector.end(), 0);
  std::fill(_maxDigit.begin(), _maxDigit.end(), 0);

  for (int i = 0; i < _skip; ++i) HaltonSet::next();
}

void HaltonSet::next(void) {
  long double radInverse = 0;
  long double bInv = 0;
  int digit = 0;
  int next = 0;

  for (int i = 0; i < _M; ++i) {
    digit = 0;
    next = 1;
    // update index by adding 1 to each of the _M bAdic representations
    while (next) {
      if (_bAdic[i][digit] + 1 < _primes[i]) {
        ++_bAdic[i][digit];
        next = 0;
      } else {
        _bAdic[i][digit] = 0;
        ++digit;
        next = 1;
      }
    }
    // update length of bAdic representation
    if (_maxDigit[i] < digit) _maxDigit[i] = digit;
    // evaluate the radical inverse by the Horner scheme
    bInv = (long double)1 / (long double)_primes[i];
    radInverse = _bAdic[i][_maxDigit[i]];
    for (int j = _maxDigit[i] - 1; j >= 0; --j)
      radInverse = bInv * radInverse + _bAdic[i][j];
    radInverse *= bInv;
    // store the current point to the interfacing vector
    _HaltonVector[i] = (double)radInverse;
  }
}

void HaltonSet::reset(void) { HaltonSet::init_HaltonVector(); }


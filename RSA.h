#ifndef RSA_H
#define RSA_H
#include <cmath>
#include <thread>
#include <vector>
#include <boost\multiprecision\cpp_int.hpp>
#define MAX_THREADS 16
#define MAX_RANGE 10000000
#define millerRabinTimes 10
// #define MAX_RANGE 1000000000000000000000000000
using namespace boost::multiprecision;
// typedef unsigned long long int ulli;
typedef cpp_int ulli;
ulli hcf(ulli a, ulli b);
bool isPrime(ulli n);
ulli modInverse(ulli a, ulli m);
ulli Encrypt(ulli e, ulli n, ulli plainText);
ulli Decrypt(ulli d, ulli n, ulli cipherText);
class Rsa
{
public:
    Rsa() {}
    Rsa(ulli range1_, ulli range2_)
    /*
    range1_ and range2_ are the range of the prime numbers to be used for generating p and q.
    range1_ != range2_
    range1_ > 2
    range2_ > 2
    */
    {
        if (range1_ != range2_ && range1_ > 2 && range2_ > 2)
        {
            range1 = std::min(range1_, range2_);
            range2 = std::max(range1_, range2_);
        }
    }
    ulli CalNextQ();
    ulli CalNextP();
    ulli CalNextE();
    ulli CalNextD();
    bool CalPQED();

    ulli GetP() { return p; }
    ulli GetQ() { return q; }
    ulli GetE() { return e; }
    ulli GetD() { return d; }
    ulli GetN()
    {
        n = p * q;
        return n;
    }
    ulli GetPhi()
    {
        phi = (p - 1) * (q - 1);
        return phi;
    }

private:
    ulli range1 = -1, range2 = -1;
    ulli p = 0, q = 0, e = 0, d = 0, n = 0, phi = 0;
};
#include "RSA.cpp"
#endif
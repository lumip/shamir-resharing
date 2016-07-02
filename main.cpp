/**
  * shamir-resharing -- a prototype implementation of shamir secret sharing demonstrating resharing
  * Copyright (C) 2016  Lukas Prediger <lukas.prediger@rwth-aachen.de>
  *
  * This program is free software: you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation, either version 3 of the License, or
  * (at your option) any later version.
  *
  * This program is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU General Public License for more details.
  *
  * You should have received a copy of the GNU General Public License
  * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/


#include <stdlib.h>
#include <time.h>
#include <gmpxx.h>
#include <gmp.h>
#include <vector>
#include <iostream>

typedef std::pair<mpz_class, mpz_class> share_t;


mpz_class computeRandomPrime(gmp_randclass& randomGenerator, unsigned int security)
{
    mpz_class prime = randomGenerator.get_z_bits(security - 1);
    mpz_setbit(prime.get_mpz_t(), security - 1);
    mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
    return prime;
}

mpz_class invertMultiplicative(const mpz_class& x, const mpz_class& p)
{
    if (x == p - 1 || x == 1)
        return x;
    mpz_t gcd, s, t;
    mpz_init(gcd);
    mpz_init(s);
    mpz_init(t);
    mpz_gcdext(gcd, s, t, x.get_mpz_t(), p.get_mpz_t());
    mpz_class inverse(s);
    if (inverse < 0)
        inverse += p;
    return inverse;
}

mpz_class invertAdditive(const mpz_class& x, const mpz_class& p)
{
    return p - x;
}

std::vector<mpz_class> getRandomPolynomial(unsigned int degree, const mpz_class& p, gmp_randclass& randomGenerator, unsigned int security)
{
    std::vector<mpz_class> coefficients;
    coefficients.resize(degree + 1);
    for (size_t i = 1; i < coefficients.size(); ++i)
    {
        coefficients[i] = randomGenerator.get_z_bits(security) % p;
    }
    return coefficients;
}

mpz_class evaluatePolynomial(unsigned int x, const std::vector<mpz_class>& coefficients, const mpz_class& p)
{
    mpz_class result(0);
    mpz_class expX(1);
    for (size_t i = 0; i < coefficients.size(); ++i)
    {
        result = (result + (coefficients[i] * expX) % p) % p;
        expX = (expX * x) % p;
    }
    return result;
}


std::vector< share_t > share(const mpz_class& secret, unsigned int t, unsigned int n, const mpz_class& p, gmp_randclass& randomGenerator, unsigned int security)
{
    std::vector<mpz_class> coefficients = getRandomPolynomial(t - 1, p, randomGenerator, security);
    coefficients[0] = secret;
    std::vector< share_t > shares;
    shares.resize(n);
    for (size_t i = 0; i < shares.size(); ++i)
    {
        shares[i].first = mpz_class(i + 1);
        shares[i].second = evaluatePolynomial(i + 1, coefficients, p);
    }
    return shares;
}

mpz_class evaluateLagrangePolynomial(const mpz_class& x, const mpz_class& x_j, const std::vector< mpz_class >& x_ms, const mpz_class& p)
{
    mpz_class result(1);
    for (size_t i = 0; i < x_ms.size(); ++i)
    {
        const mpz_class& x_m = x_ms[i];
        if (x_j == x_m)
            continue;
        const mpz_class x_m_inv = invertAdditive(x_m, p);
        result = (result * ((x + x_m_inv) % p) * invertMultiplicative((x_j + x_m_inv) % p , p)) % p;
    }
    return result;
}

mpz_class recover(const std::vector< share_t >& shares, const mpz_class& p)
{
    const mpz_class zero(0);
    std::vector<mpz_class> indices;
    indices.resize(shares.size());
    for (size_t i = 0; i < shares.size(); ++i)
    {
        indices[i] = shares[i].first;
    }
    mpz_class secret(0);
    for (size_t i = 0; i < shares.size(); ++i)
    {
        const share_t& share = shares[i];
        secret = (secret + evaluateLagrangePolynomial(zero, share.first, indices, p) * share.second) % p;
    }
    return secret;
}


int main()
{
    const unsigned int security = 512;

    gmp_randclass randomGenerator(gmp_randinit_default);
    randomGenerator.seed(time(0));
    //randomGenerator.seed(72535); // yields p = 12742741321067026999

    mpz_class p = computeRandomPrime(randomGenerator, security);
    //mpz_class p(23);
    gmp_printf("%Zu has length %d and is prime? %u\n", p.get_mpz_t(), mpz_sizeinbase(p.get_mpz_t(), 2), mpz_probab_prime_p(p.get_mpz_t(), 25));


    /*std::vector<mpz_class> coefficients;
    coefficients.push_back(mpz_class(4));
    coefficients.push_back(mpz_class(2345));
    coefficients.push_back(mpz_class(17));
    std::vector<mpz_class> samples;
    samples.push_back(mpz_class(0));
    samples.push_back(mpz_class(1));
    samples.push_back(mpz_class(2));
    samples.push_back(mpz_class(3));
    samples.push_back(mpz_class(17));
    for (size_t i = 0; i < samples.size(); ++i)
    {
        std::cout << "f( " << samples[i].get_str() << " ) = " << evaluatePolynomial(samples[i].get_ui(), coefficients, p).get_str() << std::endl;
    }*/

    mpz_class secret(20160207);
    secret = secret % p;
    std::cout << "the secret is " << secret.get_str() << std::endl;

    std::vector< share_t > shares = share(secret, 4, 30, p, randomGenerator, security);
    for (size_t i = 0; i < shares.size(); ++i)
    {
        gmp_printf("share %Zu is %Zx\n", shares[i].first.get_mpz_t(), shares[i].second.get_mpz_t());
        //std::cout << "share " << shares[i].first.get_str() << " is " << shares[i].second.get_str() << std::endl;
    }
    std::vector< share_t > recover_shares;
    recover_shares.resize(4);
    std::copy(shares.begin() + 1, shares.begin() + 5, recover_shares.begin());

    std::vector<mpz_class> indices;
    indices.resize(recover_shares.size());
    for (unsigned int i = 0; i < recover_shares.size(); ++i)
    {
        indices[i] = recover_shares[i].first;
    }

    mpz_class recovered = recover(recover_shares, p);
    std::cout << "have recovered " << recovered.get_str() << std::endl;

    return 0;
}

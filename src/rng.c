/*
 * xoshiro256++ and SplitMix64 are adapted from the public-domain reference
 * implementations by David Blackman and Sebastiano Vigna:
 * https://prng.di.unimi.it/
 */
#include "fastLISA.h"

// Rotate left helper for xoshiro256++
static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

// splitmix64 PRNG generator:
// Used to initialize/seed the 256-bit state arrays of the xoshiro256++ generator
// from a single 64-bit seed value. Highly robust for seeding.
uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

// xoshiro256++ 1.0 PRNG generator:
// This is a state-of-the-art fast pseudorandom number generator with 256 bits of state.
// Period is 2^256 - 1. Extremely fast, lightweight, and thread-independent,
// making it ideal for parallelized OpenMP simulation loops.
uint64_t xoshiro256plusplus(uint64_t s[4]) {
    const uint64_t result = rotl(s[0] + s[3], 23) + s[0];

    const uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;

    s[3] = rotl(s[3], 45);

    return result;
}

// Draw a random integer uniformly distributed in [0, n - 1].
// Plain modulo (val % n) is biased whenever n does not divide 2^64: the low
// residues occur slightly more often. We remove that bias with rejection
// sampling (OpenBSD/Java style): reject any draw falling in the first
// (2^64 mod n) values so the accepted range is an exact multiple of n.
// Callers always pass n >= 1, so no zero guard is needed. Rejection consumes
// extra RNG outputs only on the rare reject, keeping the stream deterministic.
// Because each kernel re-seeds the generator per observation, results are
// reproducible for any thread count (not just a fixed one).
int rng_int(uint64_t s[4], int n) {
    uint64_t range     = (uint64_t)n;
    uint64_t threshold = (-range) % range;   /* == 2^64 mod range */
    uint64_t r;
    do {
        r = xoshiro256plusplus(s);
    } while (r < threshold);
    return (int)(r % range);
}

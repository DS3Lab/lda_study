//
//  utils.h
//  cpp_test
//
//  Created by yulele on 16/9/6.
//  Copyright (c) 2016年 yulele. All rights reserved.
//

#ifndef cpp_test_utils_h
#define cpp_test_utils_h

#include <iostream>
#include <vector>
#include <sstream>
#include <random>
#include <sys/time.h>
#include "csr.h"
#include "xorshift.h"

namespace lda {

    void split(const std::string &s, char delim, std::vector<std::string> &elems);
    
    inline static std::mt19937 & urng() {
        static std::mt19937 u{};
        return u;
    }
    
    inline static double unif01() {
        static std::uniform_real_distribution<double> d(0.0, 1.0);
        return d(urng());
    }
    
    inline static int rand_int(int from, int to) {
        static std::uniform_int_distribution<> d{};
        using param_t = std::uniform_int_distribution<>::param_type;
        return d(urng(), param_t{from, to});
    }

    inline static long get_current_ms() {
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv, &tz);
        
        long timestamp = (long) tv.tv_sec * 1000 + tv.tv_usec / 1000;
        return timestamp;
    }

    inline static xorshift & gen() {
      static xorshift gen;
      return gen;
    }

    inline static double rand_double() {
      return gen().Rand32() * 1.0 / std::numeric_limits<uint32_t>::max();
    }

    inline static int rand_int(int max) {
      return gen().Rand32() % max;
    }
    
    inline static long murmurHash3(long x) {
        x ^= ((uint64_t) x) >> 33;
        x *= -49064778989728563L;
        x ^= ((uint64_t) x) >> 33;
        x *= -4265267296055464877L;
        x ^= ((uint64_t) x) >> 33;
        return x;
    }
    
    inline static int murmur_hash3(int x) {
        x ^= ((uint32_t)x) >> 16;
        x *= -2048144789;
        x ^= ((uint32_t)x) >> 13;
        x *= -1028477387;
        x ^= ((uint32_t)x) >> 16;
        return x;
    }
    
    inline static int next_power_of_two(int x) {
        if(x == 0) {
            return 1;
        } else {
            --x;
            x |= x >> 1;
            x |= x >> 2;
            x |= x >> 4;
            x |= x >> 8;
            return (x | x >> 16) + 1;
        }
    }
    
    double llhw(int D, int V, int K, int** nwk, int** ndk, int* nk,
                       double alpha, double beta, double vbeta, corpus* cp);
    
    
}

#endif

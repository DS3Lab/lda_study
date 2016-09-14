//
//  traverse_map.h
//  cpp_test
//
//  Created by yulele on 16/9/9.
//  Copyright (c) 2016年 yulele. All rights reserved.
//

#ifndef __cpp_test__traverse_map__
#define __cpp_test__traverse_map__

#include <stdio.h>
#include <iostream>
#include "utils.h"

namespace lda {

class traverse_map {
public:
    uint16_t * key;
    uint16_t * value;
    bool * used;
    uint16_t * idx;
    uint16_t * poss;
    int mask;
    int n;
    int size;
    
    traverse_map() {
        
    }
    
    void init(uint16_t expected) {
        n = next_power_of_two(expected);
        mask = n - 1;
        key = new uint16_t[n];
        value = new uint16_t[n];
        used = new bool[n];
        idx = new uint16_t[n];
        poss = new uint16_t[n];
        size = 0;
    }
    
    ~traverse_map() {
        delete key;
        delete value;
        delete used;
        delete idx;
        delete poss;
    }
    
    uint16_t get(uint32_t k) {
        return get((uint16_t) k);
    }
    
    uint16_t get(uint16_t k) {
        // The starting point
        int pos = ( murmur_hash3(k)) & mask;
        
        // There's always an unused entry.
        int cnt = 0;
        while (used[pos]) {
            if (key[pos] == k) {
                return value[pos];
            }
            pos = (pos + 1) & mask;
            cnt ++;
            
            if (cnt > n) {
                rehash();
                return get(k);
            }
        }
        return 0;
    }
    
    void put(const int k, const int v) {
        put((uint16_t) k, (uint16_t) v);
    }
    
    void put(uint16_t k, uint16_t v) {
        if (v == 0)
            return;
        
        // The starting point
        uint16_t pos = ( murmur_hash3(k)) & mask;
        
        // There's always an unused entry.
        while (used[pos]) {
            if (key[pos] == k) {
                value[pos] = v;
                return;
            }
            pos = (pos + 1) & mask;
        }
        
        used[pos] = true;
        key[pos] = k;
        value[pos] = v;
        idx[size] = pos;
        poss[pos] = size;
        size ++;
    }
    
    void rehash() {
        
        uint16_t * kkey = key;
        uint16_t * vvalue = value;
        
        key = new uint16_t[n];
        value = new uint16_t[n];
        
        memset(used, 0, n * sizeof(bool));
        
        int temp = size;
        size = 0;
        
        for (int i = 0; i < temp; i ++) {
            uint16_t k = kkey[idx[i]];
            uint16_t v = vvalue[idx[i]];
            put(k, v);
        }
        
        delete kkey;
        delete vvalue;
        
    }
    
    void dec(const int k) {
        dec((uint16_t) k);
    }
    
    void dec(uint16_t k) {
        int pos = ( murmur_hash3(k)) & mask;
        
        while (used[pos]) {
            if (key[pos] == k) {
                value[pos] --;
                if (value[pos] == 0) {
                    size --;
                    idx[poss[pos]] = idx[size];
                    poss[idx[size]] = poss[pos];
                }
                return;
            }
            
            pos = (pos + 1) & mask;
        }
    }
    
    void inc(int k) {
        inc((uint16_t) k);
    }
    
    void inc(uint16_t k) {
        uint16_t pos = ( murmur_hash3(k)) & mask;
        
        int cnt = 0;
        while (used[pos]) {
            if (key[pos] == k) {
                value[pos] ++;
                if (value[pos] == 1) {
                    idx[size] = pos;
                    poss[pos] = size;
                    size ++;
                }
                
                return;
            }
            
            cnt ++;
            if (cnt > n) {
                rehash();
                inc(k);
                return;
            }
            pos = (pos + 1) & mask;
        }
        
        key[pos] = k;
        value[pos] = 1;
        used[pos] = true;
        idx[size] = pos;
        poss[pos] = size;
        size ++;
    }
};
    
}

#endif /* defined(__cpp_test__traverse_map__) */

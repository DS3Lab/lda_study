//
//  ftreelda.cpp
//  cpp_test
//
//  Created by yulele on 16/9/8.
//  Copyright (c) 2016年 yulele. All rights reserved.
//

#include "ftreelda.h"
#include "utils.h"
#include <unordered_map>
#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include "traverse_map.h"

void handler(int sig) {
    void *array[10];
    size_t size;
    
    // get void*'s for all entries on the stack
    size = backtrace(array, 10);
    
    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, (int)size, STDERR_FILENO);
    exit(1);
}



namespace lda {
    
    void ftreelda::init() {
        uint32_t did;
        uint16_t tt;
        
        memset(_nk, 0, _K * sizeof(int));
        
        for (int d = 0; d < _D; d ++)
            _ndk[d].init(_adj->_dlens[d]);
        
        for (int w = 0; w < _V; w ++) {
            for (int wi = _ws[w]; wi < _ws[w + 1]; wi ++) {
                did = _dids[wi];
                tt  = rand_int(_K);
                _ndk[did].inc(tt);
                _nk[tt] ++;
                _topics[wi] = tt;
            }
        }
    }
    
    void ftreelda::train_one_word(int w) {
        
        int N = _ws[w + 1] - _ws[w];
        
        if (N == 0)
            return;
        
        memset(_wk, 0, _K * sizeof(int));
        
        // build wk
        for (int wi = _ws[w]; wi < _ws[w + 1]; wi ++) {
            uint16_t k = _topics[wi];
            _wk[k] ++;
        }
        
        // build ftree
        for (int k = 0; k < _K; k ++) {
            _p[k] = (_wk[k] + _beta) / (_nk[k] + _vbeta);
        }
        
        _tree.build(_p);
        
        // for each token of this word
        double psum, u;
        uint32_t did;
        uint16_t tt, ttt;
        
        for (int wi = _ws[w]; wi < _ws[w + 1]; wi ++) {
            tt  = _topics[wi];
            did = _dids[wi];
            
            _wk[tt] --;
            _nk[tt] --;
            _ndk[did].dec(tt);
            
            _tree.update(tt, (_wk[tt] + _beta) / (_nk[tt] + _vbeta));
            
            psum = 0.0;
            
            for (int i = 0; i < _ndk[did].size; i ++) {
                uint16_t kk = _ndk[did].key[_ndk[did].idx[i]];
                uint16_t vv = _ndk[did].value[_ndk[did].idx[i]];
                psum += vv * _tree.get(kk);
                _p[i] = psum;
                _t_idx[i] = kk;
            }
            
            u = rand_double() * (psum + _tree.first() * _alpha);
            
            ttt = -1;
            if (u < psum) {
                double u1 = rand_double() * psum;
                for (int i = 0; i < _ndk[did].size; i ++)
                    if (u1 < _p[i]) {
                        ttt = _t_idx[i];
                        break;
                    }
            } else {
                ttt = _tree.sample(rand_double() * _tree.first());
            }
            
            _wk[ttt] ++;
            _nk[ttt] ++;
            _ndk[did].inc(ttt);
            _topics[wi] = ttt;
            _tree.update(ttt, (_wk[ttt] + _beta) / (_nk[ttt] + _vbeta));
        }
        
        double lgamma_beta = lgamma(_beta);
        
        for (int k = 0; k < _K; k ++) {
            if (_wk[k] > 0)
                _ll += lgamma(_beta + _wk[k]) - lgamma_beta;
        }
    }
    
    void ftreelda::train_one_iteration(int it) {
        
        long start, end;
        
        start = get_current_ms();
        _ll = 0;
        for (int w = 0; w < _V; w ++)
            train_one_word(w);
        
        end = get_current_ms();
        long train_tt = end - start;
        
        start = get_current_ms();
        
        _ll += _K * lgamma(_beta * _V);
        for (int k = 0; k < _K; k ++)
            _ll -= lgamma(_beta * _V + _nk[k]);
        
        double lgamma_alpha = lgamma(_alpha);
        
        for (int d = 0; d < _D; d ++) {
            for (int i = 0; i < _ndk[d].size; i ++) {
                uint16_t vv = _ndk[d].value[_ndk[d].idx[i]];
                _ll += lgamma(_alpha + vv) - lgamma_alpha;
            }
            _ll += lgamma(_alpha * _K) - lgamma(_alpha * _K + _adj->_dlens[d]);
        }
        
        end = get_current_ms();
        
        long eval_tt = end - start;
        
        printf("train_tt=%ld eval_tt=%ld ll=%f\n", train_tt, eval_tt, _ll);
    }
}

void sample_ftree_lda_one_iteration(int D, int V, int K, uint32_t *ws, uint32_t* dids, uint16_t* topics, std::unordered_map<uint16_t, uint16_t> * ndk) {
    
    
    
}

void test_traverse_map() {
    lda::traverse_map map;
    map.init(50);
    
    std::unordered_map<int, int> tmap(50);
    
    for (int i = 0; i < 100; i ++) {
        int k = lda::rand_int(50);
        map.inc(k);
        tmap[k] ++;
    }
    
    for (int i = 0; i < map.size; i ++) {
        uint16_t kk = map.key[map.idx[i]];
        uint16_t vv = map.value[map.idx[i]];
        
        if (tmap[kk] != vv) {
            printf("Error\n");
        }
    }
    
    printf("done\n");
    
}

int main(int argc, const char * argv[]) {
    
    int V = 12420;
    lda::adj_list* adj = new lda::adj_list(V);
    
    adj->load("data/nips.train");
    printf("load finish\n");
    lda::ftreelda * lda = new lda::ftreelda(adj->_D, V, 1024, 0.1, 0.1, adj);
    lda->init();
    for (int i = 0; i < 100; i ++)
        lda->train_one_iteration(i);
    
//    for (int i = 0; i < 1024; i ++)
//        printf("%d %d\n", i, lda::murmurHash3(i));
    
//    test_traverse_map();
    
}
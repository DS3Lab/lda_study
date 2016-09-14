//
//  warplda.cpp
//  cpp_test
//
//  Created by yulele on 16/9/7.
//  Copyright (c) 2016年 yulele. All rights reserved.
//

#include <random>
#include <string.h>
#include <math.h>
#include "warplda.h"
#include "utils.h"
#include "xorshift.h"

using namespace std;

namespace lda {
    
    void warplda::init() {
        int* ws = _mat->_ws;
        dv_mat::entry* data = _mat->_data;
        for (int w = 0; w < _V; w ++) {
            for (int wi = ws[w]; wi < ws[w + 1]; wi ++) {
                uint16_t tt = (uint16_t) _gen.rand_int(_K);
                data[wi].tt = tt;
                _nk[tt] ++;
                for (int m = 0; m < MH_STEPS; m ++)
                    data[wi].mh[m] = tt;
            }
        }
    }
    
    void warplda::visit_by_row() {
        uint16_t tt, ttt;
        
        int* ds = _mat->_ds;
        int* widx = _mat->_widx;
        dv_mat::entry* data = _mat->_data;
        
        int N;
        float p;
        int idx;
        
        float lgamma_alpha = lgamma(_alpha);
        
        for (int d = 0; d < _D; d ++) {
            N = ds[d + 1] - ds[d];
            
            // compute dk for doc on the fly
            memset(_tk, 0, _K * sizeof(int));
            memset(_used, 0, _K * sizeof(uint8_t));
            
            idx = 0;
            for (int di = ds[d]; di < ds[d + 1]; di ++) {
                uint16_t k = data[widx[di]].tt;
                _tk[k] ++;
                if (_used[k] == 0) {
                    _used[k] = 1;
                    _nidx[idx ++] = k;
                }
            }
            
            for (int i = 0; i < idx; i ++)
                _ll += lgamma(_alpha + _tk[_nidx[i]]) - lgamma_alpha;
            
            _ll -= lgamma(_alpha * _K + N) - lgamma(_alpha * _K);
            
            // traverse doc first time, accept word proposal
            for (int di = ds[d]; di < ds[d + 1]; di ++) {
                tt = data[widx[di]].tt;
                _nk[tt] --;
                _tk[tt] --;

                for (int m = 0; m < MH_STEPS; m ++) {
                    ttt = data[widx[di]].mh[m];
                    p = ((_tk[ttt] + _alpha) * (_nk[tt] + _vbeta))
                                    / ((_tk[tt] + _alpha) * (_nk[ttt] + _vbeta));
                    
                    if (_gen.rand_double() < p) {
                        tt = ttt;
                    }
                }
                
                _nk[tt] ++;
                _tk[tt] ++;
                data[widx[di]].tt = tt;
            }
            
            // traverse doc second time, generate doc proposal
            p = (_alpha * _K) / (_alpha * _K + N);
            for (int di = ds[d]; di < ds[d + 1]; di ++) {
                for (int m = 0; m < MH_STEPS; m ++) {
                    if (_gen.rand_double() < p)
                        data[widx[di]].mh[m] = (uint16_t) _gen.rand_int(_K);
                    else
                        data[widx[di]].mh[m] = data[widx[ds[d] + _gen.rand_int(N)]].tt;
                }
            }
        }
    }
    
    void warplda::visit_by_col() {
        uint16_t tt, ttt;
        int* ws = _mat->_ws;
        dv_mat::entry* data = _mat->_data;
        
        
        int N;
        float p;
        int idx;
        
        float lgamma_beta = lgamma(_beta);
        
        for (int w = 0; w < _V; w ++) {
            // compute wk for word on the fly
            
            memset(_tk, 0, _K * sizeof(int));
            memset(_used, 0, _K * sizeof(uint8_t));
            
            N = ws[w + 1] - ws[w];
            idx = 0;
            
            for (int wi = ws[w]; wi < ws[w + 1]; wi ++) {
                uint16_t k = data[wi].tt;
                _tk[k] ++;
                if (_used[k] == 0) {
                    _nidx[idx ++] = k;
                    _used[k] = 1;
                }
            }
            
            for (int i = 0; i < idx; i ++)
                _ll += lgamma(_beta + _tk[_nidx[i]]) - lgamma_beta;
            
            
            // traverse w first time, accept doc proposal
            for (int wi = ws[w]; wi < ws[w + 1]; wi ++) {
                tt = data[wi].tt;
                _tk[tt] --;
                _nk[tt] --;

                
                for (int m = 0; m < MH_STEPS; m ++) {
                    ttt = data[wi].mh[m];
                    
                    float b = _tk[tt]  + _beta;
                    float d = _nk[tt]  + _vbeta;
                    float a = _tk[ttt] + _beta;
                    float c = _nk[ttt] + _vbeta;
                    
                    if (rand_double() < (a * d) / (b * c)) {
                        tt = ttt;
                    }
                }
                
                _tk[tt] ++;
                _nk[tt] ++;
                data[wi].tt = tt;
            }
            
            // traverse w second time, compute word proposal
            
            p = (_K * _beta) / (_K * _beta + N);
            
            for (int wi = ws[w]; wi < ws[w + 1]; wi ++) {
                for (int m = 0; m < MH_STEPS; m ++) {
                    if (_gen.rand_double() < p) {
                        data[wi].mh[m] = (uint16_t) _gen.rand_int(_K);
                    } else {
                        data[wi].mh[m] = data[ws[w] + _gen.rand_int(N)].tt;
                    }
                }
            }
        }
        
        for (int k = 0; k < _K; k ++) {
            if (_nk[k] > 0)
              _ll += - lgamma(_nk[k] + _vbeta) + lgamma(_vbeta);
        }
    }

    void warplda::build_nk() {
    
    }
    
    void warplda::train_one_iter(int iter) {
        _ll = 0;

        long st, sst;
        long col_tt, row_tt, iter_tt;


        st = get_current_ms();
        sst= st;
        visit_by_col();
        col_tt = get_current_ms() - st;
      
        st = get_current_ms();
        visit_by_row();
        row_tt = get_current_ms() - st;
        
        iter_tt = get_current_ms() - sst;
        
        printf("iter=%d iter_tt=%ld col_tt=%ld row_tt=%ld llhw=%f %02f Mtokens/s \n",
            iter, iter_tt, col_tt, row_tt, _ll, (double) (_mat->_N)/ ((double)iter_tt/1e3) / 1e6);
    }
    
}

void sample_warp_lda_one_iteration() {
    
}


using namespace std;

void warplda_nips() {
    //    string path = "/Users/yulele/workspace/cpp_test/nips.train";
    string path = "data/nips.train";
    int V = 12420;
    int K = 1024;
    float alpha = 0.1;
    float beta  = 0.1;
    lda::dv_mat* p_mat = new lda::dv_mat(path, V);
    
    lda::warplda* lda = new lda::warplda(p_mat->_D, V, K,
                                         alpha, beta, p_mat);
    lda->init();
    for (int i = 0; i < 200; i ++)
        lda->train_one_iter(i);
    
    delete lda;
    
}

int main(int argc, const char * argv[]) {
    //    cgslda_nips();
    //    check_dv_mat();
    warplda_nips();
    return 0;
}

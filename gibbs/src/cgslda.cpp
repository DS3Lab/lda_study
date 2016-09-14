//
//  cgslda.cpp
//  cpp_test
//
//  Created by yulele on 16/9/7.
//  Copyright (c) 2016年 yulele. All rights reserved.
//

#include "cgslda.h"
#include "utils.h"

using namespace std;

namespace lda {
    
    void cgslda::init() {
        printf("D=%d V=%d K=%d corpus size=%ld\n",
               _D, _V, _K, _cp->_docs.size());
        
        document * p_doc;
        int len, did, wid, tt;
        for (int d = 0; d < _D; d ++) {
            p_doc = _cp->_docs[d];
            len = p_doc->_n;
            did = p_doc->_did;
            p_doc->_topics = new short[len];
            for (int i = 0; i < len; i ++) {
                wid = p_doc->_wids[i];
                tt  = rand_int(_K);
                _nwk[wid][tt] ++;
                _ndk[did][tt] ++;
                _nk[tt] ++;
                p_doc->_topics[i] = tt;
            }
        }
    }
    
    void cgslda::train_one_iter(int iter) {
        
        long start;
        
        start = get_current_ms();
        for (int d = 0; d < _D; d ++) {
            train_one_doc(_cp->_docs[d]);
        }
        long train_tt = get_current_ms() - start;
        
        start = get_current_ms();
        double ll = llhw(_D, _V, _K, _nwk, _ndk, _nk, _alpha, _beta, _vbeta, _cp);
        long eval_tt = get_current_ms() - start;
        
        printf("iter=%d train_tt=%ld eval_tt=%ld ll=%f\n",
               iter, train_tt, eval_tt, ll);
    }
    
    void cgslda::train_one_doc(lda::document *p_doc) {
        
        int len = p_doc->_n;
        int did = p_doc->_did;
        int wid, tt, ttt;
        double u, psum;
        
        
        for (int i = 0; i < len; i ++) {
            wid = p_doc->_wids[i];
            tt  = p_doc->_topics[i];
            
            _nwk[wid][tt] --;
            _ndk[did][tt] --;
            _nk[tt] --;
            
            psum = 0.0;
            for (int k = 0; k < _K; k ++) {
                psum += (_ndk[did][k] + _alpha) * (_nwk[wid][k] + _beta) / (_nk[k] + _vbeta);
                _p[k] = psum;
            }
            
            u = unif01() * psum;
            ttt = -1;
            for (int k = 0; k < _K; k ++) {
                if (u < _p[k]) {
                    ttt = k;
                    break;
                }
            }
            
            _nwk[wid][ttt] ++;
            _ndk[did][ttt] ++;
            _nk[ttt] ++;
            p_doc->_topics[i] = ttt;
        }

    }
}

void cgslda_nips() {
    //   string path = "/Users/yulele/workspace/cpp_test/nips.train";
    string path = "nips.train";
    int V = 12420;
    int K = 1024;
    double alpha = 0.1;
    double beta  = 0.1;
    
    lda::corpus* cp = new lda::corpus(path);
    cp->load();
    
    int D = (int) cp->_docs.size();
    
    lda::cgslda* lda = new lda::cgslda(D, V, K, alpha, beta, cp);
    lda->init();
    
    for (int i = 0; i < 100; i ++)
        lda->train_one_iter(i);
    
    delete lda;
}

int main(int argc, const char * argv[]) {
        cgslda_nips();
    //    check_dv_mat();
//    warplda_nips();
    return 0;
}



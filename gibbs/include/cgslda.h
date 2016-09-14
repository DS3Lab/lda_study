//
//  cgslda.h
//  cpp_test
//
//  Created by yulele on 16/9/7.
//  Copyright (c) 2016年 yulele. All rights reserved.
//

#ifndef __cpp_test__cgslda__
#define __cpp_test__cgslda__

#include <stdio.h>
#include <string.h>
#include "csr.h"

namespace lda {

class cgslda {
public:
    int _D, _V, _K;
    double _alpha, _beta, _vbeta;
    
    corpus* _cp;
    int*  _nk;
    int** _nwk;
    int** _ndk;
    
    double* _p;
    
    cgslda(int D, int V, int K,
           double alpha, double beta,
           corpus* cp): _D(D), _V(V), _K(K), _alpha(alpha), _beta(beta) {
        
        _vbeta = V * _beta;
        _cp = cp;
        
        _nk = new int[_K];
        memset(_nk, 0, _K * (sizeof(int)));
        
        _nwk = new int*[_V];
        for (int w = 0; w < _V; w ++) {
            _nwk[w] = new int[_K];
            memset(_nwk[w], 0, _K * sizeof(int));
        }
        
        _ndk = new int*[_D];
        for (int d = 0; d < _D; d ++) {
            _ndk[d] = new int[_K];
            memset(_ndk[d], 0, _K * sizeof(int));
        }
        
        _p = new double[_K];
    }
    
    ~cgslda() {
        delete _cp;
        delete [] _nk;
        
        for (int w = 0; w < _V; w ++)
            delete [] _nwk[w];
        delete [] _nwk;
        
        for (int d = 0; d < _D; d ++)
            delete [] _ndk[d];
        delete [] _ndk;
    }
    
    void init();
    
    void train_one_iter(int iter);
    void train_one_doc(document* p_doc);
    
};
    
}

#endif /* defined(__cpp_test__cgslda__) */

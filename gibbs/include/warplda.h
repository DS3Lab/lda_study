//
//  warplda.h
//  cpp_test
//
//  Created by yulele on 16/9/7.
//  Copyright (c) 2016年 yulele. All rights reserved.
//

#ifndef __cpp_test__warplda__
#define __cpp_test__warplda__

#include <stdio.h>
#include <string.h>
#include "csr.h"
#include "xorshift.h"

namespace lda {

class warplda {
public:
    int _D, _V, _K;
    double _alpha, _beta, _vbeta;
    
    dv_mat* _mat;
    
    int* _nk;
    int* _tk;
    uint8_t* _used;
    int* _nidx;
    
    xorshift _gen;
    
    double _ll;
    
    warplda(int D, int V, int K,
            double alpha, double beta,
            dv_mat* mat): _D(D), _V(V), _K(K), _alpha(alpha), _beta(beta) {
        _vbeta = beta * V;
        _mat = mat;
        
        _nk = new int[_K];
        _tk = new int[_K];
        _used = new uint8_t[_K];
        _nidx = new int[_K];
        
        memset(_nk, 0, _K * sizeof(int));
    }
    
    ~warplda() {
        delete _mat;
        
        delete _nk;
        delete _tk;
    }
    
    void init();
    void visit_by_row();
    void visit_by_col();
    void build_nk();
    
    void train_one_iter(int iter);

    void check();
    
};
}

#endif /* defined(__cpp_test__warplda__) */

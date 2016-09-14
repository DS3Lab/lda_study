//
//  utils.cpp
//  cpp_test
//
//  Created by yulele on 16/9/6.
//  Copyright (c) 2016年 yulele. All rights reserved.
//

#include "utils.h"

using namespace std;

namespace lda {


void split(const string &s, char delim, vector<string> &elems) {
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}
    
    double llhw(int D, int V, int K, int** nwk, int** ndk, int *nk,
                double alpha, double beta, double vbeta, corpus* cp) {
        double lgamma_beta, lgamma_alpha;
        
        lgamma_alpha = lgamma(alpha);
        lgamma_beta  = lgamma(beta);
        
        double ll = 0.0;
        
        ll += K * lgamma(beta * V);
        
        for (int k = 0; k < K; k++) {
            ll -= lgamma(beta * V + nk[k]);
        }
        
        for (int w = 0; w < V; w++) {
            int* wa = nwk[w];
            for (int k = 0; k < K; k++) {
                if (wa[k] > 0) {
                    ll += lgamma(beta + wa[k]) - lgamma_beta;
                }
            }
        }
        
        for (int d = 0; d < D; d++) {
            document * p_doc = cp->_docs[d];
            int* da = ndk[d];
            ll += (lgamma(alpha * K) - lgamma(alpha * K + p_doc->_n));
            for (int k = 0; k < K; k ++) {
                if (da[k] > 0)
                    ll += lgamma(alpha + da[k]) - lgamma_alpha;
            }
        }
        
        return ll;
    }
    
}

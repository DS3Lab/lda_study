//
//  ftree.h
//  cpp_test
//
//  Created by yulele on 16/9/8.
//  Copyright (c) 2016年 yulele. All rights reserved.
//

#ifndef cpp_test_ftree_h
#define cpp_test_ftree_h

namespace lda {

class ftree {
public:
    double* _tree;
    int _len;
    
    ftree(int len): _len(len) {
        _tree = new double[2 * len];
    }
    
    ~ftree() {
        delete _tree;
    }
    
    void build(double* p) {
        for (int i = 2 * _len - 1; i > 0; i --) {
            if (i >= _len)
                _tree[i] = p[i - _len];
            else
                _tree[i] = _tree[i << 1] + _tree[(i << 1) + 1];
        }
    }
    
    void update(int idx, double val) {
        int i = idx + _len;
        double delta = val - _tree[i];
        while (i > 0) {
            _tree[i] += delta;
            i >>= 1;
        }
    }
    
    int sample(double u) {
        int i = 1;
        while (i < _len) {
            if (u < _tree[i << 1])
                i <<= 1;
            else {
                u = u - _tree[i << 1];
                i = (i << 1) + 1;
            }
        }
        return i - _len;
    }
    
    inline double get(int idx) {
        return _tree[idx + _len];
    }
    
    inline double first() {
        return _tree[1];
    }
};
}


#endif

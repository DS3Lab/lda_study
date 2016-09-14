//
//  CSRMat.h
//  cpp_test
//
//  Created by yulele on 16/9/6.
//  Copyright (c) 2016年 yulele. All rights reserved.
//

#ifndef __cpp_test__CSRMat__
#define __cpp_test__CSRMat__

#include <stdio.h>
#include <iostream>
#include <vector>

namespace lda {
    
const static int MH_STEPS = 1;


class document {
public:
    int _did;
    int _n;
    int* _wids;
    short* _topics;
    
    document(int n): _n(n) {
        _wids = new int[n];
    }
    
    ~document() {
//        printf("call ~document for did=%d\n", _did);
        delete []_wids;
    }
    
};

class corpus {
public:
    std::vector<document*> _docs;
    std::string _path;
    
    corpus(const std::string &path) {
        _docs = std::vector<document*>();
        _path = path;
    }
    
    ~corpus() {
        printf("call ~corpus\n");
        for (document* d: _docs) {
            delete d;
        }
    }
    
    void load();
    void show(int num);
    
    
};
    
class dv_mat {
public:
    int * _ws;
    
    struct entry {
        uint16_t tt;
        uint16_t mh[MH_STEPS];
    };
    
    entry* _data;
    
    int* _ds;
    int* _widx;
    
    int _D, _V, _N;
    
    dv_mat(const std::string path, int V): _V(V) {
        this->_D = 0;
        this->_N = 0;
        
        load(path);
    }
    
    ~dv_mat() {
        printf("call ~dv_mat\n");
        delete [] _ws;
        delete [] _ds;
        delete [] _widx;
        delete [] _data;
        
    }
    
    void alloc() {
        _ws = new int[_V + 1];
        _ds = new int[_D + 1];
        _widx = new int[_N];
        
        _data = new entry[_N];
    }
    
    void load(const std::string& path);
    void build_mat(std::vector<document*> &docs);
    void check(corpus &cp);
    void show();
};
    
class adj_list {
public:
    uint32_t * _ws;
    uint32_t * _dids;
    uint16_t * _topics;
    uint32_t * _dlens;
    
    int _N;
    int _V;
    int _D;
    
    adj_list(int V) {
        _N = 0;
        _V = V;
        _D = 0;
    }
    
    ~adj_list() {
        delete []_ws;
        delete []_dids;
        delete []_topics;
        delete []_dlens;
    }
    
    void load(const std::string& path);
    void build(std::vector<document *> &docs);
    
    void alloc () {
        _ws = new uint32_t[_V + 1];
        _dids = new uint32_t[_N];
        _topics = new uint16_t[_N];
        _dlens = new uint32_t[_D];
    }
    
    void check(std::vector<document *> &docs);
};
    
class warp_data {
public:
    int * _ws;
    short * _topics;
    short * _mhs;

    int * _ds;
    int * _widx;

    int _D, _V, _N;

    warp_data(const std::string& path, int V): _V(V) {
        _D = 0;
        _N = 0;
        load(path);
    }

    ~warp_data() {
        delete [] _ws;
        delete [] _topics;
        delete [] _mhs;
        delete [] _ds;
        delete [] _widx;
    }

    void load(const std::string& path) {
        corpus cp = corpus(path);
        cp.load();
        
        _N = 0;
        _D = (int) cp._docs.size();
        int V = -1;
        for (int d = 0; d < _D; d ++) {
            _N += cp._docs[d]->_n;
            
            for (int w = 0; w < cp._docs[d]->_n; w ++)
                if (cp._docs[d]->_wids[w] > V)
                    V = cp._docs[d]->_wids[w];
        }
        
        printf("V=%d\n", V);
        
        alloc();
        
        build_mat(cp._docs);

    }
    
    void build_mat(std::vector<document *> &docs) {
        printf("D=%d V=%d N=%d\n", _D, _V, _N);
        int* wcnt = new int[_V];
        document* p_doc;
        
        
        // count word and build doc start idx
        _ds[0] = 0;
        for (int d = 0; d < _D; d ++) {
            p_doc = docs[d];
            for (int w = 0; w < p_doc->_n; w ++)
                wcnt[p_doc->_wids[w]] ++;
            _ds[d + 1] = _ds[d] + p_doc->_n;
        }
        
        int total = 0;
        for (int w = 0; w < _V; w ++)
            total += wcnt[w];
        
        // build word start idx
        _ws[0] = 0;
        for (int w = 0; w < _V; w ++)
            _ws[w + 1] = _ws[w] + wcnt[w];
        
        
        // build doc to word reverse idx
        for (int d = _D - 1; d >= 0; d --) {
            int di = _ds[d];
            p_doc = docs[d];
            int wid;
            for (int w = 0; w < p_doc->_n; w ++) {
                wid = p_doc->_wids[w];
                wcnt[wid] --;
                int pos = _ws[wid] + wcnt[wid];
                _widx[di + w] = pos;
            }
        }
        
        delete[] wcnt;
    }

    void alloc() {
        _ws = new int[_V + 1];
        _ds = new int[_D + 1];
        _widx = new int[_N];

        _topics = new short[_N];
        _mhs = new short[_N * MH_STEPS];
    }
};
    
}
#endif /* defined(__cpp_test__CSRMat__) */

//
//  CSRMat.cpp
//  cpp_test
//
//  Created by yulele on 16/9/6.
//  Copyright (c) 2016年 yulele. All rights reserved.
//

#include "csr.h"
#include "utils.h"

using namespace std;

namespace lda {

    void corpus::load() {
        
        FILE* file = fopen(_path.c_str(), "r");
        char buf[102400];
        
        int did = 0;
        if (file == nullptr) {
            printf("Cannot open file %s\n", _path.c_str());
            exit(1);
        } else {
            while (fgets(buf, 102400, file) != nullptr) {
                string str = string(buf);
                vector<string> splits = vector<string>();
                split(str, ' ', splits);
                
                int len = (int) splits.size();
                document* doc = new document(len - 1);
                for (int i = 1; i < len; i ++)
                    doc->_wids[i - 1] = stoi(splits[i].c_str());
                
                doc->_did = did;
                did ++;
                
                _docs.push_back(doc);
            }
            
            printf("did=%d\n", did);
        }
    }
    
    void corpus::show(int num) {
            document* doc = _docs[num];
            printf("doc len is %d\n", doc->_n);
            for (int j = 0; j < doc->_n; j ++)
                printf("%d,", doc->_wids[j]);
            printf("\n");
    }
    
    void dv_mat::load(const std::string &path) {
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
        
//        check(cp);
    }
    
    void dv_mat::build_mat(vector<document*> &docs) {
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
    
    void dv_mat::check(lda::corpus &cp) {
        
        printf("D=%d V=%d N=%d\n", _D, _V, _N);
        
        int* wr_idx = new int[_N];
        for (int w = 0; w < _V; w ++) {
            for (int wi = _ws[w]; wi < _ws[w + 1]; wi ++)
                wr_idx[wi] = w;
        }
        
        for (int d = 0; d < _D; d ++) {
            document* p_doc = cp._docs[d];
            int di = _ds[d];
            for (int w = 0; w < p_doc->_n; w ++) {
                int pos = _widx[di + w];
                if (p_doc->_wids[w] != wr_idx[pos])
                    printf("Error found\n");
            }
        }
        
        delete [] wr_idx;
        
        printf("done\n");
        
    }
    
    void dv_mat::show() {
        
    }
    
    void adj_list::load(const std::string &path) {
        corpus cp = corpus(path);
        cp.load();
        
        _N = 0;
        _D = (int) cp._docs.size();
        for (int d = 0; d < _D; d ++) {
            _N += cp._docs[d]->_n;
        }
        
        alloc();
        build(cp._docs);
//        check(cp._docs);
    }
    
    void adj_list::build(std::vector<document *> &docs) {
        int * wcnt = new int[_V];
        document * p_doc;
        
        for (int d = 0; d < _D; d ++) {
            p_doc = docs[d];
            _dlens[d] = p_doc->_n;
            for (int w = 0; w < p_doc->_n; w ++)
                wcnt[p_doc->_wids[w]] ++;
        }
    
        _ws[0] = 0;
        for (int w = 0; w < _V; w ++)
            _ws[w + 1] = _ws[w] + wcnt[w];
        
        
        // build did array
        for (int d = _D - 1; d >= 0; d --) {
            p_doc = docs[d];
            
            int wid;
            for (int w = 0; w < p_doc->_n; w ++) {
                wid = p_doc->_wids[w];
                wcnt[wid] --;
                int pos = _ws[wid] + wcnt[wid];
                _dids[pos] = d;
            }
        }
        
        delete wcnt;
    }
    
    void adj_list::check(std::vector<document *> &docs) {
        printf("check start\n");
        int * d_idx = new int[_D];
        memset(d_idx, 0, _D * sizeof(int));
        
        document * p_doc;
        for (int w = 0; w < _V; w ++) {
            for (int wi = _ws[w]; wi < _ws[w + 1]; wi ++) {
                int did = _dids[wi];
                p_doc = docs[did];
                if (w != p_doc->_wids[d_idx[did] ++]) {
                    printf("Error\n");
                }
            }
        }
        
        printf("done!\n");
        
        delete d_idx;
    }
}



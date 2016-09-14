//#include "csr.h"
#include "lda_jni_Utils.h"
#include "csr.h"
#include <omp.h>
#include <random>
#include <math.h>
#include <string.h>
#include <sys/time.h>

class xorshift
{
	uint64_t s[16];
	int p;
	uint64_t x; /* The state must be seeded with a nonzero value. */

	uint64_t xorshift1024star(void) {
		uint64_t s0 = s[ p ];
		uint64_t s1 = s[ p = (p+1) & 15 ];
		s1 ^= s1 << 31; // a
		s1 ^= s1 >> 11; // b
		s0 ^= s0 >> 30; // c
		return ( s[p] = s0 ^ s1 ) * UINT64_C(1181783497276652981);
	}
	uint64_t xorshift128plus(void) {
		uint64_t x = s[0];
		uint64_t const y = s[1];
		s[0] = y;
		x ^= x << 23; // a
		x ^= x >> 17; // b
		x ^= y ^ (y >> 26); // c
		s[1] = x;
		return x + y;
	}
	uint64_t xorshift64star(void) {
		x ^= x >> 12; // a
		x ^= x << 25; // b
		x ^= x >> 27; // c
		return x * UINT64_C(2685821657736338717);
	}
	public:

	using result_type=uint64_t;

	xorshift() : p(0), x((uint64_t)std::rand() * RAND_MAX + std::rand()){
		for (int i = 0; i < 16; i++)
		{
			s[i] = xorshift64star();
		}
	}
	uint64_t operator()(){
		return xorshift128plus();
		//return xorshift64star();
	}
	inline uint32_t Rand32(){
		return (uint32_t)xorshift128plus();
	}
    
    inline uint32_t rand_int(int max) {
        return Rand32() % max;
    }
    
    inline double rand_double() {
        return Rand32() * 1.0 / std::numeric_limits<uint32_t>::max();
    }
    
    
	void MakeBuffer(void *p, size_t len)
	{
		int N = (int) len / sizeof(uint32_t);
		uint32_t *arr = (uint32_t *)p;
		for (int i = 0; i < N; i++)
			arr[i] = (uint32_t)(*this)();
		int M = len % sizeof(uint32_t);
		if (M > 0)
		{
			uint32_t k = (uint32_t)(*this)();
			memcpy(arr + N, &k, M);
		}
	}
	uint64_t max() {return std::numeric_limits<uint64_t>::max();}
	uint64_t min() {return std::numeric_limits<uint64_t>::min();}
};

inline static long get_current_ms() {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    
    long timestamp = (long) tv.tv_sec * 1000 + tv.tv_usec / 1000;
    return timestamp;
}

inline static std::mt19937 & urng() {
    static std::mt19937 u{};
    return u;
}

inline static double unif01() {
    static std::uniform_real_distribution<double> d(0.0, 1.0);
    return d(urng());
}

inline static int rand_int(int from, int to) {
    static std::uniform_int_distribution<> d{};
    using param_t = std::uniform_int_distribution<>::param_type;
    return d(urng(), param_t{from, to});
}


void init_warp(int D, int V, int K,
        int * ws, int * ds, int * widx, int * nk,
        short * topics, short * mhs, int MH_STEPS,
        xorshift&  gen) {
    
    printf("D=%d, V=%d, K=%d\n", D, V, K);
    short tt;
    for (int w = 0; w < V; w ++) {
        for (int wi = ws[w]; wi < ws[w + 1]; wi ++) {
            tt = (short)gen.rand_int(K);
            topics[wi] = tt;
            nk[tt] ++;
            for (int m = 0; m < MH_STEPS; m ++)
                mhs[wi * MH_STEPS + m] = tt;
        }
    }
}

double visit_by_row(int D, int V, int K,
        int * ds, int * widx, int * nk,
        int * tk, int * nidx, int * used,
        short * topics, short * mhs, int MH_STEPS,
        double alpha, double beta, double vbeta,
        xorshift& gen) {
    short tt, ttt, k;
    int N, idx;
    double p;

    double lgamma_alpha = lgamma(alpha);
    double ll = 0;
    
    // traverse each doc in sequential order
    for (int d = 0; d < D; d ++) {
        N = ds[d + 1] - ds[d];
        if (N == 0)
            continue;

        // compute dk for doc on the fly
        memset(tk, 0, K * sizeof(int));
        //memset(used, 0, K * sizeof(int));
        
        idx = 0;
        for (int di = ds[d]; di < ds[d + 1]; di ++) {
            k = topics[widx[di]];
            tk[k] ++;
            if (used[k] == 0) {
                used[k] = 1;
                nidx[idx ++] = k;
            }
        }

        // compute llhw
        for (int i = 0; i < idx; i ++) {
            ll += lgamma(alpha + tk[nidx[i]]) - lgamma_alpha;
            used[nidx[i]] = 0;
        }


        ll -= lgamma(alpha * K + N) - lgamma(alpha * K);

        // traverse doc first time, accept word proposal
        for (int di = ds[d]; di < ds[d + 1]; di ++) {
            tt = topics[widx[di]];
            nk[tt] --;
            tk[tt] --;

            for (int m = 0; m < MH_STEPS; m ++) {
                ttt = mhs[widx[di] * MH_STEPS + m];

                p = ((tk[ttt] + alpha) * (nk[tt] + vbeta))
                        / ((tk[tt] + alpha) * (nk[ttt] + vbeta));

                if (gen.rand_double() < p)
                    tt = ttt;
            }

            nk[tt] ++;
            tk[tt] ++;
            topics[widx[di]] = tt;
        }

        // traverse doc second time, generate doc proposal
        p = (alpha * K) / (alpha * K + N);
        for (int di = ds[d]; di < ds[d + 1]; di ++) {
            for (int m = 0; m < MH_STEPS; m ++) {
                int si = widx[di] * MH_STEPS;
                if (gen.rand_double() < p)
                    mhs[si + m] = gen.rand_int(K);
                else
                    mhs[si + m] = topics[widx[ds[d] + gen.rand_int(N)]];
            }
        }
    }
    return ll;
}

double visit_by_col(int D, int V, int K,
        int * ws, int * nk,
        int * tk, int * nidx, int * used,
        short * topics, short * mhs, int MH_STEPS,
        double alpha, double beta, double vbeta,
        xorshift&  gen) {

    short tt, ttt, k;
    int N, idx;
    double p;

    double lgamma_beta = lgamma(beta);
    double ll = 0.0;

    for (int w = 0; w < V; w ++) {
        N = ws[w + 1] - ws[w];

        if (N == 0)
            continue;

        // compute wk for word on the fly
        memset(tk, 0, K * sizeof(int));
        //memset(used, 0, K * sizeof(int));
        
        idx = 0;
        for (int wi = ws[w]; wi < ws[w + 1]; wi ++) {
            k = topics[wi];
            tk[k] ++;
            if (used[k] == 0) {
                nidx[idx ++] = k;
                used[k] = 1;
            }
        }
            
        // compute llhw
        for (int i = 0; i < idx; i ++) {
            ll += lgamma(beta + tk[nidx[i]]) - lgamma_beta;
            used[nidx[i]] = 0; 
        }

        // traverse w first time, accept doc proposal
        for (int wi = ws[w]; wi < ws[w + 1]; wi ++) {
            tt = topics[wi];
            tk[tt] --;
            nk[tt] --;
            
            for (int m = 0; m < MH_STEPS; m ++) {
                ttt = mhs[wi * MH_STEPS + m];
                
                double b = tk[tt]  + beta;
                double d = nk[tt]  + vbeta;
                double a = tk[ttt] + beta;
                double c = nk[ttt] + vbeta;
                
                if (gen.rand_double() < (a * d) / (b * c)) {
                    tt = ttt;
                }
            }
            
            tk[tt] ++;
            nk[tt] ++;
            topics[wi] = tt;
        }

        // traverse w second time, compute word proposal
        p = (K * beta) / (K * beta + N);
        for (int wi = ws[w]; wi < ws[w + 1]; wi ++) {
            for (int m = 0; m < MH_STEPS; m ++) {
                if (gen.rand_double() < p) {
                    mhs[wi * MH_STEPS + m] = gen.rand_int(K);
                } else {
                    mhs[wi * MH_STEPS + m] = topics[ws[w] + gen.rand_int(N)];
                }
            }
        }
    }

    for (int k = 0; k < K; k ++) {
        if (nk[k] > 0)
          ll += - lgamma(nk[k] + vbeta) + lgamma(vbeta);
    }

    return ll;

}

double train_for_one_iter(int iter, int D, int V, int K, int N,
        int * ws, int * ds, int * widx, int * nk,
        int * tk, int * nidx, int * used,
        short * topics, short * mhs, int MH_STEPS,
        double alpha, double beta, double vbeta) {

    printf("train one iter with params D=%d V=%d K=%d N=%d alpha=%f beta=%f\n",
            D, V, K, N, alpha, beta);
    
    xorshift gen;

    double ll = 0;
    long start, end, stt;
    long row_tt, col_tt, iter_tt;

    start = get_current_ms();
    stt = start;
    ll += visit_by_col(D, V, K, ws, nk, tk, nidx, used,
            topics, mhs, MH_STEPS, alpha, beta, vbeta, gen);
    end = get_current_ms();
    col_tt = end - start;

    start = get_current_ms();
    ll += visit_by_row(D, V, K, ds, widx, nk, tk, nidx, used,
            topics, mhs, MH_STEPS, alpha, beta, vbeta, gen);
    end = get_current_ms();
    row_tt = end - start;
    iter_tt = end - stt;

    printf("iter=%d iter_tt=%ld col_tt=%ld row_tt=%ld ll=%f %f M tokens/s\n", 
            iter, iter_tt, col_tt, row_tt, ll, ((double) N/1e6) / (iter_tt / 1e3));
    return ll;
}

const static int THREAD_NUM = 4;

struct thread_params {
    xorshift gen;
    int * tk;
    int * nidx;
    int * used;
    int * nk;
    int * nk_new;
    double ll;
    int K;

    thread_params() {}
    
    void init(int K) {
        tk = new int[K];
        nidx = new int[K];
        used = new int[K];
        nk = new int[K];
        nk_new = new int[K];
        ll = 0;
        this->K = K;
        memset(used, 0, K * sizeof(int));
        memset(tk, 0, K * sizeof(int));
    }

    ~thread_params() {
        delete [] tk;
        delete [] nidx;
        delete [] used;
        delete [] nk;
        delete [] nk_new;
    }

};

void init_warp_parallel(int D, int V, int K,
        int * ws, int * ds, int * widx, int * nk,
        short * topics, short * mhs, int MH_STEPS,
        thread_params * params) {

#pragma omp parallel for num_threads(THREAD_NUM) 
    for (int w = 0; w < V; w ++) {
        int tn = omp_get_thread_num();
        for (int wi = ws[w]; wi < ws[w + 1]; wi ++) {
            short tt = (short)params[tn].gen.rand_int(K);
            topics[wi] = tt;
            //nk[tt] ++;
            __sync_fetch_and_add(&nk[tt], 1);
            for (int m = 0; m < MH_STEPS; m ++)
                mhs[wi * MH_STEPS + m] = tt;
        }
    }
}

void parallel_reduce_nk(int * global_nk, thread_params * params, int K) {
#pragma omp parallel for num_threads(THREAD_NUM)
    for (int k = 0; k < K; k++) {
        global_nk[k] = 0;
        for (int i = 0; i < THREAD_NUM; i ++)
            global_nk[k] += params[i].nk_new[k];
    }
}

void parallel_copy_nk(int * global_nk, thread_params * params, int K) {
#pragma omp parallel for num_threads(THREAD_NUM) 
    for (int k = 0; k < K; k ++) {
        for (int i = 0; i < THREAD_NUM; i ++)
            params[i].nk[k] = global_nk[k];
    }
    
    for (int i = 0; i < THREAD_NUM; i ++)
        memset(params[i].nk_new, 0, sizeof(int) * K);
}

double visit_by_row_parallel(int D, int V, int K,
        int * ds, int * widx, int * global_nk,
        short * topics, short * mhs, int MH_STEPS,
        float alpha, float beta, float vbeta,
        thread_params * params) {

    float lgamma_alpha = lgamma(alpha);
    float global_ll = 0;

    parallel_copy_nk(global_nk, params, K);
    
    // traverse each doc
#pragma omp parallel for num_threads(THREAD_NUM)
    for (int d = 0; d < D; d ++) {
        int tn = omp_get_thread_num();
        int N = ds[d + 1] - ds[d];
        if (N == 0)
            continue;
        short tt, ttt, k;
        int  idx;
        float p;


        int * tk = params[tn].tk;
        int * nidx = params[tn].nidx;
        int * used = params[tn].used;
        int * nk = params[tn].nk;
        int * nk_new = params[tn].nk_new;

        // compute dk for doc on the fly
        memset(tk, 0, K * sizeof(int));
        //memset(used, 0, K * sizeof(int));
        
        idx = 0;
        for (int di = ds[d]; di < ds[d + 1]; di ++) {
            k = topics[widx[di]];
            tk[k] ++;
            if (used[k] == 0) {
                used[k] = 1;
                nidx[idx ++] = k;
            }
        }

        // compute llhw
        for (int i = 0; i < idx; i ++) {
            params[tn].ll += lgamma(alpha + tk[nidx[i]]) - lgamma_alpha;
            used[nidx[i]] = 0;
        }


        params[tn].ll -= lgamma(alpha * K + N) - lgamma(alpha * K);

        // traverse doc first time, accept word proposal
        for (int di = ds[d]; di < ds[d + 1]; di ++) {
            tt = topics[widx[di]];
            nk[tt] --;
            tk[tt] --;

            for (int m = 0; m < MH_STEPS; m ++) {
                ttt = mhs[widx[di] * MH_STEPS + m];

                p = ((tk[ttt] + alpha) * (nk[tt] + vbeta))
                        / ((tk[tt] + alpha) * (nk[ttt] + vbeta));

                if (params[tn].gen.rand_double() < p)
                    tt = ttt;
            }

            nk[tt] ++;
            tk[tt] ++;
            topics[widx[di]] = tt;
            nk_new[tt] ++;
        }

        // traverse doc second time, generate doc proposal
        p = (alpha * K) / (alpha * K + N);
        for (int di = ds[d]; di < ds[d + 1]; di ++) {
            for (int m = 0; m < MH_STEPS; m ++) {
                int si = widx[di] * MH_STEPS;
                if (params[tn].gen.rand_double() < p)
                    mhs[si + m] = params[tn].gen.rand_int(K);
                else
                    mhs[si + m] = topics[widx[ds[d] + params[tn].gen.rand_int(N)]];
            }
        }
    }

    parallel_reduce_nk(global_nk, params, K);

    for (int i = 0; i < THREAD_NUM; i ++)
        global_ll += params[i].ll;
    return global_ll;
}

double visit_by_col_parallel(int D, int V, int K,
        int * ws, int * global_nk,
        short * topics, short * mhs, int MH_STEPS,
        float alpha, float beta, float vbeta,
        thread_params * params) {

    double lgamma_beta = lgamma(beta);
    float global_ll = 0.0;

    parallel_copy_nk(global_nk, params, K);

#pragma omp parallel for num_threads(THREAD_NUM)
    for (int w = 0; w < V; w ++) {
        int N = ws[w + 1] - ws[w];
        short tt, ttt, k;
        int  idx;
        float p;

        if (N == 0)
            continue;

        int tn = omp_get_thread_num();
        int * tk = params[tn].tk;
        int * nidx = params[tn].nidx;
        int * used = params[tn].used;
        int * nk = params[tn].nk;
        int * nk_new = params[tn].nk_new;

        // compute wk for word on the fly
        memset(tk, 0, K * sizeof(int));
        //memset(used, 0, K * sizeof(int));
        
        idx = 0;
        for (int wi = ws[w]; wi < ws[w + 1]; wi ++) {
            k = topics[wi];
            tk[k] ++;
            if (used[k] == 0) {
                nidx[idx ++] = k;
                used[k] = 1;
            }
        }
            
        // compute llhw
        for (int i = 0; i < idx; i ++) {
            params[tn].ll += lgamma(beta + tk[nidx[i]]) - lgamma_beta;
            used[nidx[i]] = 0; 
        }

        // traverse w first time, accept doc proposal
        for (int wi = ws[w]; wi < ws[w + 1]; wi ++) {
            tt = topics[wi];
            tk[tt] --;
            nk[tt] --;
            
            for (int m = 0; m < MH_STEPS; m ++) {
                ttt = mhs[wi * MH_STEPS + m];
                
                float b = tk[tt]  + beta;
                float d = nk[tt]  + vbeta;
                float a = tk[ttt] + beta;
                float c = nk[ttt] + vbeta;
                
                if (params[tn].gen.rand_double() < (a * d) / (b * c)) {
                    tt = ttt;
                }
            }
            
            tk[tt] ++;
            nk[tt] ++;
            nk_new[tt] ++;
            topics[wi] = tt;
        }

        // traverse w second time, compute word proposal
        p = (K * beta) / (K * beta + N);
        for (int wi = ws[w]; wi < ws[w + 1]; wi ++) {
            for (int m = 0; m < MH_STEPS; m ++) {
                if (params[tn].gen.rand_double() < p) {
                    mhs[wi * MH_STEPS + m] = params[tn].gen.rand_int(K);
                } else {
                    mhs[wi * MH_STEPS + m] = topics[ws[w] + params[tn].gen.rand_int(N)];
                }
            }
        }
    }

    parallel_reduce_nk(global_nk, params, K);

    for (int k = 0; k < K; k ++) {
        if (global_nk[k] > 0)
          global_ll += - lgamma(global_nk[k] + vbeta) + lgamma(vbeta);
    }
    
    for (int i = 0; i < THREAD_NUM; i ++)
        global_ll += params[i].ll;

    return global_ll;
}


void train_for_one_iter_parallel(int iter, int D, int V, int K, int N,
        int * ws, int * ds, int * widx, int * nk,
        short * topics, short * mhs, int MH_STEPS,
        float alpha, float beta, float vbeta,
        thread_params * params) {
    float ll = 0;

    long start, end, stt;
    long row_tt, col_tt, iter_tt;

    for (int i = 0; i < THREAD_NUM; i ++)
        params[i].ll = 0;

    start = get_current_ms();
    stt = start;
    ll += visit_by_col_parallel(D, V, K, ws, nk,
            topics, mhs, MH_STEPS, alpha, beta, vbeta, params);
    end = get_current_ms();
    col_tt = end - start;

    for (int i = 0; i < THREAD_NUM; i ++)
        params[i].ll = 0;

    start = get_current_ms();
    ll += visit_by_row_parallel(D, V, K, ds, widx, nk,
            topics, mhs, MH_STEPS, alpha, beta, vbeta, params);
    end = get_current_ms();
    row_tt = end - start;
    iter_tt = end - stt;

    printf("iter=%d iter_tt=%ld col_tt=%ld row_tt=%ld ll=%f %f M tokens/s\n", 
            iter, iter_tt, col_tt, row_tt, ll, ((double) N/1e6) / (iter_tt / 1e3));

}

void train_warp_nips_parallel() {
    std::string path = "data/nips.train";
    int V = 12420;
    int K = 1024;
    float alpha = 0.1;
    float beta  = 0.1;
    float vbeta = V * beta;

    lda::warp_data * data = new lda::warp_data(path, V);

    printf("load finished\n");
    int * nk = new int[K];
    memset(nk, 0, K * sizeof(int));

    thread_params params[THREAD_NUM];
    for (int i = 0; i < THREAD_NUM; i ++)
        params[i].init(K);

    
    init_warp_parallel(data->_D, data->_V, K, 
            data->_ws, data->_ds, data->_widx, nk,
            data->_topics, data->_mhs, lda::MH_STEPS, params);

    printf("init finished\n");

    for (int i = 0; i < 200; i ++) {
        train_for_one_iter_parallel(i, data->_D, data->_V, K, data->_N,
                data->_ws, data->_ds, data->_widx, nk,
                data->_topics, data->_mhs, lda::MH_STEPS,
                alpha, beta, vbeta, params);
    }

    delete data;
    delete [] nk;

}

void train_warp_nips() {
    std::string path = "data/nips.train";
    int V = 12420;
    int K = 1024;
    float alpha = 0.1;
    float beta  = 0.1;
    float vbeta = V * beta;

    lda::warp_data * data = new lda::warp_data(path, V);

    printf("load finished\n");
    int * nk = new int[K];
    int * tk = new int[K];
    int * nidx = new int[K];
    int * used = new int[K];

    memset(nk, 0, K * sizeof(int));
    memset(used, 0, K * sizeof(int));
    
    xorshift gen;
    init_warp(data->_D, data->_V, K, 
            data->_ws, data->_ds, data->_widx, nk,
            data->_topics, data->_mhs, lda::MH_STEPS, gen);

    printf("init finished\n");

    for (int i = 0; i < 100; i ++) {
        train_for_one_iter(i, data->_D, data->_V, K, data->_N,
                data->_ws, data->_ds, data->_widx, nk,
                tk, nidx, used,
                data->_topics, data->_mhs, lda::MH_STEPS,
                alpha, beta, vbeta);
    }

    delete data;
    delete [] nk;
}

JNIEXPORT jdouble JNICALL Java_lda_jni_Utils_warpOneIter
  (JNIEnv * env, jclass obj, 
   jint jiter, jint jD, jint jV, jint jK, jint jN, 
   jintArray jws, jintArray jds, jintArray jwidx, jintArray jnk, 
   jintArray jtk, jintArray jnidx, jintArray jused, 
   jshortArray jtopics, jshortArray jmhs, 
   jint jMH_STEPS, 
   jdouble jalpha, jdouble jbeta, jdouble jvbeta) {
    jboolean is_copy;
    
    int * ws = (int*)env->GetPrimitiveArrayCritical(jws, &is_copy);
    int * ds = (int*)env->GetPrimitiveArrayCritical(jds, &is_copy);
    int * widx = (int*)env->GetPrimitiveArrayCritical(jwidx, &is_copy);
    int * nk = (int*)env->GetPrimitiveArrayCritical(jnk, &is_copy);
    int * tk = (int*)env->GetPrimitiveArrayCritical(jtk, &is_copy);
    int * nidx = (int*)env->GetPrimitiveArrayCritical(jnidx, &is_copy);
    int * used = (int*)env->GetPrimitiveArrayCritical(jused, &is_copy);
    short * topics = (short*)env->GetPrimitiveArrayCritical(jtopics, &is_copy);
    short * mhs = (short*)env->GetPrimitiveArrayCritical(jmhs, &is_copy);
    
    /*
    jsize wslen = env->GetArrayLength(jws);
    jsize dslen = env->GetArrayLength(jds);
    jsize widxlen = env->GetArrayLength(jwidx);
    jsize nklen = env->GetArrayLength(jnk);
    jsize tklen = env->GetArrayLength(jtk);
    jsize nidxlen = env->GetArrayLength(jnidx);
    jsize usedlen = env->GetArrayLength(jused);
    jsize topicslen = env->GetArrayLength(jtopics);
    jsize mhslen = env->GetArrayLength(jmhs);
    printf("wslen=%d dslen=%d widxlen=%d nklen=%d tklen=%d nidxlen=%d usedlen=%d topicslen=%d mhslen=%d\n", 
            wslen, dslen, widxlen, nklen, tklen, nidxlen,
            usedlen, topicslen, mhslen);
    */

    
    double ll = train_for_one_iter(jiter, jD, jV, jK, jN,
            ws, ds, widx, nk,
            tk, nidx, used,
            topics, mhs,
            jMH_STEPS,
            jalpha, jbeta, jvbeta);

    env->ReleasePrimitiveArrayCritical(jws, ws, 0);
    env->ReleasePrimitiveArrayCritical(jds, ds, 0);
    env->ReleasePrimitiveArrayCritical(jwidx, widx, 0);
    env->ReleasePrimitiveArrayCritical(jnk, nk, 0);
    env->ReleasePrimitiveArrayCritical(jtk, tk, 0);
    env->ReleasePrimitiveArrayCritical(jnidx, nidx, 0);
    env->ReleasePrimitiveArrayCritical(jused, used, 0);
    env->ReleasePrimitiveArrayCritical(jtopics, topics, 0);
    env->ReleasePrimitiveArrayCritical(jmhs, mhs, 0);

    return ll;
        
  }

JNIEXPORT jdouble JNICALL Java_lda_jni_Utils_warpInit
  (JNIEnv * env, jclass obj, 
   jint jD, jint jV, jint jK, 
   jintArray jws, jintArray jds, jintArray jwidx, jintArray jnk,
   jshortArray jtopics, jshortArray jmhs, 
   jint jMH_STEPS) {
    jboolean is_copy;
    int * ws = (int*) env->GetPrimitiveArrayCritical(jws, &is_copy);
    int * ds = (int*) env->GetPrimitiveArrayCritical(jds, &is_copy);
    int * widx = (int*) env->GetPrimitiveArrayCritical(jwidx, &is_copy);
    int * nk = (int*) env->GetPrimitiveArrayCritical(jnk, &is_copy);
    short * topics = (short*) env->GetPrimitiveArrayCritical(jtopics, &is_copy);
    short * mhs = (short *) env->GetPrimitiveArrayCritical(jmhs, &is_copy);
    int D = (int) jD;
    int V = (int) jV;
    int K = (int) jK;
    int MH_STEPS = (int) jMH_STEPS;
    
    /*
    jsize wslen = env->GetArrayLength(jws);
    jsize dslen = env->GetArrayLength(jds);
    jsize widxlen = env->GetArrayLength(jwidx);
    jsize nklen = env->GetArrayLength(jnk);
    jsize topicslen = env->GetArrayLength(jtopics);
    jsize mhslen = env->GetArrayLength(jmhs);
    printf("wslen=%d dslen=%d widxlen=%d nklen=%d topicslen=%d mhslen=%d\n", 
            wslen, dslen, widxlen, nklen, topicslen, mhslen);
    */

    xorshift gen;
    init_warp(D, V, K, ws, ds, widx, nk, topics, mhs, MH_STEPS, gen);
    printf("init finished\n");

    env->ReleasePrimitiveArrayCritical(jws, ws, 0);
    env->ReleasePrimitiveArrayCritical(jds, ds, 0);
    env->ReleasePrimitiveArrayCritical(jwidx, widx, 0);
    env->ReleasePrimitiveArrayCritical(jnk, nk, 0);
    env->ReleasePrimitiveArrayCritical(jtopics, topics, 0);
    env->ReleasePrimitiveArrayCritical(jmhs, mhs, 0);
    return 0.0;
  }

int main(int argc, char * argv[]) {

    //train_warp_nips();
    train_warp_nips_parallel();
}


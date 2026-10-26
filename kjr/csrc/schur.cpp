// ODM-KJR bundle adjustment: reduced camera system by point elimination (Schur complement),
// accumulated point by point (the way Ceres does it) instead of with sparse matrix products.
//
// Observations sorted by point; obs o of point p has K global columns cols[o, :] (-1 = unused),
// weighted global Jacobian G[o] (2 x K), weighted point Jacobian P[o] (2 x 3), weighted residual r[o] (2).
//   V_p   = sum_o P_o^T P_o + 1e-6 I,      bp_p = sum_o P_o^T r_o
//   S    += sum_o G_o^T G_o - sum_{o,o'} G_o^T (P_o V_p^-1 P_o'^T) G_o'
//   b    += sum_o G_o^T r_o - sum_o G_o^T P_o V_p^-1 bp_p
// S (NG x NG, row major) and b (NG) are ADDED to. Vinv (npt x 3 x 3) and bp (npt x 3) are written
// for the back-substitution. Thread-private S, reduced at the end.
// Build: g++ -O3 -march=native -fopenmp -shared -fPIC schur.cpp -o ../_schur.so
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <omp.h>

static inline void inv3(const double* a, double* o) {
    double c00 = a[4] * a[8] - a[5] * a[7], c01 = a[5] * a[6] - a[3] * a[8], c02 = a[3] * a[7] - a[4] * a[6];
    double det = a[0] * c00 + a[1] * c01 + a[2] * c02, id = 1.0 / det;
    o[0] = c00 * id; o[1] = (a[2] * a[7] - a[1] * a[8]) * id; o[2] = (a[1] * a[5] - a[2] * a[4]) * id;
    o[3] = c01 * id; o[4] = (a[0] * a[8] - a[2] * a[6]) * id; o[5] = (a[2] * a[3] - a[0] * a[5]) * id;
    o[6] = c02 * id; o[7] = (a[1] * a[6] - a[0] * a[7]) * id; o[8] = (a[0] * a[4] - a[1] * a[3]) * id;
}

extern "C" int schur_points(int64_t npt, const int64_t* ptr, int K, const int64_t* cols,
                            const double* G, const double* P, const double* r, int64_t NG,
                            double* S, double* b, double* Vinv, double* bpo, int nthreads) {
    if (nthreads < 1) nthreads = 1;
    std::vector<double*> St(nthreads, nullptr), bt(nthreads, nullptr);
    for (int t = 0; t < nthreads; ++t) {
        St[t] = t == 0 ? S : (double*)calloc((size_t)NG * NG, sizeof(double));
        bt[t] = t == 0 ? b : (double*)calloc((size_t)NG, sizeof(double));
        if (!St[t] || !bt[t]) return -1;
    }
#pragma omp parallel num_threads(nthreads)
    {
        int tid = omp_get_thread_num();
        double* Sl = St[tid]; double* bl = bt[tid];
        std::vector<double> GP, Q;           // GP[o] = G_o^T P_o (K x 3), Q = GP V^-1
#pragma omp for schedule(dynamic, 2048)
        for (int64_t p = 0; p < npt; ++p) {
            int64_t o0 = ptr[p], o1 = ptr[p + 1], m = o1 - o0;
            double V[9] = {1e-6, 0, 0, 0, 1e-6, 0, 0, 0, 1e-6}, bp[3] = {0, 0, 0}, Vi[9];
            for (int64_t o = o0; o < o1; ++o) {
                const double* Po = P + o * 6; const double* ro = r + o * 2;
                for (int i = 0; i < 3; ++i) {
                    bp[i] += Po[i] * ro[0] + Po[3 + i] * ro[1];
                    for (int j = 0; j < 3; ++j) V[i * 3 + j] += Po[i] * Po[j] + Po[3 + i] * Po[3 + j];
                }
            }
            inv3(V, Vi);
            memcpy(Vinv + p * 9, Vi, sizeof Vi); memcpy(bpo + p * 3, bp, sizeof bp);
            GP.assign((size_t)m * K * 3, 0.0); Q.assign((size_t)m * K * 3, 0.0);
            for (int64_t k = 0; k < m; ++k) {
                int64_t o = o0 + k;
                const double* Go = G + o * 2 * K; const double* Po = P + o * 6;
                const int64_t* co = cols + o * K;
                double* gp = &GP[k * K * 3]; double* q = &Q[k * K * 3];
                for (int i = 0; i < K; ++i) {
                    double g0 = Go[i], g1 = Go[K + i];
                    for (int j = 0; j < 3; ++j) gp[i * 3 + j] = g0 * Po[j] + g1 * Po[3 + j];
                    for (int j = 0; j < 3; ++j)
                        q[i * 3 + j] = gp[i * 3] * Vi[j] + gp[i * 3 + 1] * Vi[3 + j] + gp[i * 3 + 2] * Vi[6 + j];
                    if (co[i] < 0) continue;
                    // b: G^T r - G^T P V^-1 bp
                    bl[co[i]] += g0 * r[o * 2] + g1 * r[o * 2 + 1]
                                 - (q[i * 3] * bp[0] + q[i * 3 + 1] * bp[1] + q[i * 3 + 2] * bp[2]);
                    // S: G_o^T G_o
                    double* row = Sl + co[i] * NG;
                    for (int j = 0; j < K; ++j) if (co[j] >= 0) row[co[j]] += g0 * Go[j] + g1 * Go[K + j];
                }
            }
            // S: - sum_{o,o'} (G_o^T P_o V^-1) (P_o'^T G_o')
            for (int64_t k = 0; k < m; ++k) {
                const int64_t* co = cols + (o0 + k) * K; const double* q = &Q[k * K * 3];
                for (int64_t l = 0; l < m; ++l) {
                    const int64_t* cl = cols + (o0 + l) * K; const double* gp = &GP[l * K * 3];
                    for (int i = 0; i < K; ++i) {
                        if (co[i] < 0) continue;
                        double* row = Sl + co[i] * NG;
                        double q0 = q[i * 3], q1 = q[i * 3 + 1], q2 = q[i * 3 + 2];
                        for (int j = 0; j < K; ++j)
                            if (cl[j] >= 0) row[cl[j]] -= q0 * gp[j * 3] + q1 * gp[j * 3 + 1] + q2 * gp[j * 3 + 2];
                    }
                }
            }
        }
    }
    // reduce thread copies into S, b
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int64_t i = 0; i < NG; ++i)
        for (int t = 1; t < nthreads; ++t) {
            double* a = S + i * NG; const double* c = St[t] + i * NG;
            for (int64_t j = 0; j < NG; ++j) a[j] += c[j];
        }
    for (int t = 1; t < nthreads; ++t) {
        for (int64_t i = 0; i < NG; ++i) b[i] += bt[t][i];
        free(St[t]); free(bt[t]);
    }
    return 0;
}

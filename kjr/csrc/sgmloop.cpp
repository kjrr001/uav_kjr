// The step loop of the UAV library's SGM/MGM pass (matching_cost/sgm_cuda.py _run) in
// C++: the same `sgm_step` kernel with the same arguments, launched back to back on
// the caller's stream. The python loop spent ~15 us per launch against ~4.6 us of GPU
// work per step (67,890 launches per KJ-3 pair).
// Build: g++ -O2 -shared -fPIC sgmloop.cpp -o _sgmloop.so -lcuda
#include <cstdint>
#include <utility>

typedef int CUresult;
extern "C" CUresult cuLaunchKernel(void* f, unsigned gx, unsigned gy, unsigned gz,
                                   unsigned bx, unsigned by, unsigned bz, unsigned shmem,
                                   void* stream, void** params, void** extra);

// returns 0 and writes into *last the pointer that holds the last step's row
extern "C" int sgm_loop(void* fn, void* stream, unsigned grid, unsigned block, unsigned shmem,
                        void* vol, void* agg, void* prev, void* cur, int n_steps, int reverse,
                        long long s_step, long long s_perp, int N, int D, float p1, float p2,
                        int mode, int dx, int first, void** last) {
    int s = 0;
    long sStep = (long)s_step, sPerp = (long)s_perp;
    for (int i = 0; i < n_steps; ++i) {
        s = reverse ? n_steps - 1 - i : i;
        void* args[14] = {&vol, &agg, &prev, &cur, &s, &sStep, &sPerp, &N, &D, &p1, &p2,
                          &mode, &dx, &first};
        CUresult e = cuLaunchKernel(fn, grid, 1, 1, block, 1, 1, shmem, stream, args, nullptr);
        if (e != 0) return (int)e;
        first = 0;
        std::swap(prev, cur);
    }
    *last = prev;
    return 0;
}

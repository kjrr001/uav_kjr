"""Faster drop-in for the UAV library's census sub-pixel refiner (matching_cost/subpixel_cuda.py).

Profile of fsgbm_match on a KJ-3 pair (2636 x 5400 canvas, +-8 px, 2026-10-26): the
sub-pixel refiner `subpix_refine` was 79 % of the GPU time (2.17 of 2.74 s); the census
cost 2.5 %, the MGM passes 11.5 %. Per point the UAV kernel makes 8 E/G/H evaluations,
each with 25 support x (4 + 8 x 4) bilinear samples and 25 x 8 x 4 soft-ternary
descriptors (2 expf + 2 divisions each), i.e. ~13 k sigmoids per point.

This kernel computes EXACTLY the same numbers (same IEEE operations in the same order,
--fmad=false; checked bit-identical on KJ-3 pairs) with less work:
  * the points sit on the integer pixel grid and tau / eps are one value for the frame,
    so the LEFT descriptors sdesc(L(p + o_k) - L(p)) do not depend on the disparity:
    they are computed once per pixel (16 offset maps) instead of 8 times per point
    and 25 times per pixel (every point whose 5x5 support covers it);
  * the evaluations whose G / H are never read (E_0, the two validation-set E's) skip
    the +-h right samples and their descriptors; only their finiteness, which enters
    the `good` gate, is kept (bounds + finite samples, as sdesc of a finite difference
    is finite).
Falls back to the library kernel when those preconditions do not hold (non-integer
points, a per-point tau map, half != 2).

install() patches matching_cost.subpixel_cuda.refine, which census_subpixel.refine
imports at call time.
"""
import numpy as np

SOURCE = r"""
__device__ __forceinline__ float sigm(float x) { return 1.0f / (1.0f + expf(-x)); }
__device__ __forceinline__ float sdesc(float a, float tau, float eps) {
    return sigm((a - tau) / eps) - sigm((-a - tau) / eps);
}
// = kernels.py bsamp (bilinear, NaN outside)
__device__ __forceinline__ float bsamp(const float* img, int H, int W, float xx,
                                       float yy, bool* ok) {
    float x0 = floorf(xx), y0 = floorf(yy);
    bool o = (x0 >= 0.0f) && (x0 + 1.0f < (float)W) && (y0 >= 0.0f) && (y0 + 1.0f < (float)H);
    float xcf = fminf(fmaxf(x0, 0.0f), (float)(W - 2));
    float ycf = fminf(fmaxf(y0, 0.0f), (float)(H - 2));
    int xc = (int)xcf, yc = (int)ycf;
    float tx = xx - (float)xc, ty = yy - (float)yc;
    const float* r0 = img + (long)yc * W + xc;
    const float* r1 = r0 + W;
    float v = ((r0[0] * (1.0f - tx)) * (1.0f - ty) + (r0[1] * tx) * (1.0f - ty))
              + (r1[0] * (1.0f - tx)) * ty + (r1[1] * tx) * ty;
    *ok = o;
    return o ? v : __int_as_float(0x7fc00000);
}

// left descriptor maps: DL[k][i][j] = sdesc(L(j+dj_k, i+di_k) - L(j, i)), NaN where a
// sample is outside (as bsamp + sdesc give it)
extern "C" __global__ void subpix_ldesc(const float* L, int H, int W, float tau, float eps,
                                        const int* offs, int nk, float* DL) {
    long p = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= (long)H * W) return;
    int i = (int)(p / W), j = (int)(p % W);
    float px = (float)j, py = (float)i;
    bool o1, a1;
    float cL = bsamp(L, H, W, px, py, &o1);
    for (int k = 0; k < nk; ++k) {
        int di = offs[2 * k], dj = offs[2 * k + 1];
        float nL = bsamp(L, H, W, px + (float)dj, py + (float)di, &a1);
        float s = sdesc(nL - cL, tau, eps);
        DL[(long)k * H * W + p] = (o1 && a1 && isfinite(s)) ? s : __int_as_float(0x7fc00000);
    }
}

// 1-D form of bsamp for an integer row yy (ty = 0) and a finite image: the second-row
// terms of bsamp are then (r1 * ..) * 0 = +-0 and adding them leaves the sum unchanged,
// so this is bit-identical to bsamp and reads 2 samples instead of 4
__device__ __forceinline__ float bs1(const float* img, int H, int W, float xx, int row) {
    float x0 = floorf(xx);
    float xcf = fminf(fmaxf(x0, 0.0f), (float)(W - 2));
    int xc = (int)xcf;
    float tx = xx - (float)xc;
    const float* r0 = img + (long)row * W + xc;
    return (r0[0] * (1.0f - tx)) + (r0[1] * tx);
}
__device__ __forceinline__ bool in1(int H, int W, float xx, int row) {
    float x0 = floorf(xx);
    return (x0 >= 0.0f) && (x0 + 1.0f < (float)W) && (row >= 0) && (row + 1 < H);
}

// one evaluation at disparity d over offsets [k0, k0 + 8); gh = 0: E only.
// Same terms in the same summation order as kernels.py subpix_egh. With finite images
// a term is finite exactly when its samples are inside, so `good` is decided first from
// the bounds and the left maps (NaN = outside), then the terms are summed directly.
__device__ void egh2(const float* DL, const float* R, int H, int W, int xi, int yi,
                     float x, float d, float tau, float eps, const float* w,
                     const int* offs, int k0, int gh, float* E, float* G, float* Hs,
                     int* ngood) {
    const float h = 0.25f;
    const long HW = (long)H * W;
    float e_ = 0.0f, g_ = 0.0f, h_ = 0.0f;
    int ng = 0;
    int v = 0;
    for (int vy = -2; vy <= 2; ++vy) {
        int row = yi + vy;
        for (int vx = -2; vx <= 2; ++vx, ++v) {
            float px = x + (float)vx;
            float qx = px + d;
            const float* dl = DL + (long)k0 * HW + (long)row * W + (xi + vx);
            bool good = in1(H, W, qx, row) && in1(H, W, qx + h, row) && in1(H, W, qx - h, row);
            for (int k = 0; k < 8 && good; ++k) {
                int di = offs[2 * (k0 + k)], dj = offs[2 * (k0 + k) + 1];
                good = isfinite(dl[(long)k * HW]) && in1(H, W, qx + (float)dj, row + di)
                       && in1(H, W, (qx + h) + (float)dj, row + di)
                       && in1(H, W, (qx - h) + (float)dj, row + di);
            }
            if (!good) continue;
            ++ng;
            float wv = w[v];
            float cR = bs1(R, H, W, qx, row);
            float cP = 0.0f, cM = 0.0f;
            if (gh) { cP = bs1(R, H, W, qx + h, row); cM = bs1(R, H, W, qx - h, row); }
            for (int k = 0; k < 8; ++k) {
                int di = offs[2 * (k0 + k)], dj = offs[2 * (k0 + k) + 1];
                float sL = dl[(long)k * HW];
                float sR = sdesc(bs1(R, H, W, qx + (float)dj, row + di) - cR, tau, eps);
                float er = sL - sR;
                float r2 = sqrtf(er * er + 0.01f);
                e_ += (r2 - 0.1f) * wv;
                if (gh) {
                    float sP = sdesc(bs1(R, H, W, (qx + h) + (float)dj, row + di) - cP, tau, eps);
                    float sM = sdesc(bs1(R, H, W, (qx - h) + (float)dj, row + di) - cM, tau, eps);
                    float de = -(sP - sM) * 2.0f;
                    g_ += ((er / r2) * de) * wv;
                    h_ += ((de * de) / r2) * wv;
                }
            }
        }
    }
    *E = e_; *G = g_; *Hs = h_; *ngood = ng;
}

extern "C" __global__ void subpix_refine2(const float* DL, const float* R, int H, int W,
                                          const float* xs, const float* ys,
                                          const float* d0, float tau, float eps,
                                          const float* w, const int* offs, int n, int iters,
                                          float bound, float curv_min, float val_delta,
                                          float* d_out, unsigned char* ok_out,
                                          float* raw_out) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    float x = xs[i], y = ys[i], base = d0[i];
    int xi = (int)x, yi = (int)y;
    float E0, G, Hs, Ev0, Ef, Hf, Evf;
    int g0, gf, gv;
    egh2(DL, R, H, W, xi, yi, x, base, tau, eps, w, offs, 0, 0, &E0, &G, &Hs, &g0);
    egh2(DL, R, H, W, xi, yi, x, base, tau, eps, w, offs, 8, 0, &Ev0, &G, &Hs, &gv);
    float dv = base;
    for (int it = 0; it < iters; ++it) {
        float E;
        int gg;
        egh2(DL, R, H, W, xi, yi, x, dv, tau, eps, w, offs, 0, 1, &E, &G, &Hs, &gg);
        float step = Hs > curv_min ? -G / fmaxf(Hs, curv_min) : 0.0f;
        step = fminf(fmaxf(step, -0.25f), 0.25f);
        dv = fminf(fmaxf(dv + step, base - bound), base + bound);
    }
    egh2(DL, R, H, W, xi, yi, x, dv, tau, eps, w, offs, 0, 1, &Ef, &G, &Hf, &gf);
    egh2(DL, R, H, W, xi, yi, x, dv, tau, eps, w, offs, 8, 0, &Evf, &G, &Hs, &gv);
    const float inv25 = 1.0f / 25.0f;
    bool ok = ((float)g0 * inv25 > 0.75f) && ((float)gf * inv25 > 0.75f)
              && (Hf > curv_min) && (Ef <= E0 + 1e-6f)
              && (Evf <= (Ev0 + val_delta) + 1e-5f) && isfinite(dv)
              && (fabsf(dv - base) <= bound + 1e-6f);
    d_out[i] = ok ? dv : base;
    ok_out[i] = ok ? 1 : 0;
    raw_out[i] = dv;
}
"""

_ST = {"mod": None, "fns": {}}
_NP = {"p": np.uint64, "i": np.int32, "l": np.int64, "f": np.float32}


def _chk(res, what):
    if int(res[0]) != 0:
        raise RuntimeError(f"{what} failed: {res[0]}")
    return res[1] if len(res) == 2 else res[1:]


def _load():
    if _ST["mod"] is not None:
        return
    import torch
    from cuda.bindings import driver, nvrtc
    torch.empty(1, device="cuda")
    ma, mi = torch.cuda.get_device_capability()
    prog = _chk(nvrtc.nvrtcCreateProgram(SOURCE.encode(), b"kjr_fastmatch.cu", 0, [], []), "create")
    opts = [f"--gpu-architecture=sm_{ma}{mi}".encode(), b"--fmad=false", b"-std=c++17",
            b"--prec-div=true", b"--prec-sqrt=true"]          # = matching_cost/_nvrtc.py
    if int(nvrtc.nvrtcCompileProgram(prog, len(opts), opts)[0]) != 0:
        n = _chk(nvrtc.nvrtcGetProgramLogSize(prog), "log size")
        log = b" " * n
        nvrtc.nvrtcGetProgramLog(prog, log)
        raise RuntimeError("NVRTC: " + log.decode(errors="replace"))
    n = _chk(nvrtc.nvrtcGetCUBINSize(prog), "cubin size")
    cubin = b" " * n
    nvrtc.nvrtcGetCUBIN(prog, cubin)
    _ST["mod"] = _chk(driver.cuModuleLoadData(cubin), "cuModuleLoadData")
    _ST["drv"] = driver


def _launch(name, grid, block, sig, *args):
    import torch
    _load()
    f = _ST["fns"].get(name)
    if f is None:
        f = _ST["fns"][name] = _chk(_ST["drv"].cuModuleGetFunction(_ST["mod"], name.encode()), name)
    vals = [np.array([0 if a is None else a.data_ptr()], np.uint64) if k == "p" else np.array([a], _NP[k])
            for k, a in zip(sig, args)]
    ptrs = np.array([v.ctypes.data for v in vals], np.uint64)
    err = _ST["drv"].cuLaunchKernel(f, grid, 1, 1, block, 1, 1, 0,
                                    torch.cuda.current_stream().cuda_stream, ptrs.ctypes.data, 0)
    if int(err[0]) != 0:
        raise RuntimeError(f"cuLaunchKernel {name}: {err[0]}")


_orig = None


def refine(L, R, x, y, d_in, cw=7, tau=None, half=2, sigma=1.2, bound=0.5, iters=4,
           curv_min=1e-6, chunk=None, log=None, eps_k=0.5, val_delta=0.0, return_raw=False):
    """matching_cost.subpixel_cuda.refine, same interface and output."""
    import torch
    x = np.ascontiguousarray(np.asarray(x, np.float32))
    y = np.ascontiguousarray(np.asarray(y, np.float32))
    d0 = np.ascontiguousarray(np.asarray(d_in, np.float32))
    n = d0.size
    scalar_tau = tau is None or np.asarray(tau).size == 1
    if (int(half) != 2 or n == 0 or not scalar_tau or not np.isfinite(R).all()
            or not (np.array_equal(x, np.round(x)) and np.array_equal(y, np.round(y)))):
        return _orig(L, R, x, y, d_in, cw=cw, tau=tau, half=half, sigma=sigma, bound=bound,
                     iters=iters, curv_min=curv_min, chunk=chunk, log=log, eps_k=eps_k,
                     val_delta=val_delta, return_raw=return_raw)
    from matching_cost.subpixel_torch import _FIT8, _VAL8
    dev = torch.device("cuda")
    Lt = torch.from_numpy(np.ascontiguousarray(L, np.float32)).to(dev)
    Rt = torch.from_numpy(np.ascontiguousarray(R, np.float32)).to(dev)
    # tau / eps exactly as the library prepares them (one value here)
    if tau is None:
        _fb = 0.2 * float(np.nanstd(L))
        t = _fb if np.isfinite(_fb) and _fb > 0 else 1e-6
    else:
        t = float(np.asarray(tau, np.float32).ravel()[0])
        if not np.isfinite(t) or t <= 0:
            _fb = 0.2 * float(np.nanstd(L))
            t = _fb if np.isfinite(_fb) and _fb > 0 else 1e-6
    t = np.float32(t)
    e = np.maximum(float(eps_k) * np.full(1, t, np.float32), 1e-6).astype(np.float32)[0]
    vv, uu = np.mgrid[-half:half + 1, -half:half + 1].astype(np.float32)
    w = np.exp(-0.5 * (uu ** 2 + vv ** 2) / sigma ** 2).ravel(); w /= w.sum()
    offs = np.asarray(_FIT8 + _VAL8, np.int32).ravel()
    H, W = Lt.shape
    T = lambda a: torch.from_numpy(np.ascontiguousarray(a)).to(dev)
    wt, ot = T(w.astype(np.float32)), T(offs)
    DL = torch.empty((16, H, W), dtype=torch.float32, device=dev)
    _launch("subpix_ldesc", (H * W + 255) // 256, 256, "piiffpip", Lt, H, W, float(t), float(e), ot, 16, DL)
    xs, ys, dd = T(x), T(y), T(d0)
    d_out = torch.empty(n, dtype=torch.float32, device=dev)
    raw = torch.empty_like(d_out)
    ok = torch.empty(n, dtype=torch.uint8, device=dev)
    _launch("subpix_refine2", (n + 127) // 128, 128, "ppiipppffppiifffppp",
            DL, Rt, H, W, xs, ys, dd, float(t), float(e), wt, ot, int(n), int(iters),
            float(np.float32(bound)), float(np.float32(curv_min)), float(np.float32(val_delta)),
            d_out, ok, raw)
    del DL
    d_out = d_out.cpu().numpy(); okn = ok.cpu().numpy().astype(bool); raw = raw.cpu().numpy()
    if log is not None:
        mv = np.abs(d_out - d0)[okn]
        log(f"    census-subpixel (kjr cuda): {okn.sum()}/{n} refined ({100*okn.mean():.1f} %), "
            f"|delta| median {np.median(mv) if mv.size else 0:.4f} px")
    return (d_out, okn, raw) if return_raw else (d_out, okn)


def install():
    """route the UAV library's CUDA sub-pixel refiner through refine() above"""
    global _orig
    import matching_cost.subpixel_cuda as SC
    if SC.refine is not refine:
        _orig = SC.refine
        SC.refine = refine


# ---------------------------------------------------------------- SGM step loop in C++
def _sgmloop():
    import ctypes, os
    if "sgm" in _ST:
        return _ST["sgm"]
    so = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_sgmloop.so")
    if not os.path.exists(so):
        src = os.path.join(os.path.dirname(so), "csrc", "sgmloop.cpp")
        if os.system(f"g++ -O2 -shared -fPIC {src} -o {so} -lcuda") != 0:
            _ST["sgm"] = None
            return None
    lib = ctypes.CDLL(so)
    c = ctypes
    lib.sgm_loop.argtypes = [c.c_void_p, c.c_void_p, c.c_uint, c.c_uint, c.c_uint, c.c_void_p, c.c_void_p,
                             c.c_void_p, c.c_void_p, c.c_int, c.c_int, c.c_longlong, c.c_longlong, c.c_int,
                             c.c_int, c.c_float, c.c_float, c.c_int, c.c_int, c.c_int,
                             c.POINTER(c.c_void_p)]
    lib.sgm_loop.restype = c.c_int
    _ST["sgm"] = lib
    return lib


def _sgm_run(vol, agg, p1, p2, n_steps, N, s_step, s_perp, mode, dx, reverse, state):
    """matching_cost.sgm_cuda._run with the step loop in C++ (same kernel, same order)"""
    import ctypes
    import torch
    from matching_cost import _nvrtc as NV
    from matching_cost import sgm_cuda as SC
    H, W, D = vol.shape
    prev = (state.float().contiguous().clone() if state is not None
            else torch.empty((N, D), dtype=torch.float32, device=vol.device))
    cur = torch.empty_like(prev)
    if not NV.available():
        raise RuntimeError(f"CUDA kernels unavailable: {NV.error()}")
    f = NV.fn("sgm_step")
    last = ctypes.c_void_p()
    shmem = ((SC.TK + 2) * D + (SC.TK + 2)) * 4
    e = _sgmloop().sgm_loop(int(f), torch.cuda.current_stream().cuda_stream, (N + SC.TK - 1) // SC.TK,
                            SC._BLOCK, shmem, vol.data_ptr(), agg.data_ptr(), prev.data_ptr(),
                            cur.data_ptr(), int(n_steps), 1 if reverse else 0, int(s_step), int(s_perp),
                            int(N), int(D), float(p1), float(p2), int(mode), int(dx),
                            1 if state is None else 0, ctypes.byref(last))
    if e != 0:
        raise RuntimeError(f"sgm_loop: cuLaunchKernel error {e}")
    return prev if last.value == prev.data_ptr() else cur


def install_sgm():
    from matching_cost import sgm_cuda as SC
    if _sgmloop() is not None and SC._run is not _sgm_run:
        SC._run = _sgm_run

#!/usr/bin/env python
"""parallel_kernels.py -- multi-core versions of JUICE's compiled reaction kernels (2026-11-14).

The serial kernels (`scenario_runner._stoich_eval_numba`, `twostep_solver.pl_kernel`) walk the stoichiometry
entries in REACTION order and scatter into the species arrays.  For large mechanisms (GECKO-A: 4e5 - 1e6
reactions) that loop is limited by memory access and takes > 50 % of a TWOSTEP run.

Here the same entries are stored ordered by SPECIES (CSR): the production, loss or net rate of one species is
then an independent sum, and the species loop (and the reaction loop for the velocities) runs on several
cores (numba prange).  The stable sort keeps the order of the terms of each sum, and every term is computed
with the same operations as in the serial kernels, so the results are bit-identical to them.

    lay = build_layout(n_species, n_rxn, ode_r_rxn, ode_r_sp, ode_r_coeff, env_r_rxn, env_r_name_idx,
                       env_r_coeff, stoich_rxn, stoich_sp, stoich_sc)
    set_threads(4)
    rhs_par(k, y, env, *lay['v'], *lay['s'], v, dydt)        # dydt assigned (not accumulated)
    vel_par(k, y, env, *lay['v'], v)                         # velocities only (Jacobian)
    pl_par(k, y, env, *lay['v'], *lay['p'], *lay['l'], v, P, L)

    python parallel_kernels.py --selftest [threads]          # serial vs parallel on a random mechanism
"""
from __future__ import annotations

import numpy as np

try:
    import numba as _numba
    from numba import prange as _prange
    NUMBA = True
except ImportError:          # pragma: no cover
    NUMBA = False


def set_threads(n):
    """number of cores used by the parallel kernels; returns the number actually set"""
    n = max(1, min(int(n), int(_numba.config.NUMBA_NUM_THREADS)))
    _numba.set_num_threads(n)
    return n


def _ptr(idx, n):
    """CSR pointer for the sorted index array idx (values 0..n-1)"""
    p = np.zeros(n + 1, dtype=np.int64)
    if len(idx):
        np.add.at(p, np.asarray(idx, dtype=np.int64) + 1, 1)
    return np.cumsum(p)


def _bins(work, n_bins):
    """Deal the species out to n_bins cores so that every core gets the same amount of `work`.
    -> (order, borders): core c handles the species order[borders[c]:borders[c+1]] (ascending inside a core).
    Equal NUMBERS of species per core would give the core with the inorganic species (OH, NO, HO2, NO2: tens of
    thousands of terms each) most of the work: 62 % on the first of 4 cores for a GECKO-A limonene mechanism.
    Largest species first, each to the core with the least work so far."""
    import heapq
    work = np.asarray(work, dtype=np.float64)
    n = len(work)
    n_bins = max(1, min(int(n_bins), max(n, 1)))
    members = [[] for _ in range(n_bins)]
    heap = [(0.0, c) for c in range(n_bins)]
    for sidx in np.argsort(-work, kind="stable").tolist():
        load, c = heapq.heappop(heap)
        members[c].append(sidx)
        heapq.heappush(heap, (load + work[sidx], c))
    order = np.array([x for mem in members for x in sorted(mem)], dtype=np.int64)
    borders = np.concatenate(([0], np.cumsum([len(mem) for mem in members]))).astype(np.int64)
    return order, borders


def build_layout(n_species, n_rxn, ode_r_rxn, ode_r_sp, ode_r_coeff, env_r_rxn, env_r_name_idx, env_r_coeff,
                 stoich_rxn, stoich_sp, stoich_sc, n_threads=4):
    """Species-ordered copies of the runner's stoichiometry arrays.
    Returns {'v': (rp, r_sp, r_co, ep, e_idx, e_co), 's': (sp, s_rxn, s_sc, s_order, s_borders),
             'p': (pp, p_rxn, p_sc), 'l': (lp, l_rxn, l_sc, pl_order, pl_borders)}  (order / borders: which species
    each of the n_threads cores handles, equal work per core, see _bins)  or None if the reactant arrays are not in reaction
    order (then the serial kernels must be used)."""
    ode_r_rxn = np.asarray(ode_r_rxn); env_r_rxn = np.asarray(env_r_rxn)
    if (len(ode_r_rxn) > 1 and np.any(np.diff(ode_r_rxn) < 0)) or (len(env_r_rxn) > 1 and np.any(np.diff(env_r_rxn) < 0)):
        return None
    v = (_ptr(ode_r_rxn, n_rxn), np.ascontiguousarray(ode_r_sp, dtype=np.int32), np.ascontiguousarray(ode_r_coeff, dtype=np.float64),
         _ptr(env_r_rxn, n_rxn), np.ascontiguousarray(env_r_name_idx, dtype=np.int32), np.ascontiguousarray(env_r_coeff, dtype=np.float64))
    stoich_rxn = np.asarray(stoich_rxn); stoich_sp = np.asarray(stoich_sp); stoich_sc = np.asarray(stoich_sc, dtype=np.float64)

    def by_species(mask, sign):
        sp, rx, sc = stoich_sp[mask], stoich_rxn[mask], stoich_sc[mask]
        o = np.argsort(sp, kind="stable")                 # stable: terms of one species stay in reaction order
        return _ptr(sp[o], n_species), np.ascontiguousarray(rx[o], dtype=np.int32), np.ascontiguousarray(sign * sc[o])
    all_ = np.ones(len(stoich_sp), dtype=bool)
    lay_s, lay_p, lay_l = by_species(all_, 1.0), by_species(stoich_sc > 0.0, 1.0), by_species(~(stoich_sc > 0.0), -1.0)
    # work per species: one unit per term (a loss term has a division: two units) plus one for the species itself
    w_s = np.diff(lay_s[0]) + 1.0
    w_pl = np.diff(lay_p[0]) + 2.0 * np.diff(lay_l[0]) + 1.0
    return {'v': v, 's': lay_s + _bins(w_s, n_threads), 'p': lay_p, 'l': lay_l + _bins(w_pl, n_threads)}


if NUMBA:
    @_numba.njit(cache=True, parallel=True)
    def vel_par(k_buf, y, env_vals, rp, r_sp, r_co, ep, e_idx, e_co, v_buf):
        """v[i] = k[i] * prod(y^c) * prod(env^c)   (same operations and order as the serial kernels)"""
        for i in _prange(len(k_buf)):
            x = k_buf[i]
            if x != 0.0:
                for j in range(rp[i], rp[i + 1]):
                    c = r_co[j]
                    x *= y[r_sp[j]] if c == 1.0 else y[r_sp[j]] ** c      # x**1.0 == x exactly
                    if x == 0.0:
                        break
            if x != 0.0:
                for j in range(ep[i], ep[i + 1]):
                    c = e_co[j]
                    x *= env_vals[e_idx[j]] if c == 1.0 else env_vals[e_idx[j]] ** c
                    if x == 0.0:
                        break
            v_buf[i] = x

    @_numba.njit(cache=True, parallel=True)
    def rhs_par(k_buf, y, env_vals, rp, r_sp, r_co, ep, e_idx, e_co, sp, s_rxn, s_sc, so, sb, v_buf, dydt_out):
        """dydt_out[s] = sum(v * signed coefficient)   (assigned; replaces _stoich_eval_numba)"""
        for i in _prange(len(k_buf)):
            x = k_buf[i]
            if x != 0.0:
                for j in range(rp[i], rp[i + 1]):
                    c = r_co[j]
                    x *= y[r_sp[j]] if c == 1.0 else y[r_sp[j]] ** c      # x**1.0 == x exactly
                    if x == 0.0:
                        break
            if x != 0.0:
                for j in range(ep[i], ep[i + 1]):
                    c = e_co[j]
                    x *= env_vals[e_idx[j]] if c == 1.0 else env_vals[e_idx[j]] ** c
                    if x == 0.0:
                        break
            v_buf[i] = x
        for c in _prange(len(sb) - 1):          # one bin of species per core, equal work, see _bins
            for q in range(sb[c], sb[c + 1]):
                s = so[q]
                a = 0.0
                for j in range(sp[s], sp[s + 1]):
                    a += v_buf[s_rxn[j]] * s_sc[j]
                dydt_out[s] = a

    @_numba.njit(cache=True, parallel=True)
    def pl_par(k_buf, y, env_vals, rp, r_sp, r_co, ep, e_idx, e_co, pp, p_rxn, p_sc, lp, l_rxn, l_sc, plo, plb, v_buf, P_out, L_out):
        """P_out (molec cm-3 s-1) and L_out (s-1), assigned; replaces twostep_solver.pl_kernel (y > 0)"""
        for i in _prange(len(k_buf)):
            x = k_buf[i]
            if x != 0.0:
                for j in range(rp[i], rp[i + 1]):
                    c = r_co[j]
                    x *= y[r_sp[j]] if c == 1.0 else y[r_sp[j]] ** c      # x**1.0 == x exactly
                    if x == 0.0:
                        break
            if x != 0.0:
                for j in range(ep[i], ep[i + 1]):
                    c = e_co[j]
                    x *= env_vals[e_idx[j]] if c == 1.0 else env_vals[e_idx[j]] ** c
                    if x == 0.0:
                        break
            v_buf[i] = x
        for c in _prange(len(plb) - 1):         # one bin of species per core, equal work, see _bins
            for q in range(plb[c], plb[c + 1]):
                s = plo[q]
                a = 0.0
                for j in range(pp[s], pp[s + 1]):
                    a += v_buf[p_rxn[j]] * p_sc[j]
                P_out[s] = a
                b = 0.0
                ys = y[s]
                for j in range(lp[s], lp[s + 1]):
                    vv = v_buf[l_rxn[j]]
                    if vv != 0.0:
                        b += l_sc[j] * vv / ys
                L_out[s] = b


def _selftest(threads=4):
    import time
    from juice.core.scenario_runner import _stoich_eval_numba
    from juice.core.twostep_solver import pl_kernel
    rng = np.random.default_rng(1)
    NS, NR = 92341, 444980
    nrea = rng.choice([1, 2], NR, p=[0.45, 0.55]); npro = rng.integers(1, 5, NR)
    r_rxn = np.repeat(np.arange(NR), nrea).astype(np.int32); r_sp = rng.integers(0, NS, len(r_rxn)).astype(np.int32)
    hub = rng.random(len(r_sp)) < 0.25                      # a quarter of all reactant terms on 20 hub species at the front
    r_sp[hub] = rng.integers(0, 20, int(hub.sum())).astype(np.int32)
    r_co = rng.choice([1.0, 2.0], len(r_rxn), p=[0.97, 0.03])
    e_rxn = np.sort(rng.choice(NR, NR // 10, replace=False)).astype(np.int32); e_idx = rng.integers(0, 3, len(e_rxn)).astype(np.int32)
    e_co = np.ones(len(e_rxn)); env = np.array([2.5e19, 5e18, 1.0])
    p_rxn = np.repeat(np.arange(NR), npro); p_sp = rng.integers(0, NS, len(p_rxn)); p_sc = rng.uniform(0.1, 1, len(p_rxn))
    hubp = rng.random(len(p_sp)) < 0.2
    p_sp[hubp] = rng.integers(0, 20, int(hubp.sum()))
    s_rxn = np.concatenate([r_rxn, p_rxn]); s_sp = np.concatenate([r_sp, p_sp]); s_sc = np.concatenate([-r_co, p_sc])
    o = np.argsort(s_rxn, kind="stable"); s_rxn, s_sp, s_sc = s_rxn[o].astype(np.int32), s_sp[o].astype(np.int32), s_sc[o]
    k = rng.uniform(1e-14, 1e-10, NR); k[rng.choice(NR, 2000, replace=False)] = 0.0; y = rng.uniform(1e3, 1e10, NS)
    n = set_threads(threads)
    lay = build_layout(NS, NR, r_rxn, r_sp, r_co, e_rxn, e_idx, e_co, s_rxn, s_sp, s_sc, n_threads=n)
    wpl = np.diff(lay['p'][0]) + 2.0 * np.diff(lay['l'][0]) + 1.0
    eq = np.linspace(0, NS, n + 1).astype(int)
    oo, bb = lay['l'][3], lay['l'][4]
    per = [wpl[oo[bb[c]:bb[c + 1]]].sum() / wpl.sum() for c in range(len(bb) - 1)]
    print(f"work share of the busiest core: equal species counts {100 * max(wpl[eq[c]:eq[c + 1]].sum() for c in range(n)) / wpl.sum():.0f} %, "
          f"dealt out by work {100 * max(per):.1f} % (ideal {100 / n:.1f} %)")
    v0 = np.zeros(NR); d0 = np.zeros(NS); P0 = np.zeros(NS); L0 = np.zeros(NS)
    v1 = np.zeros(NR); d1 = np.zeros(NS); P1 = np.zeros(NS); L1 = np.zeros(NS); v2 = np.zeros(NR)

    def bench(f, nrep=15):
        f(); t = time.perf_counter()
        for _ in range(nrep): f()
        return (time.perf_counter() - t) / nrep * 1e3

    def ser_rhs():
        d0[:] = 0.0; _stoich_eval_numba(k, y, env, r_rxn, r_sp, r_co, e_rxn, e_idx, e_co, s_rxn, s_sp, s_sc, v0, d0)

    def ser_pl():
        P0[:] = 0.0; L0[:] = 0.0; pl_kernel(k, y, env, r_rxn, r_sp, r_co, e_rxn, e_idx, e_co, s_rxn, s_sp, s_sc, v0, P0, L0)
    t1 = bench(ser_rhs); t2 = bench(lambda: rhs_par(k, y, env, *lay['v'], *lay['s'], v1, d1))
    print(f"RHS      serial {t1:6.2f} ms   {n} threads {t2:6.2f} ms   identical: {np.array_equal(d0, d1)} (v {np.array_equal(v0, v1)})")
    t1 = bench(ser_pl); t2 = bench(lambda: pl_par(k, y, env, *lay['v'], *lay['p'], *lay['l'], v1, P1, L1))
    print(f"P/L      serial {t1:6.2f} ms   {n} threads {t2:6.2f} ms   identical: P {np.array_equal(P0, P1)}  L {np.array_equal(L0, L1)}")
    t2 = bench(lambda: vel_par(k, y, env, *lay['v'], v2))
    print(f"velocities only        {n} threads {t2:6.2f} ms   identical: {np.array_equal(v0, v2)}")


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        _selftest(int(sys.argv[-1]) if sys.argv[-1].isdigit() else 4)

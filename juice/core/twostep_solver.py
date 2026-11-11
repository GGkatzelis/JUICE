#!/usr/bin/env python
"""twostep_solver.py -- TWOSTEP (Verwer 1994) as a JUICE ODE backend.

Written by F. Rohrer (MCM_SAPHIR_Tool branch easy-start-gui; 2026-11-09, adaptive iteration 279df43) and taken over
unchanged into JUICE; the runner (scenario_runner._solve, solver_method "TWOSTEP") and the solver profile "twostep"
(scenario_prep) drive it.  Integration notes: Franz's TWOSTEP_BACKEND.md; here README "Solver profiles".

    python -m juice.core.twostep_solver --selftest          # mini SAPHIR-like mechanism vs scipy BDF reference

TWOSTEP integrates   dy/dt = P(t, y) - L(t, y) * y   with the variable-step 2-step BDF formula whose implicit
stage is approximated by a FIXED number (numit) of fixed-point iterations instead of Newton + LU:

        y_new <- max(smallc, (w + g*P(y_new)) / (1 + g*L(y_new)))

(w = BDF2 history combination, g = scaled step).  No Jacobian, no linear algebra: the cost per step is numit
production/loss evaluations, i.e. ~numit x the cost of ONE RHS call.  This is the solver of GECKO-A's
boxmodel4gecko (solver.f90 + iteration.f90, transcribed here line by line, see design note
`reports/TWOSTEP Solver Modul fuer JUICE.md`).

Two parts:
  1. `pl_kernel`  -- @njit production/loss kernel on EXACTLY the stoichiometry arrays scenario_runner._solve
     already builds for `_stoich_eval_numba` (ode_r_*, env_r_*, stoich_*), so the runner can call it with
     its existing k_buf / env_vals.  P_i = sum of product terms, L_i = sum of reactant terms / y_i.
  2. `twostep`    -- the integrator (startup implicit Euler, BDF2 with variable ratio, Verwer's error
     estimate, NEWDT/FIT step control, restart logic, fixed-step mode when dtmin == dtmax).

Iteration type: boxmodel4gecko's `iterfirst` evaluates ALL rates from the previous iterate and then updates
all species at once -- that is a JACOBI sweep, not Gauss-Seidel (despite the comment in the Fortran).  This
prototype does the same (vectorised, order-independent).  A true Gauss-Seidel sweep (species-ordered,
incremental rate updates) is a possible later refinement, see the design note.

Exchange pairs (JUICE, 10 Nov 2026): a Jacobi sweep cannot converge on a pair A <-> B that interconverts by
first-order reactions much faster than 1/step (complete MCM: CH3SO + O2 <-> CH3SOO2, ~3.5e5 s-1 each way); the
error test then fails at every step and TWOSTEP sits at dtmin.  `find_exchange_pairs` detects such pairs in the
stoichiometry and `PLSystem.pair_solve` solves each pair's 2x2 implicit system exactly inside every sweep.
Mechanisms without such pairs are integrated bit for bit as before; PLSystem(..., pairs=False) restores the plain
sweep everywhere.
"""
from __future__ import annotations

import math
import time

import numpy as np

try:
    import numba as _numba
    _NJIT = _numba.njit(cache=True)
    NUMBA = True
except ImportError:          # pragma: no cover
    def _NJIT(f):
        return f
    NUMBA = False

SMALLC = 1.0e-30            # concentration floor (boxmodel4gecko: 1e-40); molecules cm-3


# ----------------------------------------------------------------------------------------------------------
# production / loss kernel on JUICE's stoichiometry arrays
# ----------------------------------------------------------------------------------------------------------
@_NJIT
def pl_kernel(k_buf, y, env_vals,
              ode_r_rxn, ode_r_sp, ode_r_coeff,
              env_r_rxn, env_r_name_idx, env_r_coeff,
              stoich_rxn, stoich_sp, stoich_sc,
              v_buf, P_out, L_out):
    """Fill P_out (molec cm-3 s-1) and L_out (s-1) for state y (must be > 0).  Same inputs as
    scenario_runner._stoich_eval_numba; v_buf is scratch, P_out/L_out must be pre-zeroed."""
    for i in range(len(k_buf)):
        v_buf[i] = k_buf[i]
    for j in range(len(ode_r_rxn)):
        i = ode_r_rxn[j]
        if v_buf[i] != 0.0:
            v_buf[i] *= y[ode_r_sp[j]] ** ode_r_coeff[j]
    for j in range(len(env_r_rxn)):
        i = env_r_rxn[j]
        if v_buf[i] != 0.0:
            v_buf[i] *= env_vals[env_r_name_idx[j]] ** env_r_coeff[j]
    for j in range(len(stoich_rxn)):
        sc = stoich_sc[j]
        v = v_buf[stoich_rxn[j]]
        s = stoich_sp[j]
        if sc > 0.0:
            P_out[s] += v * sc
        elif v != 0.0:
            L_out[s] += (-sc) * v / y[s]        # reactant: linearised loss frequency (coeff 2 -> 2 k y)


def find_exchange_pairs(arrays, n_species):
    """Reversible first-order pairs A <-> B in the stoichiometry arrays.

    A reaction r is a first-order exchange A -> B when its only ODE reactant is A with exponent 1 (env
    reactants such as O2 or M are allowed), A's net coefficient is -1 and B's is +1.  A pair needs both
    directions.  Each species joins at most one pair (greedy, in species order).  Returns
    (ia, ib, rx_ab, seg_ab, rx_ba, seg_ba): pair species, and the reactions of each direction with the
    pair index they belong to."""
    (ode_r_rxn, ode_r_sp, ode_r_coeff, _er, _en, _ec, stoich_rxn, stoich_sp, stoich_sc) = arrays
    reac = {}
    for r, sp, c in zip(np.asarray(ode_r_rxn).tolist(), np.asarray(ode_r_sp).tolist(), np.asarray(ode_r_coeff).tolist()):
        reac.setdefault(r, []).append((sp, c))
    net = {}
    for r, sp, sc in zip(np.asarray(stoich_rxn).tolist(), np.asarray(stoich_sp).tolist(), np.asarray(stoich_sc).tolist()):
        net.setdefault(r, {})
        net[r][sp] = net[r].get(sp, 0.0) + sc
    edges = {}                                     # (a, b) -> [reactions a -> b]
    for r, lst in reac.items():
        if len(lst) != 1 or lst[0][1] != 1.0:
            continue
        a = lst[0][0]
        st = net.get(r, {})
        if st.get(a, 0.0) != -1.0:
            continue
        for b, sc in st.items():
            if b != a and sc == 1.0:
                edges.setdefault((a, b), []).append(r)
    used = set()
    ia, ib, rx_ab, seg_ab, rx_ba, seg_ba = [], [], [], [], [], []
    for (a, b) in sorted(edges):
        if a >= b or (b, a) not in edges or a in used or b in used:
            continue
        k = len(ia)
        ia.append(a); ib.append(b); used.update((a, b))
        rx_ab += edges[(a, b)]; seg_ab += [k] * len(edges[(a, b)])
        rx_ba += edges[(b, a)]; seg_ba += [k] * len(edges[(b, a)])
    as_i = lambda v: np.asarray(v, dtype=np.int64)
    return as_i(ia), as_i(ib), as_i(rx_ab), as_i(seg_ab), as_i(rx_ba), as_i(seg_ba)


@_NJIT
def pair_kernel(ia, ib, rx_ab, seg_ab, rx_ba, seg_ba, v_buf, P, L, y, w, g, ynext, smallc, vab, vba):
    """Exact 2x2 implicit solve of the exchange pairs (see PLSystem.pair_solve); vab/vba are scratch."""
    n = ia.size
    for k in range(n):
        vab[k] = 0.0
        vba[k] = 0.0
    for j in range(rx_ab.size):
        vab[seg_ab[j]] += v_buf[rx_ab[j]]
    for j in range(rx_ba.size):
        vba[seg_ba[j]] += v_buf[rx_ba[j]]
    for k in range(n):
        a = ia[k]
        b = ib[k]
        kab = vab[k] / y[a]
        kba = vba[k] / y[b]
        if kab <= 0.0 and kba <= 0.0:
            continue
        if L[a] < kab * (1.0 - 1e-9) or L[b] < kba * (1.0 - 1e-9):
            continue                                   # held species (P = L = 0) or inconsistent: leave alone
        ra = w[a] + g * max(P[a] - vba[k], 0.0)
        rb = w[b] + g * max(P[b] - vab[k], 0.0)
        Aa = 1.0 + g * L[a]
        Ab = 1.0 + g * L[b]
        D = Aa * Ab - (g * g) * kab * kba
        ynext[a] = max(smallc, (ra * Ab + g * kba * rb) / D)
        ynext[b] = max(smallc, (rb * Aa + g * kab * ra) / D)


class PLSystem:
    """Convenience wrapper: production/loss function from JUICE-style arrays + optional extra terms.

    extra_P(t, y, P) / extra_L(t, y, L) may add continuous sources, bursts, analytical dilution ...
    k_buf and env_vals are OWNED by the caller (the runner refreshes them per ENZ row)."""

    def __init__(self, n_species, k_buf, env_vals, arrays, extra_P=None, extra_L=None, pairs=True):
        (self.ode_r_rxn, self.ode_r_sp, self.ode_r_coeff,
         self.env_r_rxn, self.env_r_name_idx, self.env_r_coeff,
         self.stoich_rxn, self.stoich_sp, self.stoich_sc) = arrays
        self.k_buf = k_buf
        self.env_vals = env_vals
        self.v_buf = np.zeros(len(k_buf))
        self.P = np.zeros(n_species)
        self.L = np.zeros(n_species)
        self.extra_P = extra_P
        self.extra_L = extra_L
        self.ncall = 0
        # reversible first-order pairs A <-> B, solved together inside each TWOSTEP iteration (pair_solve)
        self.pairs = find_exchange_pairs(arrays, n_species) if pairs else None
        if self.pairs is not None and self.pairs[0].size == 0:
            self.pairs = None
        if self.pairs is not None:
            self._vab = np.zeros(self.pairs[0].size)
            self._vba = np.zeros(self.pairs[0].size)

    def __call__(self, t, y):
        self.P[:] = 0.0
        self.L[:] = 0.0
        pl_kernel(self.k_buf, y, self.env_vals,
                  self.ode_r_rxn, self.ode_r_sp, self.ode_r_coeff,
                  self.env_r_rxn, self.env_r_name_idx, self.env_r_coeff,
                  self.stoich_rxn, self.stoich_sp, self.stoich_sc,
                  self.v_buf, self.P, self.L)
        if self.extra_P is not None:
            self.extra_P(t, y, self.P)
        if self.extra_L is not None:
            self.extra_L(t, y, self.L)
        self.ncall += 1
        return self.P, self.L

    def pair_solve(self, y, w, g, ynext):
        """Overwrite ynext for the exchange pairs with the exact solution of the 2x2 implicit system.

        The Jacobi sweep y_new = (w + g P(y)) / (1 + g L(y)) takes A's production from B at the OLD iterate;
        for a pair exchanging at k >> 1/g each sweep only moves the pair by a factor ~1/(g k), so it never
        converges (DMS: CH3SO + O2 <-> CH3SOO2, ~3.5e5 s-1 each way).  With P and L from the last call at
        state y, the pair's own terms are taken out of P and the two linear equations
            a (1 + g La) - g kba b = w_a + g (Pa - vba)
            b (1 + g Lb) - g kab a = w_b + g (Pb - vab)
        are solved exactly (La >= kab and Lb >= kba, so the determinant is positive).  Every other term
        stays at the iterate, as in the Jacobi sweep.  A pair whose L does not contain its exchange rate
        (held species: P = L = 0) is left alone."""
        if self.pairs is None:
            return
        ia, ib, rx_ab, seg_ab, rx_ba, seg_ba = self.pairs
        pair_kernel(ia, ib, rx_ab, seg_ab, rx_ba, seg_ba, self.v_buf, self.P, self.L, y,
                    np.asarray(w, dtype=np.float64), float(g), ynext, SMALLC, self._vab, self._vba)

    def rhs(self, t, y):
        """dy/dt for scipy (reference runs)"""
        y = np.maximum(y, SMALLC)
        P, L = self(t, y)
        return P - L * y


# ----------------------------------------------------------------------------------------------------------
# the integrator (solver.f90 TWOSTEP / FIT / NEWDT)
# ----------------------------------------------------------------------------------------------------------
def _fit(t, te, dt):
    """shrink dt so that the remaining interval is an integer multiple of dt (as soon as <= 10 steps remain)"""
    rns = (te - t) / dt
    if rns > 10.0:
        return dt
    ns = int(rns) + 1
    return (te - t) / ns


def twostep(pl_fn, y0, t0, t1, dtmin, dtmax, atol, rtol, numit=2, smallc=SMALLC,
            max_steps=2_000_000, stats=None, trace=None,
            numit_max=None, conv=0.1, max_growth=2.0, restart_growth=None, restart_steps=0):
    """Integrate dy/dt = P - L*y from t0 to t1 with Verwer's TWOSTEP.

    pl_fn(t, y) -> (P, L)   arrays (molec cm-3 s-1, s-1); y is always >= smallc.
    atol        float or array (molec cm-3), rtol float; error test = max-norm of the BDF2 local error estimate.
    dtmin == dtmax -> fixed step (no error control, last step shortened to hit t1).
    Returns y(t1).  `stats` (dict) collects nfcn, naccpt, nrejec, nstart.
    `trace(t_abs, y, dt, kind)` is called after every accepted step (kind 'euler' or 'bdf2') -- debugging only.

    Iteration control: `numit` sweeps are always done; when `numit_max` > numit, further sweeps follow until the
    change of the iterate is below `conv` x (atol + rtol|y|) in the max-norm or numit_max is reached (an unconverged
    Jacobi iterate pollutes the BDF2 error estimate and forces sub-second steps -- intermediate-lifetime RO2).
    Step growth: the factor per step is limited to `max_growth` (Verwer: 2); `restart_growth` (if set) applies for
    the first `restart_steps` BDF2 steps after a (re)start -- the parasitic root of variable-step BDF2 is
    omega^2/(1+2 omega) (0.8 for doubling, 1/3 for constant steps), so fast growth right after the implicit-Euler
    start carries the transient of the start step into the following steps (floor species end at ~20 x atol)."""
    y = np.maximum(np.asarray(y0, dtype=np.float64), smallc)
    n = y.size
    atol_v = np.broadcast_to(np.asarray(atol, dtype=np.float64), (n,))
    if numit_max is None or numit_max < numit:
        numit_max = numit
    # Integrate in time RELATIVE to t0: JUICE times are seconds since 2000 (~8e8 s), where an absolute
    # end test such as t >= te - 1e-9*|te| would stop up to ~1 s early (1 % of a 60-s ENZ row).
    t0_abs = float(t0)
    t = 0.0
    te = float(t1) - t0_abs
    nfcn = naccpt = nrejec = nstart = 0
    failer = restart = False
    eps_end = 1e-9 * max(1e-3, te)

    pair_solve = getattr(pl_fn, "pair_solve", None)

    def iterate(w, g, ynew, t_eval):
        nonlocal nfcn
        for it in range(numit_max):
            P, L = pl_fn(t0_abs + t_eval, ynew)
            nfcn += 1
            ynext = np.maximum(smallc, (w + g * P) / (1.0 + g * L))
            if pair_solve is not None:
                pair_solve(ynew, w, g, ynext)      # fast reversible pairs A <-> B solved exactly
            if it + 1 >= numit and it + 1 < numit_max:
                # converged? (max-norm change relative to the tolerance)
                if float(np.max(np.abs(ynext - ynew) / (atol_v + rtol * np.abs(ynext)))) < conv:
                    return ynext
            ynew = ynext
        return ynew

    fixed = (dtmin == dtmax)
    if fixed:
        nstart = 1
        dt = min(dtmin, (te - t) / 2.0)
    else:
        # initial step size from the explicit derivative (Verwer: dt = min_i ytol_i / |dy_i|)
        P, L = pl_fn(t0_abs + t, y)
        nfcn += 1
        dy = P - y * L
        ytol = atol_v + rtol * np.abs(y)
        nz = dy != 0.0
        dt = te - t
        if np.any(nz):
            dt = min(dt, float(np.min(ytol[nz] / np.abs(dy[nz]))))

    while True:                                   # label 25: (re)start
        if not fixed:
            nstart += 1
            if restart:
                dt /= 10.0
            restart = True
            dt = max(dtmin, min(dt, dtmax))
            dt = _fit(t, te, dt)
            dt = min(dt, (te - t) / 2.0)
        # label 28: implicit Euler starting step
        yold = y.copy()
        ynew = iterate(y.copy(), dt, y.copy(), t + dt)
        naccpt += 1
        t += dt
        y = ynew
        if trace is not None:
            trace(t0_abs + t, y, dt, 'euler')
        dtold = dt
        ratio = 1.0
        back_to_restart = False
        nsteps = 0
        n_since_restart = 0
        while True:                               # label 60: BDF2 steps
            nsteps += 1
            if nsteps > max_steps:
                raise RuntimeError(f"twostep: more than {max_steps} steps between {t0_abs:.1f} and {t0_abs + te:.1f} s")
            c = 1.0 / ratio
            cp1 = c + 1.0
            a1 = (c + 1.0) ** 2 / (c * c + 2.0 * c)
            a2 = -1.0 / (c * c + 2.0 * c)
            dtg = dt * (1.0 + c) / (2.0 + c)
            w = a1 * y + a2 * yold
            ynew0 = np.maximum(smallc, y + ratio * (y - yold))     # extrapolated predictor
            ynew = iterate(w, dtg, ynew0, t + dt)
            if fixed:
                t += dtold
                naccpt += 1
                yold, y = y, ynew
                if dt != dtold:                    # the shortened last step was just taken
                    t = t - dtold + dt
                    break
                dt = min(dtold, te - t)
                ratio = dt / dtold
                if t >= te - eps_end:
                    break
                continue
            # error estimate (local truncation error of BDF2, max-norm, scaled)
            ytol = atol_v + rtol * np.abs(y)
            errlte = float(np.max(np.abs(c * ynew - cp1 * y + yold) / ytol))
            errlte = 2.0 * errlte / (c + c * c)
            # NEWDT
            if errlte > 1.0 and dt > dtmin:
                accept = False
                ts = t
            else:
                accept = True
                dtold = dt
                ts = t + dtold
            fac = 0.8 / math.sqrt(errlte) if errlte > 0.0 else max_growth
            growth = max_growth
            if restart_growth is not None and n_since_restart < restart_steps:
                growth = min(max_growth, restart_growth)
            dt = max(0.5, min(growth, fac)) * dt
            dt = max(dtmin, min(dt, dtmax))
            dt = _fit(ts, te, dt)
            ratio = dt / dtold
            if accept:
                failer = False
                restart = False
                t += dtold
                naccpt += 1
                n_since_restart += 1
                yold, y = y, ynew
                if trace is not None:
                    trace(t0_abs + t, y, dtold, 'bdf2')
                if t >= te - eps_end:
                    break
                continue
            if failer:                            # second failure in a row -> go back one step and restart
                nrejec += 1
                failer = False
                naccpt -= 1
                t -= dtold
                y = yold
                back_to_restart = True
                break
            nrejec += 1                           # first failure: retry the step with the smaller dt
            failer = True
        if not back_to_restart:
            break
    if stats is not None:
        stats.update(nfcn=nfcn, naccpt=naccpt, nrejec=nrejec, nstart=nstart)
    return y


# ----------------------------------------------------------------------------------------------------------
# self-test: mini SAPHIR-like photochemistry in JUICE array format vs scipy BDF (tight tolerance)
# ----------------------------------------------------------------------------------------------------------
def _mini_mechanism():
    """species, reactions (reactants, env reactants, products with coefficients, k) -- a stiff but small system"""
    species = ['O3', 'NO', 'NO2', 'OH', 'HO2', 'CO', 'HCHO', 'H2O2', 'HNO3', 'CH4', 'CH3O2', 'CH3OOH']
    M = 2.46e19
    rx = [  # (reactants[(name, coeff)], env[(name, coeff)], products[(name, coeff)], k)
        ([('NO2', 1)],              [('HV', 1)], [('NO', 1), ('O3', 1)],            8.0e-3),
        ([('NO', 1), ('O3', 1)],    [],          [('NO2', 1)],                      1.8e-14),
        ([('O3', 1)],               [('HV', 1)], [('OH', 2)],                       4.0e-6),   # lumped O1D+H2O
        ([('OH', 1), ('CO', 1)],    [],          [('HO2', 1)],                      2.4e-13),
        ([('HO2', 1), ('NO', 1)],   [],          [('OH', 1), ('NO2', 1)],           8.5e-12),
        ([('OH', 1), ('NO2', 1)],   [],          [('HNO3', 1)],                     1.0e-11),
        ([('HO2', 2)],              [],          [('H2O2', 1)],                     3.0e-12),  # self reaction
        ([('OH', 1), ('HCHO', 1)],  [],          [('HO2', 1), ('CO', 1)],           8.5e-12),
        ([('HCHO', 1)],             [('HV', 1)], [('HO2', 2), ('CO', 1)],           3.0e-5),
        ([('OH', 1), ('O3', 1)],    [],          [('HO2', 1)],                      7.0e-14),
        ([('HO2', 1), ('O3', 1)],   [],          [('OH', 1)],                       2.0e-15),
        ([('OH', 1), ('CH4', 1)],   [],          [('CH3O2', 1)],                    6.4e-15),
        ([('CH3O2', 1), ('NO', 1)], [],          [('HCHO', 1), ('HO2', 1), ('NO2', 1)], 7.7e-12),
        ([('CH3O2', 1), ('HO2', 1)], [],         [('CH3OOH', 1)],                   5.2e-12),
        ([('OH', 1), ('H2O2', 1)],  [],          [('HO2', 1)],                      1.7e-12),
        ([('O3', 1)],               [('DIL', 1)], [],                               3.0e-5),   # dilution via env reactant
        ([('CO', 1)],               [('DIL', 1)], [],                               3.0e-5),
    ]
    env_names = ['HV', 'DIL']
    sp_idx = {s: i for i, s in enumerate(species)}
    env_idx = {e: i for i, e in enumerate(env_names)}
    ode_r_rxn, ode_r_sp, ode_r_coeff = [], [], []
    env_r_rxn, env_r_name_idx, env_r_coeff = [], [], []
    st_rxn, st_sp, st_sc = [], [], []
    k = np.zeros(len(rx))
    for i, (reac, env, prod, kk) in enumerate(rx):
        k[i] = kk
        for name, c in reac:
            ode_r_rxn.append(i); ode_r_sp.append(sp_idx[name]); ode_r_coeff.append(float(c))
            st_rxn.append(i); st_sp.append(sp_idx[name]); st_sc.append(-float(c))
        for name, c in env:
            env_r_rxn.append(i); env_r_name_idx.append(env_idx[name]); env_r_coeff.append(float(c))
        for name, c in prod:
            st_rxn.append(i); st_sp.append(sp_idx[name]); st_sc.append(float(c))
    arrays = (np.array(ode_r_rxn, np.int32), np.array(ode_r_sp, np.int32), np.array(ode_r_coeff),
              np.array(env_r_rxn, np.int32), np.array(env_r_name_idx, np.int32), np.array(env_r_coeff),
              np.array(st_rxn, np.int32), np.array(st_sp, np.int32), np.array(st_sc))
    y0 = np.full(len(species), 1.0)
    for name, ppb in (('O3', 40.0), ('NO', 1.0), ('NO2', 2.0), ('CO', 200.0), ('HCHO', 3.0), ('CH4', 1900.0), ('H2O2', 0.5)):
        y0[sp_idx[name]] = ppb * 1e-9 * M
    return species, k, env_names, arrays, y0, M


def selftest(verbose=True):
    from scipy.integrate import solve_ivp
    species, k, env_names, arrays, y0, M = _mini_mechanism()
    n = len(species)
    env = np.array([1.0, 1.0])                       # HV on, DIL on
    sysm = PLSystem(n, k, env, arrays)
    # ---- two "ENZ rows": light 0..1800 s, dark 1800..3600 s (HV -> 0 = roof closure) ----
    rows = [(0.0, 1800.0, 1.0), (1800.0, 3600.0, 0.0)]
    t_ref = time.perf_counter()
    y_ref = y0.copy()
    ref_states = []
    for (ta, tb, hv) in rows:
        env[0] = hv
        sol = solve_ivp(sysm.rhs, (ta, tb), y_ref, method='BDF', rtol=1e-10, atol=1e-4, max_step=60.0)
        y_ref = np.maximum(sol.y[:, -1], SMALLC)
        ref_states.append(y_ref.copy())
    t_ref = time.perf_counter() - t_ref
    if verbose:
        print(f"reference scipy BDF (rtol 1e-10): {t_ref*1000:.0f} ms")
        print(f"{'case':28s} {'maxerr%':>8s} {'OH err%':>8s} {'HO2 err%':>9s} {'NO err%':>8s} {'steps':>6s} {'rej':>4s} {'nfcn':>6s} {'ms':>7s}")
    iOH, iHO2, iNO = species.index('OH'), species.index('HO2'), species.index('NO')
    results = {}
    for numit in (1, 2, 3, 4):
        for rtol in (1e-2, 1e-3, 1e-4):
            y = y0.copy()
            stats_tot = dict(nfcn=0, naccpt=0, nrejec=0, nstart=0)
            t0 = time.perf_counter()
            errs = []
            for (ta, tb, hv), yr in zip(rows, ref_states):
                env[0] = hv
                st = {}
                y = twostep(sysm, y, ta, tb, dtmin=1e-6, dtmax=60.0, atol=100.0, rtol=rtol, numit=numit, stats=st)
                for kk in stats_tot:
                    stats_tot[kk] += st[kk]
                mask = yr > 1e6                       # compare species above 1e6 molec cm-3 (~0.04 ppt)
                errs.append(np.abs(y[mask] / yr[mask] - 1.0) * 100.0)
            ms = (time.perf_counter() - t0) * 1000.0
            yr = ref_states[-1]
            maxerr = max(float(e.max()) for e in errs)
            rel = lambda i: abs(y[i] / yr[i] - 1.0) * 100.0
            results[(numit, rtol)] = maxerr
            if verbose:
                print(f"numit={numit} rtol={rtol:<7g}            {maxerr:8.3f} {rel(iOH):8.3f} {rel(iHO2):9.3f} {rel(iNO):8.3f} "
                      f"{stats_tot['naccpt']:6d} {stats_tot['nrejec']:4d} {stats_tot['nfcn']:6d} {ms:7.1f}")
    # fixed-step mode
    y = y0.copy()
    for (ta, tb, hv), yr in zip(rows, ref_states):
        env[0] = hv
        y = twostep(sysm, y, ta, tb, dtmin=10.0, dtmax=10.0, atol=100.0, rtol=1e-3, numit=3)
    mask = ref_states[-1] > 1e6
    if verbose:
        print(f"fixed step 10 s, numit=3: max error {float(np.max(np.abs(y[mask]/ref_states[-1][mask]-1)))*100:.3f} %")
        print("reference end state (ppb): " + ", ".join(f"{s} {ref_states[-1][i]/M*1e9:.4g}" for i, s in enumerate(species)))
    return results


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--selftest', action='store_true', help='run the mini-mechanism comparison against scipy BDF')
    a = ap.parse_args()
    if a.selftest:
        print(f"numba kernel: {NUMBA}")
        selftest()
    else:
        ap.print_help()

# scene_T1: a temperature step in an ENZ file that starts at t = 0

Test scenario by F. Rohrer (16 Nov 2026). X + X → XX with a rate constant that follows X
(k = c·X), T = 280 K for the first hour and 300 K for the second, P = 1013 hPa. At 3600 s the
air density falls by 280/300 and every concentration with it; the mixing ratios stay. The header
of `scene_T1.easy` gives the closed-form solution and the three wrong answers a model gives
when k is held per row, the density step is missed, or both.

| | X at 3600 s | XX at 3600 s | X at 7200 s | XX at 7200 s |
|---|---|---|---|---|
| exact (k follows X, density step) | 1.458318e8 | 4.270841e8 | 1.000000e8 | 4.166667e8 |
| JUICE precise | 1.458335e8 | 4.270833e8 | 0.999993e8 | 4.166670e8 |
| JUICE fast | 1.458581e8 | 4.270710e8 | 0.999857e8 | 4.166738e8 |
| JUICE twostep | 1.457647e8 | 4.271176e8 | 0.999555e8 | 4.166889e8 |
| JUICE before 17 Nov 2026 (any profile) | as above | as above | 1.0367e8 | 4.4816e8 |
| F0AM, k from the previous 10-s call | 1.456597e8 | 4.271701e8 | 0.999277e8 | 4.167028e8 |

Before 17 Nov 2026 JUICE read every ENZ line from the time of the next line on when the times
start near 0 (a 1-ns shift of the breakpoints that rounds away at SAPHIR's absolute times), so the
temperature step was never applied: the "k follows X, no density step" answer. Fixed by F. Rohrer
in MCM_SAPHIR_Tool 9e4fc04 and ported. `tests/regression/test_scene_T1.py` checks the three solver profiles
against the closed form.

With live settings off JUICE gives the same answer: a species written directly into a rate
expression is always read from the current state; only named settings are held per row. The
header expects EASY 6.00 to hold k for the row (X 3.756546e7 at 7200 s).

## F0AM

F0AM cannot write this rate constant. It computes the rate constants from the meteorology before
it integrates, a reaction has at most two reactants, and only the RO2 sum follows the state during
the integration. `f0am/run_scene_T1.m` therefore calls F0AM in a loop of short runs. Each run holds
k = c·X with X from F0AM's previous run, and the mixing ratios are passed on from run to run. When T
steps to 300 K, F0AM converts them with its own density, so the density step is F0AM's own.

```
matlab -batch "addpath(genpath('<F0AM>')); cd('validation/scene_T1/f0am'); run('run_scene_T1.m')"
```

| run length | X at 3600 s | X at 7200 s | XX at 7200 s |
|---|---|---|---|
| 3600 s (k held per row) | 4.151e7 | 3.745e7 | 4.479e8 |
| 600 s | 1.252e8 | 0.912e8 | 4.211e8 |
| 60 s | 1.447e8 | 0.9953e8 | 4.1690e8 |
| 10 s | 1.4566e8 | 0.99928e8 | 4.16703e8 |

With one run per row F0AM gives the "k held per row, density step" answer: 3.756546e7 in the
header, and 0.3 % lower in F0AM at its default tolerance. The error then shrinks in proportion to
the run length, towards the exact answer. F0AM applies the density step correctly, but it reaches
the "k follows X" solution only by this loop, which is not a set-up an F0AM user would write.

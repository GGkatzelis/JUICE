# scene_00: the clock scenario

By F. Rohrer (14 November 2026; `scene_00.easy`, `scene_00.enz`, and `scene_00.cs` = his FACSIMILE result with
EASY 7.05). It shows whether a model re-evaluates settings that depend on a species inside an output step.

- `TICK --> zeit` makes `zeit` a species that counts seconds; `T = 280 + 100*zeit/86400`; `M` follows `T`.
- `X + X --> XX` with a constant rate constant; one ENZ interval and one output step of 24 h.

Expected at 86 400 s (closed form in the header of `scene_00.easy`):

| | zeit | T / K | X | XX |
|---|---|---|---|---|
| (A) correct: the density change acts on every species | 76395.94 | 368.4212 | 0.3046754 | 0.2276621 |
| (B) settings held for the step | 86400 | 380 | 0.3665689 | 0.3167155 |

JUICE, 15 November 2026: with live settings (the default) twostep gives X = 0.3046688 and the BDF profiles
0.3050353 (0.3046752 with `atol_floor = 1e-6`: the species are of order 1 here, far below the default
absolute tolerance of 100 molecules cm-3); with `live_settings = false` every profile gives (B).
An independent integration of the three equations with scipy (LSODA, rtol 1e-11) gives (A) to 8 digits.
`tests/test_scene_00.py` checks all of this. The scene is kept out of `IDL Input`, because `run_benchmark`
runs with the settings held per step and would report (B).

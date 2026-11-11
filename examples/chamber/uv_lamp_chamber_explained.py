"""
UV-lamp chamber in code
=======================

Ethanol oxidation in a small chamber with UV-A blacklights: NO2 and H2O2 in the dark, lights on, ethanol
added, lights off. Set up step by step in the same order as the experiment window (1 Experiment,
2 Chemistry, 3 Injections, 4 Chamber, 5 Run), like saphir_isoprene_explained.py (the base case). This file
explains what is different in a lamp chamber: where the light comes from, the OH source, and the walls.

    python examples/chamber/uv_lamp_chamber_explained.py             # about 20 s
    python examples/chamber/uv_lamp_chamber_explained.py --scan      # plus three H2O2 levels (off-road A)

Copy this file and change it for your own chamber. The profile 'generic_lamps' is a starting point, not a
real chamber: put your chamber's numbers into a profile of your own (off-road D).
"""
import os
import sys

import juice              # found after "pip install -e ." in the JUICE folder (once; README, Install)

# ── where the run goes ───────────────────────────────────────────────────────────────────────────────

# Every run gets its own folder. JUICE writes into it:
#   <name>.enz                the time table: conditions per output step, lights, injection pulses
#   <name>.easy               the FACSIMILE master file: includes, lamp J values, chamber processes, injections
#   chemistry/                the MCM subset for your compounds (EASY include files)
#   <name>.experiment.json    this experiment as data (the window opens it)
#   <name>.json               the solver input (written by the prepare step)
#   <name>.cs                 the results: one column per species, molecules cm-3, every output step
#   <name>_run.log            what the solver printed
# The .enz and .easy files are ordinary FACSIMILE input: the same run works in FACSIMILE.
# docs/run_path.md follows these files through the code.
OUT_ROOT = os.path.join(juice.paths()["output_root"], "examples")    # output_root comes from the settings file


def build(name="uv_lamp_ethanol", h2o2_ppb=1000.0):
    """The experiment description. Change the values here; nothing else needs touching."""

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 1 — EXPERIMENT: which chamber, when, and where the conditions come from
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # The chamber profile sets the volume, the light and which chamber processes are on.
    # 'generic_lamps' (setups/chambers/generic_lamps.chamber.json): 10 m3 FEP bag, 350 nm blacklights.
    #
    # THE LIGHT IS SET IN THE PROFILE, NOT HERE. Its "light" block holds:
    #   jno2_full_s1       jNO2 at full lamp power, 1/s (your actinometry): the absolute scale. Here 0.002.
    #   spectrum           your lamp spectrum: a two-column file (nm, intensity) in setups/chambers/ or setups/chambers/spectra/
    #   spectrum_builtin   "blacklight_350" | "blacklight_365": synthetic stand-ins until you have a spectrum
    #   solar_like         true for a filtered xenon arc: the sun's J ratios at effective_sza
    #   j_ratios           measured J/jNO2 ratios, {"J1": 1.0e-3, ...}: they win over everything above
    # Every J value is jNO2 x the ratio that the spectrum gives (lamp_photolysis.py, cross sections and
    # quantum yields at 298 K). To change the lamps, edit the profile (chamber editor, or off-road D).
    exp = juice.new_experiment(name, chamber="generic_lamps")         # juice.chambers() lists the others
    exp.title = "Ethanol + OH under UV-A blacklights (lamp-chamber example)"

    c = exp.conditions
    c.date = "2024-07-09"              # any day: indoors, the date and place do not change the light
    c.start, c.end = "08:00", "15:00"  # the modelled window, HH:MM (UTC)
    c.step_s = 60                      # output step in seconds (one ENZ row and one result row per step)

    # Where the physical conditions come from:
    #   "manual"   the constants below (this example)
    #   "file"     a table you prepare: CSV or Excel, one row per time, columns T, p, H2O or RH, flow,
    #              lights (0-1: 0.5 = half the lamps), jNO2 (your actinometer), NO/NO2/O3 (any subset).
    #              Write an empty, documented table with
    #                  import measured_inputs; measured_inputs.write_template("my_conditions.csv")
    #              keep it anywhere (next to this script is a good place) and give its path below.
    #   ("archive" is the SAPHIR data archive and does not apply to other chambers.)
    c.source = "manual"
    c.input_file = ""                  # c.source = "file": e.g. os.path.join(os.path.dirname(__file__), "my_conditions.csv")

    # Constants (used for anything a file does not give)
    c.temperature_K = 298.0
    c.temperature_mode = "const"       # "const" = temperature_K | "enz" = the T column of the file (lamps warm the bag)
    c.pressure_hPa = 1013.0
    c.h2o_percent = 1.0                # water vapour, vol-% (dry chambers: 0.0; then add an RH column if you have one)
    c.flow_m3h = 0.5                   # flow that replaces sampled air, m3/h: dilution = flow / volume (here 5 %/h)

    # Light: when the lamps are on, as [(on, off), ...]. A 'lights' column of the file wins.
    # (c.lights_on is the same field as c.roof_open, which the SAPHIR examples use.)
    c.lights_on = [("09:00", "13:00")]
    c.j_source = "clear_sky"           # lamps: the J values always come from the profile. "measured" adds the
                                       # file's jNO2, jO1D, ... where it has them (needs c.source = "file")

    # Holding species to measurements: NO, NO2 and/or O3 set to the file's column at every step.
    c.constrain_species = []           # e.g. ["NO2"] (needs c.source = "file" with an NO2 column)

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 2 — CHEMISTRY: which mechanism, which compounds
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    exp.base_mechanism = "MCM_v3.3.1"  # or "CRACMM2", "CRACMM1_aq", "RACM2" (lumped, used whole: precursors = [])
    exp.precursors = ["C2H5OH"]        # MCM names of the compounds you inject (C2H5OH = ethanol)
    exp.addons = []                    # juice.mechanisms()["addons"]

    # Initial values: what is in the bag before the experiment starts. Expressions inside CONST();
    # M is the air density, so "20e-9*M" is 20 ppb.
    # THE OH SOURCE. A UV-A chamber has no O(1D) worth speaking of, so OH comes from a precursor:
    #   H2O2      about 1 ppm, photolysed at 254-350 nm (this example; the usual choice)
    #   HONO      tens of ppb: OH + NO, also adds NOx
    #   CH3ONO    methyl nitrite: needs an add-on with its chemistry (not in the MCM)
    exp.initials = [
        ("NO2", "20e-9*M"),                     # 20 ppb NO2
        ("H2O2", f"{h2o2_ppb:g}e-9*M"),         # the OH source
    ]

    # Rate-constant changes: (reaction as the mechanism writes it, new rate, comment); "0" removes it.
    exp.reaction_adjustments = []      # e.g. [("C2H5OH + OH --> CH3CHO + HO2", "3.0D-12*EXP(20/(T))*0.9*1.2",
                                       #         "test: main channel +20 %")]   (chemistry/*.rate shows the MCM rates)

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 3 — INJECTIONS: what goes into the chamber, how much, when
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # One Injection per event: total_ppb spread over duration_min minutes from start (HH:MM). Something that
    # is there from the start goes in exp.initials (Step 2) instead.
    # A liquid evaporated into the bag (injection_calc does the conversion with this chamber's volume):
    #     import injection_calc as ic, chambers
    #     ppb, note = ic.ppb_liquid(vol_ul=1.2, mw=46.07, density=0.789, T_K=298, P_hPa=1013,
    #                               ch=chambers.chamber("generic_lamps"))   # also ic.ppb_mfc, ic.ppb_gas_syringe
    exp.injections = [
        juice.Injection(name="ethanol", products="C2H5OH", total_ppb="200", start="10:30", duration_min=10,
                        comment="ethanol, 200 ppb over 10 min"),
        # juice.Injection(name="NO", products="NO", total_ppb="10", start="11:30"),
    ]

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 4 — CHAMBER: what the chamber adds on its own
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # generic_lamps switches off every SAPHIR wall process: they were measured in SAPHIR and are not known
    # for your chamber. Only the dilution by the flow is on. Switch on what your chamber characterisation
    # shows, with your own strengths (the expressions are EASY and use HV_in = lights, T_in, jno2, FLOW).
    ch = exp.chamber
    ch.dilution = True                 # dilution by the flow (ch.dilute_expr, from the profile's volume)
    ch.hono = False                    # HONO from the walls under light: common in FEP bags ("chamber radical source")
    ch.hcho = False                    # HCHO from the walls under light
    ch.bgreac = False                  # background OH reactivity
    ch.walldep = False                 # SAPHIR's wall deposition of HNO3, HONO, H2O2, ... (saphir.walldep.incl)
    ch.dark_hono = False               # HONO in the dark
    ch.light_source = False            # NO2 under light
    # First-order wall losses of single species, 1/s (e.g. O3 on FEP). The profile's "wall_loss_s1" sets them
    # for every experiment in this chamber; here only for this one:
    ch.wall_loss = False
    ch.wall_loss_s1 = {}               # e.g. {"O3": 1e-6, "H2O2": 2e-5}   (and ch.wall_loss = True)
    # e.g. a wall HONO source scaled to your chamber:  ch.hono = True;  ch.hono_const = "0.5"

    exp.hmax = "30"                    # EASY: the longest internal time step, s (rarely changed)
    return exp


def main():
    exp = build()
    out = os.path.join(OUT_ROOT, exp.name)

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 5 — RUN
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # The review panel's checks: errors stop the run. For a lamp chamber the notes say which J values come
    # from the lamp table and which from a file.
    for level, message in juice.check(exp, out):
        print(f"{level:8s} {message}")

    # The solver: juice.solvers() prints what each one costs.
    #   "fast"     BDF with the exact Jacobian, the default
    #   "precise"  the same at a tighter tolerance; the closest agreement with FACSIMILE
    #   "twostep"  Jacobian-free; for very large mechanisms (GECKO-A) and the complete MCM
    solver = "fast"
    kv_pd = False                      # True: also rates (.kv) and production/destruction (.pd) for Process Analysis

    # The four steps juice.run() makes, one by one, so you can stop in between:
    files = juice.write(exp, out)                                   # 1  .enz, .easy, chemistry/, .experiment.json
    #    <-- here you may edit <name>.easy or <name>.enz by hand (off-road C); the next step reads them from disk.
    #        The lamp J values are the "k[... + hv --> ...]" and "J<n> = CONST(jno2*...)" lines of the .easy.
    json_path = juice.prepare(exp, out, solver=solver, kv_pd=kv_pd, files=files)   # 2  the solver input <name>.json
    cs_path = juice.solve(json_path)                                # 3  the solver; its log is <name>_run.log
    r = juice.results(cs_path, exp)                                 # 4  the results as a table
    # (all four in one call:  run = juice.run(exp, out, solver=solver, kv_pd=kv_pd);  r = run.results())

    # Reading the results: r["OH"] molecules cm-3, r.ppb("O3"), r.hour (hours of the day), r.species,
    # r.df (pandas). The viewer: python -m juice.gui.results_viewer <cs_path>
    lit = (r.hour >= 9.0) & (r.hour < 13.0)
    print(f"\n{exp.name}: {len(r.species)} species, results in {cs_path}")
    print(f"lights on: mean OH {r['OH'][lit].mean():.2e} cm-3, NO2 {r.ppb('NO2')[lit][0]:.1f} -> {r.ppb('NO2')[lit][-1]:.1f} ppb")
    print(f"end: ethanol {r.ppb('C2H5OH')[-1]:.1f} of 200 ppb left, CH3CHO {r.ppb('CH3CHO')[-1]:.1f} ppb, O3 {r.ppb('O3')[-1]:.1f} ppb")

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # OFF-ROAD — what the window does not do in one click
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # A  Loops: any value of build() can become a loop variable, here the OH source.
    if "--scan" in sys.argv:
        for h2o2 in (200.0, 1000.0, 5000.0):
            e = build(f"uv_lamp_ethanol_H2O2_{h2o2:g}", h2o2)
            rr = juice.run(e, os.path.join(OUT_ROOT, e.name), solver=solver, kv_pd=False).results()
            print(f"H2O2 {h2o2:5g} ppb: ethanol left {rr.ppb('C2H5OH')[-1]:6.1f} ppb, max OH {rr['OH'].max():.2e} cm-3")
    # B  Your own chemistry: an EASY include file of extra reactions, in a folder of your choice,
    #        exp.include_dirs = [r"D:\my chemistry"];  exp.extra_includes = ["my_reactions.incl"]
    #    or EASY statements written verbatim after everything else (last definition wins):
    #        exp.overrides = ["k[H2O2+hv-->OH+OH] = CONST(0.0)"]   # e.g. no H2O2 photolysis: what OH is left?
    # C  Edit the input files: juice.write() leaves .easy and .enz in the run folder; change them in an editor
    #    and run juice.prepare() and juice.solve() as above. Note what you changed: the experiment description
    #    no longer tells the whole story.
    # D  Your own chamber: a profile with your volume, lamps and walls. In the window: the chamber editor on
    #    the Experiment step. In code (writes setups/chambers/my_lamps.chamber.json; the window lists it too):
    #        import copy, chambers
    #        p = copy.deepcopy(chambers.chamber("generic_lamps"))
    #        p.id, p.name, p.volume_m3 = "my_lamps", "My lamp chamber", 8.0
    #        p.light["jno2_full_s1"] = 0.0035                  # your actinometry
    #        p.light["spectrum"] = "my_lamps.csv"              # nm, intensity; in setups/chambers/ or setups/chambers/spectra/
    #        p.wall_loss_s1 = {"O3": 1e-6}                     # your measured wall losses
    #        chambers.save_profile(p)
    #        exp = juice.new_experiment("my_run", chamber="my_lamps")
    #    import lamp_photolysis as lp; print(lp.summary(lp.lamp_j(p))) shows the J values your lamps give.


if __name__ == "__main__":
    main()

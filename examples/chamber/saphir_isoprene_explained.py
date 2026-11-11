"""
SAPHIR experiment in code: the base case
========================================

Isoprene photo-oxidation under the open roof of SAPHIR, set up step by step in the same order as the
experiment window (1 Experiment, 2 Chemistry, 3 Injections, 4 Chamber, 5 Run). Every step says what
it does, what can be chosen instead, and where files go. Nothing is read from the SAPHIR wiki or
archive: every condition is typed here (or read from a file you prepare, see Step 1).

    python examples/chamber/saphir_isoprene_explained.py             # the base case, about 30 s
    python examples/chamber/saphir_isoprene_explained.py --scan      # plus the same day at three NO levels (off-road A)

Copy this file and change it: it is the template for your own experiment. The same structure appears
in every <name>_run.py that the experiment window writes next to a prepared run.
"""
import os
import sys

import juice              # found after "pip install -e ." in the JUICE folder (once; README, Install)

# ── where the run goes ───────────────────────────────────────────────────────────────────────────────

# Every run gets its own folder. JUICE writes into it:
#   <name>.enz                the time table: conditions per output step, injection pulses
#   <name>.easy               the FACSIMILE master file: includes, chamber processes, injections, initials
#   chemistry/                the MCM subset for your compounds (EASY include files)
#   <name>.experiment.json    this experiment as data (the window opens it)
#   <name>.json               the solver input (written by the prepare step)
#   <name>.cs                 the results: one column per species, molecules cm-3, every output step
#   <name>_run.log            what the solver printed
# The .enz and .easy files are ordinary FACSIMILE input: the same run works in FACSIMILE.
OUT_ROOT = os.path.join(juice.paths()["output_root"], "examples")    # output_root comes from the settings file


def build(name="saphir_isoprene", no_ppb=5.0):
    """The experiment description. Change the values here; nothing else needs touching."""

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 1 — EXPERIMENT: which chamber, when, and where the conditions come from
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # The chamber profile sets the volume, the dilution expression and which chamber processes are on.
    # juice.chambers() lists them: SAPHIR, generic_lamps (UV-lamp chamber), EUPHORE, ... or your own
    # profile made with the chamber editor (setups/chambers/<id>.chamber.json).
    exp = juice.new_experiment(name, chamber="SAPHIR")
    exp.title = "Isoprene photo-oxidation, open roof, 5 ppb NO (base case)"

    c = exp.conditions
    c.date = "2024-07-09"              # the day (UTC); with typed conditions it only sets the sun's position
    c.start, c.end = "06:00", "16:00"  # the modelled window, HH:MM UTC
    c.step_s = 60                      # output step in seconds (one ENZ row and one result row per step)

    # Where the physical conditions come from:
    #   "manual"   the constants below (this example)
    #   "file"     a table you prepare: CSV or Excel, one row per time, columns T, p, H2O or RH, flow,
    #              roof, jNO2/jO1D/..., NO/NO2/O3 (any subset; missing ones fall back to the constants).
    #              Write an empty, documented table with
    #                  import measured_inputs; measured_inputs.write_template("my_conditions.csv")
    #              keep it anywhere (next to this script is a good place) and give its path below.
    #              An ENZ file of an earlier FACSIMILE run is read the same way.
    #   "archive"  the SAPHIR data archive of that day (needs the archive drive, c.archive_root)
    c.source = "manual"
    c.input_file = ""                  # c.source = "file": e.g. os.path.join(os.path.dirname(__file__), "my_conditions.csv")

    # Constants (used for anything the file or archive does not give)
    c.temperature_K = 298.0
    c.temperature_mode = "const"       # "const" = temperature_K | "enz" = the T column of the file/archive
    c.pressure_hPa = 1013.0
    c.h2o_percent = 1.0                # water vapour, vol-%
    c.flow_m3h = 0.5                   # replenishment flow, m3/h: sets the dilution (SAPHIR: 276 m3)

    # Light. Roof open (lamp chamber: lights on) as [(open, close), ...]; a "roof" column of the file wins.
    c.roof_open = [("08:00", "15:00")]
    c.j_source = "clear_sky"           # "clear_sky" = the mechanism's clear-sky J for the sun's position
                                       # "measured"  = jNO2, jO1D, ... from the file/archive (needs c.source != "manual")

    # Holding species to measurements: NO, NO2 and/or O3 are set to the file's (or archive's) column at
    # every step instead of being computed. Needs c.source = "file" (or "archive") with those columns.
    c.constrain_species = []           # e.g. ["NO2", "O3"]

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 2 — CHEMISTRY: which mechanism, which compounds
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # The base mechanism. juice.mechanisms()["bases"] lists them:
    #   "MCM_v3.3.1"                   explicit; JUICE cuts the subset of the precursors below
    #   "CRACMM2", "CRACMM1_aq", "RACM2"   lumped, used whole: leave precursors empty
    #   an imported GECKO-A base       see docs/tutorial.md, Sect. 7
    exp.base_mechanism = "MCM_v3.3.1"
    exp.precursors = ["C5H8"]          # MCM names of the compounds you inject (C5H8 = isoprene)
    exp.addons = []                    # juice.mechanisms()["addons"], e.g. ["NOAA_BB"] (furans, phenolics)
    # Other routes: exp.inorganic_only = True (CO/NOx/O3 tests, no VOC), or exp.mechanism_stem = "<include set>"
    # for an existing EASY include set (with exp.include_dirs = [its folder]).

    # Initial values: everything starts at exp.initial_all molecules cm-3, these at their own value
    # (an expression inside CONST(); M is the air density, so "20e-9*M" is 20 ppb).
    exp.initials = []                  # e.g. [("O3", "30e-9*M"), ("CO", "150e-9*M")]

    # Rate-constant changes: the reaction exactly as the mechanism writes it, the new rate (EASY notation;
    # "0" removes the reaction), a comment. Written after the mechanism, so they replace its rate.
    exp.reaction_adjustments = []      # e.g. [("OH + NO2 --> HNO3", "1.1*KMT08", "test: +10 %")]

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 3 — INJECTIONS: what goes into the chamber, how much, when
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # One Injection per event. total_ppb is the whole amount; it is spread over duration_min minutes
    # from start (HH:MM UTC). products is the EASY product string: a species, or a split such as
    # "0.9*X + 0.1*Y". Leave products empty to record an event in the ENZ without modelling it.
    # From a liquid volume or a flow instead of ppb (injection_calc does the conversion):
    #     import injection_calc as ic, chambers
    #     ppb, note = ic.ppb_liquid(vol_ul=10, mw=68.12, density=0.681, T_K=298, P_hPa=1013,
    #                               ch=chambers.chamber("SAPHIR"))        # also ic.ppb_mfc, ic.ppb_gas_syringe
    exp.injections = [
        juice.Injection(name="NO", products="NO", total_ppb=str(no_ppb), start="07:30",
                        comment="NO from the gas bottle"),
        juice.Injection(name="isoprene", products="C5H8", total_ppb="10", start="08:30", duration_min=5,
                        comment="isoprene, 10 ppb"),
        # juice.Injection(name="O3", products="O3", total_ppb="50", start="07:45"),
        # impurities that come with an injection:   extra_channels=[("HONO+XNOy", "0.002")]
        # a delayed release through a holding species: chain=[("X_HOLD", "X", "1/600")]
    ]

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 4 — CHAMBER: what the chamber adds on its own
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # The profile picked in Step 1 set these. Each process can be switched off or changed here; the
    # rate expressions are EASY and use the ENZ inputs (HV_in = roof, T_in, jno2, FLOW ...).
    ch = exp.chamber
    ch.hono = True                     # HONO wall source under light       (strength: ch.hono_const, rate: ch.hono_rate)
    ch.hcho = True                     # HCHO wall source under light       (ch.hcho_const, ch.hcho_rate)
    ch.bgreac = True                   # background OH reactivity "Y"       (ch.y_ppb CO-equivalents -> ch.y_product)
    ch.dilution = True                 # dilution by the replenishment flow (ch.dilute_expr)
    ch.walldep = True                  # wall deposition of HNO3, HONO, H2O2, HO2, O3, N2O5 (saphir.walldep.incl)
    ch.dark_hono = True                # HONO source in the dark            (ch.dark_hono_rate)
    ch.light_source = False            # NO2 source under light             (ch.light_source_rate)
    # e.g. a stronger HONO source:  ch.hono_const = "2.6*0.7*1.5"

    # Numerical setting of the EASY file (rarely changed): the longest internal time step, s
    exp.hmax = "30"
    return exp


def main():
    exp = build()
    out = os.path.join(OUT_ROOT, exp.name)

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 5 — RUN
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # The review panel's checks: errors stop the run, warnings and notes are worth reading.
    for level, message in juice.check(exp, out):
        print(f"{level:8s} {message}")

    # The solver: juice.solvers() prints what each one costs.
    #   "fast"     BDF with the exact Jacobian, the default; agrees with FACSIMILE within 0.25 % for O3/NOx/HOx
    #   "precise"  the same at a tighter tolerance; the closest agreement with FACSIMILE
    #   "twostep"  Jacobian-free; for very large mechanisms (GECKO-A) and the complete MCM
    solver = "fast"
    kv_pd = True                       # also write rates (.kv) and production/destruction (.pd) for Process Analysis

    # The four steps juice.run() makes, one by one, so you can stop in between:
    files = juice.write(exp, out)                                   # 1  .enz, .easy, chemistry/, .experiment.json
    #    <-- here you may edit <name>.easy or <name>.enz by hand (off-road C); the next step reads them from disk
    json_path = juice.prepare(exp, out, solver=solver, kv_pd=kv_pd, files=files)   # 2  the solver input <name>.json
    cs_path = juice.solve(json_path)                                # 3  the solver; its log is <name>_run.log
    r = juice.results(cs_path, exp)                                 # 4  the results as a table
    # (all four in one call:  run = juice.run(exp, out, solver=solver, kv_pd=kv_pd);  r = run.results())

    # Reading the results: r["OH"] molecules cm-3, r.ppb("O3"), r.hour (hours of the day, UTC),
    # r.species (every column), r.df (pandas). The viewer: python -m juice.gui.results_viewer <cs_path>
    t9 = int(abs(r.hour - 9.0).argmin())
    print(f"\n{exp.name}: {len(r.species)} species, results in {cs_path}")
    print(f"09:00 UTC: O3 {r.ppb('O3')[t9]:.1f} ppb, isoprene {r.ppb('C5H8')[t9]:.2f} ppb, OH {r['OH'][t9]:.2e} cm-3")
    print(f"end:       O3 {r.ppb('O3')[-1]:.1f} ppb, NO2 {r.ppb('NO2')[-1]:.2f} ppb")

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # OFF-ROAD — what the window does not do in one click
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # A  Loops: any value of build() can become a loop variable.
    if "--scan" in sys.argv:
        for no in (1.0, 5.0, 20.0):
            e = build(f"saphir_isoprene_NO{no:g}", no)
            rr = juice.run(e, os.path.join(OUT_ROOT, e.name), solver=solver, kv_pd=False).results()
            print(f"NO {no:4g} ppb: O3 at the end {rr.ppb('O3')[-1]:6.1f} ppb, max OH {rr['OH'].max():.2e} cm-3")
    # B  Your own chemistry: an EASY include file of extra reactions, in a folder of your choice,
    #        exp.include_dirs = [r"D:\my chemistry"];  exp.extra_includes = ["my_reactions.incl"]
    #    or EASY statements written verbatim after everything else (last definition wins):
    #        exp.overrides = ["WALLDEP = CONST(0.5)"]          # e.g. half the wall deposition
    # C  Edit the input files: juice.write() leaves .easy and .enz in the run folder; change them in an
    #    editor and run juice.prepare() and juice.solve() as above. The experiment description is then no
    #    longer the whole story; note what you changed.
    # D  Another chamber or design: juice.new_experiment(name, chamber="generic_lamps") for a lamp chamber,
    #    or change exp.chamber.dilute_expr, the injections and the conditions freely; juice.check() tells
    #    you when something does not fit together.


if __name__ == "__main__":
    main()

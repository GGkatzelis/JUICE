"""
Field observations in code: a constrained site
==============================================

One summer day at a Los Angeles-like site: NO, NO2, O3, CO, isoprene and toluene are held to the measured
values, everything else (OH, HO2, RO2, HCHO, PAN, ...) is computed, and the modelled HCHO is compared with
the measured HCHO. Set up step by step in the same order as the field window (1 Field case,
2 Chemistry, 3 Emissions and constraints, 4 Mixing and deposition, 5 Run), like the chamber examples.

    python examples/field/field_site.py            # about 1 min (48 h: a spin-up day and the measured day)
    python examples/field/field_site.py --scan     # plus two other generic lifetimes (off-road A)

The numbers of the table are INVENTED (examples/field/site_day/site_day.csv). Replace the table with your
measurements and this file is your case.
"""
import os
import sys

import juice              # found after "pip install -e ." in the JUICE folder (once; README, Installation)

# ── where the run goes ───────────────────────────────────────────────────────────────────────────────

# Every run gets its own folder. JUICE writes into it:
#   <name>.enz                the table on the model's time grid (conditions, J values, held species, ...)
#   <name>.easy               the FACSIMILE master file: chemistry, emissions, held species, dilution, deposition
#   chemistry/                the MCM subset for your compounds (EASY include files)
#   <name>.field.json         this case as data (the field window opens it)
#   <day>.observations.csv    the C_ and O_ columns, which the results viewer plots next to the model
#   <name>.json               the solver input (written by the prepare step)
#   <name>.cs                 the results: one column per species, molecules cm-3, every output step
#   <name>_skill.csv          the model against the O_ columns: bias, NMB, RMSE, R2 (Step 5)
#   <name>_run.log            what the solver printed
# The .enz and .easy files are ordinary FACSIMILE input: the same run works in FACSIMILE.
OUT_ROOT = os.path.join(juice.paths()["output_root"], "examples")    # output_root comes from the settings file

# The input table. Keep yours anywhere; next to this script is a good place.
TABLE = os.path.join(juice.ROOT, "examples", "field", "site_day", "site_day.csv")


def build(name="site_day", lifetime_h=6.0):
    """The case description. Change the values here; nothing else needs touching."""

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 1 — FIELD CASE: which kind of run, which table
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # The kind of run decides what the table must hold and how the air is mixed:
    #   "constrained"  a site: measured species held, the rest computed (this example)
    #   "plume"        a plume or transect: air diluted towards the background (examples/field/plume_ageing)
    #   "lagrangian"   an air parcel along a trajectory, with emissions (examples/field/la_parcel)
    #
    # THE TABLE: one row per time (UTC), comma, semicolon or tab separated, lines starting with # ignored.
    #   time        "YYYY-MM-DD HH:MM", or "HH:MM" / hours of day together with case.date below
    #   lat, lon    degrees, per row (a trajectory) or constant (a site; or type them in the case below)
    #   T, p, RH    K, hPa, % (or H2O in vol-%); degC and Pa are recognised
    #   h           boundary-layer height, m (needed for flux emissions and deposition velocities)
    #   jNO2, ...   measured photolysis frequencies, 1/s; without them the clear-sky values of the place and time
    #   E_<sp> [unit]   emission         B_<sp> [unit]   background the air mixes with
    #   C_<sp> [unit]   held to this     O_<sp> [unit]   observed only: compared, not used
    # Species names are those of the mechanism (MCM: C5H8 isoprene, TOLUENE, C2H5OH ethanol, ...).
    # An empty, documented table with the columns of a kind of run:
    #     import field_spec; field_spec.write_template("my_site.csv", "constrained")
    case = juice.new_field_case(name, run_type="constrained", table=TABLE)
    case.description = "Synthetic diurnal cycle; NO, NO2, O3, CO, isoprene and toluene held; HCHO observed for comparison."
    case.site = "Pasadena-like site (synthetic)"
    case.date = ""                     # only needed when the table's times have no date
    case.lat = None                    # degrees; only when the table has no lat column
    case.lon = None
    case.step_s = 60                   # output step in seconds; the table may be coarser (it is interpolated)

    # A constrained site starts from nothing: repeat the first day before the output, so that the species
    # that are not held (PAN, HCHO, ...) reach their daily cycle. 1 day is usually enough, 2 for long-lived ones.
    case.spin_up_days = 1.0
    # Everything that is neither held nor deposited is lost at this first-order rate: it stands for the
    # mixing and deposition the box does not resolve. The result depends on it; vary it (off-road A).
    case.generic_lifetime_h = lifetime_h

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 2 — CHEMISTRY: which mechanism, which compounds
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # The same choices as for a chamber experiment, as one dictionary:
    #   base_mechanism   "MCM_v3.3.1" (the subset of the precursors) | "CRACMM2" | "CRACMM1_aq" | "RACM2" (whole)
    #   precursors       MCM names of the VOCs in the table (held, emitted or in the background)
    #   addons           e.g. ["NOAA_BB"] for a fire plume (furans, phenolics)
    #   inorganic_only, mechanism_stem, include_dirs, reaction_adjustments: as in examples/chamber/saphir_isoprene.py
    case.chemistry = {
        "base_mechanism": "MCM_v3.3.1",
        "addons": [],
        "precursors": ["C5H8", "TOLUENE"],
        "mechanism_stem": "",
        "inorganic_only": False,
        "include_dirs": [],
        "reaction_adjustments": [],
    }

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 3 — EMISSIONS AND CONSTRAINTS
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # Both come from the table: E_ columns are emissions, C_ columns are held. Here: C_NO, C_NO2, C_O3, C_CO,
    # C_C5H8, C_TOLUENE (held) and O_HCHO (compared only). To stop holding a species, rename its column
    # C_ -> O_; it is then computed and compared.
    case.default_emission_unit = "molec/cm2/s"   # for E_ columns without [unit]: ppb/s, molec/cm3/s, molec/cm2/s,
                                                 # mol/km2/s, mol/km2/h, mol/s, molec/s
    case.area_km2 = 0.0                # footprint A for the mol/s and molec/s units (emission / (h * A))
    case.emission_scale = {}           # factors on the E_ columns: {"NOX": 0.5, "VOC": 2, "C5H8": 1.2, "ALL": 1}

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 4 — MIXING AND DEPOSITION
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # Dilution of the box towards the background air (the B_ columns of the table):
    #   "none"      no dilution (a constrained site: the generic lifetime above does that job)
    #   "pbl"       from the growth of the boundary layer, column h (case.k_collapse_per_h while it does not grow)
    #   "table"     a kdil column, 1/s
    #   "constant"  case.k_dil_per_h
    #   "tracer"    from the decay of a held tracer (case.tracer, a species with a C_ column, e.g. "CO")
    case.dilution = "none"
    case.k_dil_per_h = 0.0
    case.k_collapse_per_h = 0.3
    case.tracer = ""
    # Deposition: a velocity (loss v/h, needs column h) or a lifetime per species.
    case.deposition_velocity_cm_s = {}           # e.g. {"HNO3": 2.0, "H2O2": 1.0}
    case.deposition_lifetime_h = {}              # e.g. {"PAN": 24.0}
    # Start values of species that are not held, ppb (everything else starts at 0).
    case.initials = {}                           # e.g. {"PAN": 0.5}
    return case


def main():
    case = build()
    out = os.path.join(OUT_ROOT, case.name)

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # STEP 5 — RUN
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # The checks read the table: unknown columns, units, gaps, species the mechanism does not have.
    for level, message in juice.check(case):
        if level != "info":
            print(f"{level:8s} {message}")

    # The solver: juice.solvers() prints what each one costs.
    #   "fast"     BDF with the exact Jacobian, the default
    #   "precise"  the same at a tighter tolerance; the closest agreement with FACSIMILE
    #   "twostep"  Jacobian-free; for very large mechanisms (GECKO-A) and the complete MCM
    solver = "fast"
    kv_pd = False                      # True: also rates (.kv) and production/destruction (.pd) for Process Analysis

    # The four steps juice.run() makes, one by one, so you can stop in between:
    files = juice.write(case, out)                                  # 1  .enz, .easy, chemistry/, .field.json
    #    <-- here you may edit <name>.easy or <name>.enz by hand (off-road C); the next step reads them from disk
    json_path = juice.prepare(case, out, solver=solver, kv_pd=kv_pd, files=files)  # 2  the solver input <name>.json
    cs_path = juice.solve(json_path)                                # 3  the solver; its log is <name>_run.log
    r = juice.results(cs_path, case)                                # 4  the results as a table
    rows, skill_path = juice.skill(case, cs_path)                   # 5  the model against the O_ columns
    # (all in one call:  run = juice.run(case, out, solver=solver, kv_pd=kv_pd);  r = run.results();
    #  run.skill_path is the skill table)

    # Reading the results: r["OH"] molecules cm-3, r.ppb("HCHO") with the model's air density, r.hour (hours
    # from midnight UTC of the first output time: the spin-up day comes first), r.species, r.df (pandas).
    # The viewer, with the observations: python -m juice.gui.results_viewer <cs_path>
    measured_day = r.hour >= r.hour[0] + 24.0 * case.spin_up_days
    print(f"\n{case.name}: {len(r.species)} species, results in {cs_path}")
    print(f"measured day: mean OH {r['OH'][measured_day].mean():.2e} cm-3, max O3 {r.ppb('O3')[measured_day].max():.0f} ppb (held)")
    for row in rows:
        if row.get("n", 0) >= 2:
            print(f"{row['species']}: observed {row['obs_mean']:.2f} ppb, modelled {row['model_mean']:.2f} ppb, "
                  f"NMB {row['nmb']:+.0%}, R2 {row['r2']:.2f}   ({skill_path})")

    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # OFF-ROAD — what the window does not do in one click
    # ════════════════════════════════════════════════════════════════════════════════════════════════
    # A  Loops: the generic lifetime sets how much HCHO the box keeps. Two more values:
    if "--scan" in sys.argv:
        for hours in (3.0, 24.0):
            c = build(f"site_day_tau{hours:g}h", hours)
            run = juice.run(c, os.path.join(OUT_ROOT, c.name), solver=solver, kv_pd=False)
            hcho = next(x for x in juice.skill(c, run.cs)[0] if x["species"] == "HCHO")
            print(f"lifetime {hours:4g} h: HCHO modelled {hcho['model_mean']:.2f} ppb "
                  f"(observed {hcho['obs_mean']:.2f}), NMB {hcho['nmb']:+.0%}")
    # B  Your own chemistry: case.chemistry["include_dirs"] = [r"D:\my chemistry"] with your include files,
    #    or case.chemistry["reaction_adjustments"] = [("OH + NO2 --> HNO3", "1.1*KMT08", "test: +10 %")]
    # C  Edit the input files: juice.write() leaves .easy and .enz in the run folder; change them in an editor
    #    and run juice.prepare() and juice.solve() as above. Note what you changed.
    # D  Another kind of run: start from the other examples, e.g.
    #        case = juice.load(os.path.join(juice.ROOT, "examples", "field", "la_parcel", "la_parcel.field.json"))
    #    or a saved case of your own (<name>.field.json, which the field window also opens).


if __name__ == "__main__":
    main()

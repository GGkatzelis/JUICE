"""
Isoprene photo-oxidation in SAPHIR: open roof, 5 ppb NO, 10 ppb isoprene.

    python examples/chamber/saphir_isoprene.py          # about 30 s, then a figure

Only what this experiment sets is written here; everything else keeps the defaults of the SAPHIR chamber
profile. saphir_isoprene_explained.py is the same experiment step by step, with what every setting means
and what else can be chosen. docs/juice_from_code.md is the reference.
"""
import juice

exp = juice.new_experiment("saphir_isoprene", chamber="SAPHIR")
exp.title = "Isoprene photo-oxidation, open roof, 5 ppb NO (base case)"

# When, and under which conditions (typed here; for a measured table: c.source = "file")
c = exp.conditions
c.date, c.start, c.end = "2024-07-09", "06:00", "16:00"        # UTC
c.temperature_K, c.temperature_mode = 298.0, "const"
c.pressure_hPa, c.h2o_percent, c.flow_m3h = 1013.0, 1.0, 0.5   # hPa, vol-% water, m3/h replenishment flow
c.roof_open = [("08:00", "15:00")]
c.j_source = "clear_sky"                                        # J values from the sun's position (none measured)

# Chemistry: the MCM subset of isoprene
exp.base_mechanism = "MCM_v3.3.1"
exp.precursors = ["C5H8"]

# Injections: what, how much in total, when (and over how many minutes)
exp.injections = [
    juice.Injection(name="NO", products="NO", total_ppb="5.0", start="07:30", comment="NO from the gas bottle"),
    juice.Injection(name="isoprene", products="C5H8", total_ppb="10", start="08:30", duration_min=5,
                    comment="isoprene, 10 ppb"),
]

if __name__ == "__main__":
    run = juice.run(exp)                    # input files, solver, results: in runs/saphir_isoprene
    print(run.summary())
    run.results().plot(["O3", "NO", "NO2", "C5H8", "OH", "HO2"])

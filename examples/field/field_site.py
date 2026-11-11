"""
A constrained field site: NO, NO2, O3, CO, isoprene and toluene held to measurements, the rest computed.

    python examples/field/field_site.py          # about 1 min (a spin-up day and the measured day), then a figure

The table's numbers are invented (examples/field/site_day/site_day.csv); put your own measurements in a table
with the same columns and point TABLE at it. Only what this case sets is written here; field_site_explained.py
is the same case step by step, with the table's columns and every choice of the field window.
"""
import os

import juice

TABLE = os.path.join(juice.ROOT, "examples", "field", "site_day", "site_day.csv")

case = juice.new_field_case("site_day", run_type="constrained", table=TABLE)
case.description = "Synthetic diurnal cycle; NO, NO2, O3, CO, isoprene and toluene held; HCHO observed for comparison."
case.site = "Pasadena-like site (synthetic)"
case.spin_up_days = 1.0              # repeat the first day once, so that PAN, HCHO, ... reach their daily cycle
case.generic_lifetime_h = 6.0        # first-order loss of everything not held: stands for mixing and deposition
case.chemistry = {"base_mechanism": "MCM_v3.3.1", "precursors": ["C5H8", "TOLUENE"]}   # MCM names of the VOCs in the table

if __name__ == "__main__":
    run = juice.run(case, kv_pd=False)      # input files, solver, results: in runs/site_day
    print(run.summary())
    print("model against the observed (O_) columns:", run.skill_path)
    run.results().plot(["OH", "HO2", "HCHO", "PAN"])

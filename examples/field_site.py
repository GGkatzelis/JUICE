"""
A constrained field site from code: the site_day example, run twice with two generic lifetimes.

    python examples/field_site.py           # ~4 min (two 48-h runs: spin-up day + measured day)

The example's numbers are invented (Los Angeles-like, August).  NO, NO2, O3, CO, isoprene and toluene
are held to the table; HCHO is only compared.  The comparison (bias, NMB, RMSE, R2) is written by the
solver as <name>_skill.csv next to the results.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import juice
import field_spec

OUT = os.path.join(juice.paths()["output_root"], "examples")
path, _about = field_spec.list_examples()["site_day"]

for hours in (6.0, 24.0):
    case = juice.load(path)                            # the .field.json of the example; its table is found next to it
    case.name = f"site_day_tau{hours:g}h"
    case.generic_lifetime_h = hours                    # first-order loss of everything that is neither held nor deposited
    run = juice.run(case, os.path.join(OUT, case.name), solver="fast", kv_pd=False)
    r = run.results()
    day = r.hour >= 24.0                               # hour 0-24 is the spin-up day
    print(f"lifetime {hours:4g} h: mean modelled HCHO on the measured day {r.ppb('HCHO')[day].mean():.2f} ppb  "
          f"(skill table: {os.path.join(run.out_dir, case.name + '_skill.csv')})")

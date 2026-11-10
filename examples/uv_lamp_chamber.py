"""
A UV-lamp chamber: NO2 in the dark, lights on, ethanol, lights off.

    python examples/uv_lamp_chamber.py      # ~20 s

The chamber profile 'generic_lamps' (chambers/generic_lamps.chamber.json) is a 10 m3 chamber with
350 nm blacklights; its J values come from the lamp spectrum.  The OH source is 1 ppm H2O2 at the start,
as in most UV-chamber experiments.  For a lamp chamber, roof_open means lights on.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import juice

OUT = os.path.join(juice.paths()["output_root"], "examples")

exp = juice.new_experiment("uv_lamp_ethanol", chamber="generic_lamps")   # juice.chambers() lists the others
c = exp.conditions
c.source, c.date, c.start, c.end, c.step_s = "manual", "2024-07-09", "08:00", "15:00", 60
c.temperature_K, c.temperature_mode, c.pressure_hPa, c.h2o_percent, c.flow_m3h = 298.0, "const", 1013.0, 1.0, 0.5
c.roof_open = [("09:00", "13:00")]                 # lights on 09:00, off 13:00

exp.precursors = ["C2H5OH"]
exp.initials = [("NO2", "20e-9*M"), ("H2O2", "1000e-9*M")]   # expression inside CONST(); M = air density
exp.injections = [juice.Injection(name="ethanol", products="C2H5OH", total_ppb="200", start="10:30", duration_min=10)]

for level, message in juice.check(exp):
    print(level, message)

run = juice.run(exp, os.path.join(OUT, exp.name), solver="fast", kv_pd=False)
print(run.summary())
r = run.results()
print(f"ethanol left at the end: {r.ppb('C2H5OH')[-1]:.1f} of 200 ppb")

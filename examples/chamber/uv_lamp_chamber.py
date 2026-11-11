"""
Ethanol + OH in a small chamber with UV-A blacklights: NO2 and H2O2 in the dark, lights on, ethanol, lights off.

    python examples/chamber/uv_lamp_chamber.py          # about 20 s, then a figure

Only what this experiment sets is written here; everything else keeps the defaults of the chamber profile
'generic_lamps' (setups/chambers/generic_lamps.chamber.json: a 10 m3 bag, 350 nm blacklights). The lamps
(jNO2 at full power, spectrum) are set in that profile, not here. uv_lamp_chamber_explained.py is the same
experiment step by step: the light, the OH source, the walls, and a profile for your own chamber.
"""
import juice

exp = juice.new_experiment("uv_lamp_ethanol", chamber="generic_lamps")
exp.title = "Ethanol + OH under UV-A blacklights (lamp-chamber example)"

c = exp.conditions
c.date, c.start, c.end = "2024-07-09", "08:00", "15:00"        # UTC; indoors the date does not change the light
c.temperature_K, c.temperature_mode = 298.0, "const"
c.pressure_hPa, c.h2o_percent, c.flow_m3h = 1013.0, 1.0, 0.5   # hPa, vol-% water, m3/h replenishment flow
c.lights_on = [("09:00", "13:00")]

# Chemistry: the MCM subset of ethanol. A UV-A chamber needs an OH source: 1 ppm H2O2 from the start.
exp.base_mechanism = "MCM_v3.3.1"
exp.precursors = ["C2H5OH"]
exp.initials = [("NO2", "20e-9*M"), ("H2O2", "1000e-9*M")]     # M is the air density: 20 ppb NO2, 1000 ppb H2O2

exp.injections = [
    juice.Injection(name="ethanol", products="C2H5OH", total_ppb="200", start="10:30", duration_min=10,
                    comment="ethanol, 200 ppb over 10 min"),
]

if __name__ == "__main__":
    run = juice.run(exp)                    # input files, solver, results: in runs/uv_lamp_ethanol
    print(run.summary())
    run.results().plot(["NO2", "O3", "H2O2", "C2H5OH", "CH3CHO", "OH"])

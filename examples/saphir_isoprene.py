"""
SAPHIR, isoprene under the open roof, built in code - and the same day at three NO levels.

    python examples/saphir_isoprene.py            # one run (~30 s)
    python examples/saphir_isoprene.py --scan     # plus NO = 1, 5 and 20 ppb

Invented conditions (manual, clear-sky J): nothing here is read from the SAPHIR archive.
Results go to <output_root>/examples/ (output_root from the settings file).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import juice

OUT = os.path.join(juice.paths()["output_root"], "examples")


def isoprene_day(name="saphir_isoprene", no_ppb=5.0):
    exp = juice.new_experiment(name, chamber="SAPHIR")

    # when and under which conditions: constants, J values from the clear-sky formulas of the mechanism
    c = exp.conditions
    c.date, c.start, c.end, c.step_s = "2024-07-09", "06:00", "16:00", 60
    c.source = "manual"                   # manual | archive | file
    c.temperature_K, c.temperature_mode = 298.0, "const"
    c.pressure_hPa, c.h2o_percent, c.flow_m3h = 1013.0, 1.0, 0.5
    c.roof_open = [("08:00", "15:00")]
    c.j_source = "clear_sky"

    # chemistry: the MCM v3.3.1 subset of isoprene
    exp.precursors = ["C5H8"]

    # injections: total amount, start time, duration
    exp.injections = [
        juice.Injection(name="NO", products="NO", total_ppb=str(no_ppb), start="07:30"),
        juice.Injection(name="isoprene", products="C5H8", total_ppb="10", start="08:30", duration_min=5),
    ]
    # the SAPHIR chamber effects (HONO and HCHO wall sources, dilution, ...) keep their defaults;
    # switch one off with e.g.  exp.chamber.hono = False
    return exp


if __name__ == "__main__":
    exp = isoprene_day()
    run = juice.run(exp, os.path.join(OUT, exp.name), solver="fast", kv_pd=False)
    print(run.summary())

    r = run.results()
    t9 = int(abs(r.hour - 9.0).argmin())         # half an hour after the isoprene injection
    print(f"at 09:00 UTC: O3 {r.ppb('O3')[t9]:.1f} ppb, isoprene {r.ppb('C5H8')[t9]:.2f} ppb, OH {r['OH'][t9]:.2e} cm-3")

    if "--scan" in sys.argv:              # what a GUI cannot do in one click: a loop over anything
        for no in (1.0, 5.0, 20.0):
            e = isoprene_day(f"saphir_isoprene_NO{no:g}", no)
            rr = juice.run(e, os.path.join(OUT, e.name), solver="fast", kv_pd=False).results()
            print(f"NO {no:4g} ppb: O3 at the end {rr.ppb('O3')[-1]:6.1f} ppb, max OH {rr['OH'].max():.2e} cm-3")

"""
make_examples.py — the three synthetic field cases in this folder, written by this script.

    python field_examples/make_examples.py          # rewrites site_day/, plume_ageing/, la_parcel/

Every number is invented but plausible for the Los Angeles basin in August; nothing here is a
measurement.  The point of the examples is the format: each case is one <name>.field.json and one
<name>.csv whose comment lines explain every column, so a user can copy a case, replace the numbers
with their own and keep the structure.  The field window loads them with "Load an example"; the
copies land in the user's project folder, these originals stay as they are.

    site_day      constrained site: a 24-h diurnal cycle at a Pasadena-like site, NO, NO2, O3, CO,
                  isoprene and toluene held to the "measurements"; HCHO deliberately not held, to
                  compare the model with it afterwards
    plume_ageing  plume or transect: 3 h downwind of a source, the dilution rate from the decay of
                  CO, background air mixed in, initial values of the fresh plume
    la_parcel     Lagrangian air parcel: 4 h from Santa Monica Bay to Pasadena on the sea breeze,
                  emissions in mol/s over an 8 x 8 km footprint, boundary-layer growth as dilution,
                  background air, deposition of O3 and HNO3
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import field_spec as fs                                  # noqa: E402
from enz_builder import compute_solar_zenith             # noqa: E402

DAY = "2021-08-15"
ENZ_EPOCH_UNIX = 946684800.0


def utc_seconds(hhmm: str, day_offset: int = 0) -> float:
    """ENZ seconds of DAY hh:mm (+ day_offset days)."""
    import datetime as dt
    d = dt.datetime.strptime(DAY, "%Y-%m-%d").replace(tzinfo=dt.timezone.utc) + dt.timedelta(days=day_offset)
    h, m = (int(x) for x in hhmm.split(":"))
    return d.timestamp() - ENZ_EPOCH_UNIX + 3600 * h + 60 * m


def stamp(sec: float) -> str:
    import datetime as dt
    return dt.datetime.fromtimestamp(sec + ENZ_EPOCH_UNIX, dt.timezone.utc).strftime("%Y-%m-%d %H:%M")


def jno2_clear(sec: np.ndarray, lat, lon, cloud: float = 1.0) -> np.ndarray:
    """jNO2 with the clear-sky shape JUICE uses (1.165e-2 cos^0.244 exp(-0.267/cos)), times a cloud factor."""
    lat = np.broadcast_to(np.asarray(lat, float), sec.shape); lon = np.broadcast_to(np.asarray(lon, float), sec.shape)
    sza = np.array([compute_solar_zenith([s], la, lo)[0] for s, la, lo in zip(sec, lat, lon)])
    cz = np.cos(np.radians(np.minimum(sza, 89.9)))
    return np.where(sza < 89.9, 1.165e-2 * cz ** 0.244 * np.exp(-0.267 / cz), 0.0) * cloud


def bell(x: np.ndarray, centre: float, width: float) -> np.ndarray:
    return np.exp(-0.5 * ((x - centre) / width) ** 2)


def write_table(path: str, header: str, columns: dict, fmt: dict | None = None) -> None:
    fmt = fmt or {}
    names = list(columns)
    n = len(columns[names[0]])
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(header.rstrip("\n") + "\n")
        f.write(",".join(names) + "\n")
        for i in range(n):
            cells = []
            for k in names:
                v = columns[k][i]
                cells.append(v if isinstance(v, str) else format(v, fmt.get(k, ".4g")))
            f.write(",".join(cells) + "\n")


COMMON = """# {title}
# A SYNTHETIC JUICE field case: every number is invented (plausible for Los Angeles in August), nothing
# was measured.  Written by field_examples/make_examples.py; copy it and replace the numbers with yours.
#
# One row per time, UTC.  Lines starting with # are ignored.  Comma, semicolon or tab separated; a
# decimal comma is understood in semicolon tables.  Column names are matched ignoring case.
#   time      'YYYY-MM-DD HH:MM' (UTC); also 'HH:MM' or hours of day together with the date of the
#             case, unix seconds, or seconds since 2000-01-01
#   lat, lon  position in degrees, per row (a trajectory) or constant (a site); may also be typed
#             once in the case instead of a column
#   T         air temperature, K (degC is recognised when below 150)
#   p         pressure, hPa (Pa is recognised when above 5000)
#   RH        relative humidity, %  (or H2O in vol-% instead)
#   h         boundary-layer height, m
#   jNO2      measured NO2 photolysis frequency, 1/s: the clear-sky frequencies of the mechanism are
#             scaled by it; without the column the clear-sky values for time and position are used.
#             jO1D, jHONO, jH2O2, jHCHO_r, jHCHO_m are read the same way.
#   E_<species> [unit]   emission (units: ppb/s, molec/cm3/s, molec/cm2/s, mol/km2/s, mol/km2/h,
#                        mol/s, molec/s; the last two need the footprint area of the case)
#   B_<species> [unit]   background mixing ratio the air is diluted towards (ppb, ppt, ppm, molec/cm3)
#   C_<species> [unit]   measured mixing ratio the species is held to (same units)
#   Species names are those of the mechanism (MCM: C5H8 isoprene, TOLUENE, C2H5OH ethanol, ...).
#"""


# ── 1. constrained site ──────────────────────────────────────────────────────────────────────────
def site_day(folder: str) -> str:
    lat, lon = 34.1403, -118.1254                       # the Caltech site of SUNVEx
    sec = np.array([utc_seconds("07:00") + 1800 * i for i in range(49)])   # 00:00-24:00 PDT, every 30 min
    local_h = (np.arange(49) * 0.5) % 24
    T = 291 + 13 * bell(local_h, 15, 4)
    RH = 75 - 40 * bell(local_h, 15, 4)
    h = 350 + 900 * bell(local_h, 15, 3)
    j = jno2_clear(sec, lat, lon, cloud=0.92)
    no = 0.4 + 7 * bell(local_h, 7.5, 1.2) + 1.5 * bell(local_h, 18, 2)
    no2 = 6 + 12 * bell(local_h, 7.5, 1.5) + 6 * bell(local_h, 19, 2.5) - 3 * bell(local_h, 14, 2)
    o3 = 12 + 60 * bell(local_h, 14.5, 3)
    co = 170 + 120 * bell(local_h, 7.5, 1.5) + 40 * bell(local_h, 18, 2)
    iso = 0.1 + 2.2 * bell(local_h, 14, 3)
    tol = 0.6 + 1.2 * bell(local_h, 7.5, 1.5) + 0.5 * bell(local_h, 18, 2)
    header = COMMON.format(title="Constrained site: one diurnal cycle at a Pasadena-like site") + """
# The measured species are HELD: at every output step the model takes the value of the table (interpolated);
# radicals and everything not listed run free with the lifetime set in the case (24 h) and the first day
# is repeated once as spin-up (case: spin_up_days = 1).  HCHO is deliberately NOT held: compare the modelled
# HCHO with your measurement afterwards.  Columns h is informative here (no emissions, no dilution).
#"""
    cols = {"time": [stamp(s) for s in sec], "lat": np.full(49, lat), "lon": np.full(49, lon), "T": T, "p": np.full(49, 985.0),
            "RH": RH, "h": h, "jNO2": j, "C_NO [ppb]": no, "C_NO2 [ppb]": no2, "C_O3 [ppb]": o3, "C_CO [ppb]": co,
            "C_C5H8 [ppb]": iso, "C_TOLUENE [ppb]": tol}
    os.makedirs(folder, exist_ok=True)
    table = os.path.join(folder, "site_day.csv")
    write_table(table, header, cols, {"jNO2": ".3e", "T": ".1f", "RH": ".1f", "h": ".0f"})
    spec = fs.FieldSpec(name="site_day", run_type="constrained", table=table, site="Pasadena-like site (synthetic)",
                        description="Synthetic diurnal cycle; NO, NO2, O3, CO, isoprene and toluene held; HCHO free for comparison.",
                        chemistry={"base_mechanism": "MCM_v3.3.1", "addons": [], "precursors": ["C5H8", "TOLUENE"], "mechanism_stem": "",
                                   "inorganic_only": False, "include_dirs": [], "reaction_adjustments": []},
                        step_s=60, dilution="none", spin_up_days=1.0, generic_lifetime_h=24.0)
    return fs.save_spec(spec, os.path.join(folder, "site_day"))


# ── 2. plume or transect ─────────────────────────────────────────────────────────────────────────
def plume_ageing(folder: str) -> str:
    n = 19                                                # every 10 min, 3 h, 13:00-16:00 PDT
    sec = np.array([utc_seconds("20:00") + 600 * i for i in range(n)])
    frac = np.arange(n) / (n - 1)
    lat = 34.05 + 0.12 * frac; lon = -118.25 + 0.30 * frac
    T = 303 - 2 * frac; RH = 32 + 6 * frac
    j = jno2_clear(sec, lat, lon, cloud=1.0)
    co = 120 + 480 * np.exp(-frac * 3 * 0.36)             # the tracer: decays towards its background of 120 ppb
    header = COMMON.format(title="Plume or transect: 3 h of ageing downwind of a source") + """
# The DILUTION RATE comes from the decay of the measured tracer CO above its background (case: dilution =
# 'tracer', tracer = 'CO'); kdil could also be a column (dilution = 'table') or one constant.  B_ columns
# are the air the plume mixes with; the fresh plume's composition is given as initial values in the case.
#"""
    cols = {"time": [stamp(s) for s in sec], "lat": lat, "lon": lon, "T": T, "p": np.full(n, 1000.0), "RH": RH, "jNO2": j,
            "C_CO [ppb]": co, "B_CO [ppb]": np.full(n, 120.0), "B_O3 [ppb]": np.full(n, 45.0), "B_NO [ppb]": np.full(n, 0.1),
            "B_NO2 [ppb]": np.full(n, 0.6), "B_HCHO [ppb]": np.full(n, 1.0)}
    os.makedirs(folder, exist_ok=True)
    table = os.path.join(folder, "plume_ageing.csv")
    write_table(table, header, cols, {"jNO2": ".3e", "T": ".1f", "RH": ".1f", "lat": ".4f", "lon": ".4f"})
    spec = fs.FieldSpec(name="plume_ageing", run_type="plume", table=table, site="Los Angeles basin (synthetic)",
                        description="Synthetic plume: NOx and VOC rich air diluted towards background at the rate the CO decay gives.",
                        chemistry={"base_mechanism": "MCM_v3.3.1", "addons": [], "precursors": ["C2H4", "C3H6", "TOLUENE"], "mechanism_stem": "",
                                   "inorganic_only": False, "include_dirs": [], "reaction_adjustments": []},
                        step_s=60, dilution="tracer", tracer="CO",
                        initials={"NO": 20.0, "NO2": 30.0, "O3": 20.0, "CO": 600.0, "HCHO": 5.0, "C2H4": 10.0, "C3H6": 5.0, "TOLUENE": 5.0,
                                  "HONO": 0.5})
    return fs.save_spec(spec, os.path.join(folder, "plume_ageing"))


# ── 3. Lagrangian air parcel ─────────────────────────────────────────────────────────────────────
def la_parcel(folder: str) -> str:
    n = 17                                                # every 15 min, 4 h, 10:00-14:00 PDT
    sec = np.array([utc_seconds("17:00") + 900 * i for i in range(n)])
    frac = np.arange(n) / (n - 1)
    lat = 33.95 + (34.14 - 33.95) * frac; lon = -118.45 + (-118.13 + 118.45) * frac     # Santa Monica Bay -> Pasadena
    T = 293 + 9 * frac; p = 1010 - 25 * frac; RH = 70 - 30 * frac
    h = 400 + 800 * (1 - np.cos(np.pi * frac)) / 2
    j = jno2_clear(sec, lat, lon, cloud=0.95)
    urban = 0.3 + 0.7 * np.sin(np.pi * np.clip(frac, 0, 1))      # the parcel crosses the dense urban core mid-way
    header = COMMON.format(title="Lagrangian air parcel: 4 h from Santa Monica Bay to Pasadena on the sea breeze") + """
# EMISSIONS are the inventory totals of the 8 x 8 km box around the parcel, mol/s (case: area_km2 = 64); they
# are mixed into the boundary layer h and held from one row to the next.  The growth of h dilutes the parcel
# with the background air of the B_ columns (case: dilution = 'pbl'); when h stops growing the rate of the
# case (0.3 per hour) applies.  Deposition: O3 0.4 cm/s and HNO3 2 cm/s over the layer (case).
# Initial values (case) are the air at the coast.  lat/lon per row give the sun's position along the way.
#"""
    cols = {"time": [stamp(s) for s in sec], "lat": lat, "lon": lon, "T": T, "p": p, "RH": RH, "h": h, "jNO2": j,
            "E_NO [mol/s]": 4.0 * urban, "E_NO2 [mol/s]": 0.6 * urban, "E_CO [mol/s]": 35.0 * urban,
            "E_C5H8 [mol/s]": 0.02 + 0.10 * frac, "E_TOLUENE [mol/s]": 0.5 * urban, "E_C2H5OH [mol/s]": 2.0 * urban,
            "B_O3 [ppb]": np.full(n, 35.0), "B_CO [ppb]": np.full(n, 150.0), "B_NO2 [ppb]": np.full(n, 1.0), "B_HCHO [ppb]": np.full(n, 1.0)}
    os.makedirs(folder, exist_ok=True)
    table = os.path.join(folder, "la_parcel.csv")
    write_table(table, header, cols, {"jNO2": ".3e", "T": ".1f", "p": ".1f", "RH": ".1f", "h": ".0f", "lat": ".4f", "lon": ".4f"})
    spec = fs.FieldSpec(name="la_parcel", run_type="lagrangian", table=table, site="Santa Monica Bay to Pasadena (synthetic)",
                        description="Synthetic sea-breeze trajectory with inventory-style emissions, boundary-layer growth, background and deposition.",
                        chemistry={"base_mechanism": "MCM_v3.3.1", "addons": [], "precursors": ["C5H8", "TOLUENE", "C2H5OH"], "mechanism_stem": "",
                                   "inorganic_only": False, "include_dirs": [], "reaction_adjustments": []},
                        step_s=60, default_emission_unit="mol/s", area_km2=64.0, dilution="pbl", k_collapse_per_h=0.3,
                        initials={"O3": 25.0, "NO": 0.5, "NO2": 2.0, "CO": 160.0, "HONO": 0.015, "HCHO": 1.0, "C5H8": 0.1, "TOLUENE": 0.2,
                                  "C2H5OH": 2.0},
                        deposition_velocity_cm_s={"O3": 0.4, "HNO3": 2.0})
    return fs.save_spec(spec, os.path.join(folder, "la_parcel"))


EXAMPLES = {"site_day": site_day, "plume_ageing": plume_ageing, "la_parcel": la_parcel}


def main() -> int:
    for key, make in EXAMPLES.items():
        path = make(os.path.join(HERE, key))
        spec = fs.load_spec(path)
        _table, rep = fs.check(spec)
        print(key, "ok" if rep.ok else rep.text())
        if not rep.ok:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

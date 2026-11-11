# JUICE from code

Everything the JUICE windows do can be done from Python. The windows and the code are not two
programs: the windows fill in the same experiment description and call the same functions. An
experiment built in code and the same experiment built by clicking give byte-identical input files.

There are two ways in:

- **Start from the window.** Every experiment the window prepares also gets a `<name>_run.py` next to its
  files. That file is the experiment written out as code. Run it, change it, copy it.
- **Start from scratch.** Write a script with `import juice`, as in the examples below.

All of it is in one module, [`juice/api.py`](../juice/api.py), which `import juice` gives you. The examples are in [`examples/`](../examples).
[`run_path.md`](run_path.md) follows one run through the files it calls, and shows where you can step in.

## 1. Setup

JUICE's own Python environment (the one the windows use) has everything. Install JUICE into it once, in
place, from the JUICE folder:

```
pip install -e .
```

After that `import juice` works in any script, in any folder. Without the install, a script has to say
where the JUICE folder is:

```python
import sys
sys.path.insert(0, r"C:\path\to\JUICE")   # the folder that holds the 'juice' package
import juice
```

Run scripts with that environment's Python, e.g. `<venv>\Scripts\python.exe my_experiment.py`.

## 2. A first experiment

Isoprene in SAPHIR under the open roof, with invented conditions
([`examples/chamber/saphir_isoprene.py`](../examples/chamber/saphir_isoprene.py)):

```python
import juice

exp = juice.new_experiment("saphir_isoprene", chamber="SAPHIR")

c = exp.conditions                         # when and under which conditions
c.date, c.start, c.end = "2024-07-09", "06:00", "16:00"
c.temperature_K, c.temperature_mode = 298.0, "const"
c.pressure_hPa, c.h2o_percent, c.flow_m3h = 1013.0, 1.0, 0.5
c.roof_open = [("08:00", "15:00")]
c.j_source = "clear_sky"

exp.precursors = ["C5H8"]                  # the MCM subset of isoprene is cut from the full MCM

exp.injections = [
    juice.Injection(name="NO", products="NO", total_ppb="5", start="07:30"),
    juice.Injection(name="isoprene", products="C5H8", total_ppb="10", start="08:30", duration_min=5),
]

run = juice.run(exp)                       # writes the input files and runs the solver in runs/saphir_isoprene
print(run.summary())
run.results().plot(["O3", "NO", "NO2", "C5H8", "OH", "HO2"])
```

`new_experiment` does what picking a chamber on the first step does: it sets that chamber's defaults (its
volume, wall sources, dilution, light). Everything after it is a plain field of the experiment, and a field
you do not set keeps that default. `juice.run(exp)` puts the run in `<output_root>/<name>` (`runs/` unless
`settings.local.toml` says otherwise); `juice.run(exp, r"D:\my runs\day_1")` puts it anywhere.

Every example comes twice: the short file (as above), and `<name>_explained.py`, the same experiment step by
step with what each setting means and what else can be chosen. Both write the same input files.

The experiment has four parts, the same as the window's steps:

| part | in code | what it holds |
|---|---|---|
| Where and when | `exp.conditions` | date, time window (UTC), output step, T, p, humidity, flow, roof or lights, J-value source |
| Chemistry | `exp.base_mechanism`, `exp.addons`, `exp.precursors`, `exp.initials` | the mechanism and what is in the chamber at the start |
| Injections | `exp.injections` | one `juice.Injection` per event: what, total amount, start, duration |
| Chamber | `exp.chamber` | what the chamber adds on its own: HONO and HCHO wall sources, background reactivity, dilution, wall deposition |

The fields that take a keyword, and what each keyword means:

| field | choices |
|---|---|
| `conditions.source` | `manual` (the constants you set), `archive` (the SAPHIR data archive of that day), `file` (a CSV or Excel table with T, p, H2O, J values, ...) |
| `conditions.temperature_mode` | `enz` (temperature read per time step), `const` (`temperature_K`) |
| `conditions.j_source` | `measured` (J values per time step), `clear_sky` (the mechanism's clear-sky formulas) |
| solver | `fast`, `precise`, `twostep` (Sect. 5) |

`juice.OPTIONS` holds the same list in code.

## 3. The four steps

`juice.run` is four calls in a row. They can be made one by one, for example to look at the files before
the solver runs:

```python
files = juice.write(exp, out_dir)                    # <name>.enz, <name>.easy, chemistry/, <name>.experiment.json
json  = juice.prepare(exp, out_dir, solver="fast")   # the solver input <name>.json
cs    = juice.solve(json)                            # runs the solver; the log goes to <name>_run.log
res   = juice.results(cs, exp)                       # the results as a table
```

`juice.write` runs the checks of the review panel first and refuses an experiment with errors. Call
`juice.check(exp)` to see the messages yourself; each one has a level (`error`, `warning`, `info`).

`<name>.easy` and `<name>.enz` are ordinary FACSIMILE input files, so the same experiment runs in FACSIMILE
unchanged.

## 4. Results

```python
r = run.results()             # or juice.results(path_to_cs, exp)
r.species                     # every output column
r["OH"]                       # molecules cm-3, one value per output time
r.ppb("O3")                   # ppb, with the model's air density (the M column of the results)
r.hour                        # hours since midnight UTC of the first output time
r.utc                         # the output times as datetimes
r.df                          # everything as a pandas DataFrame
r.plot(["O3", "OH", "C5H8"])  # one panel per species (radicals in molecules cm-3, the rest in ppb);
                              # the hours with the roof open / lights on are shaded; path="fig.png" saves it
```

For a field case, `juice.skill(case, cs_path)` compares the model with the `O_` columns of the table
(bias, normalised mean bias, RMSE, R² per species) and writes `<name>_skill.csv`. `juice.run` writes it
for every field case, as the field window's Run step does (`run.skill_path`).

The results open in the results viewer like any other run: `python -m juice.gui.results_viewer <name>.cs`.

## 5. What can be chosen

```python
juice.chambers()       # {'SAPHIR': ..., 'generic_lamps': ..., 'EUPHORE': ..., ...}
juice.mechanisms()     # {'bases': {'MCM_v3.3.1': ...}, 'addons': {...}}
juice.solvers()        # {'fast': ..., 'precise': ..., 'twostep': ...}, with what each costs
juice.templates()      # the templates of the first step
```

The three solver profiles solve the same equations and agree closely with FACSIMILE. They differ in run time
and in how closely they follow FACSIMILE:

| profile | method | use it for |
|---|---|---|
| `fast` (default) | BDF with the exact Jacobian, relative tolerance 1e-4, warm continuation across ENZ rows | everyday runs |
| `precise` | the same at 1e-6, cold restart at every ENZ row | the closest agreement with FACSIMILE |
| `twostep` | Verwer's TWOSTEP, no Jacobian, relative tolerance 1e-3 | very large mechanisms (GECKO-A); also the fastest for the complete MCM |

Other chemistry is one field away:

```python
exp.base_mechanism, exp.precursors = "CRACMM2", []   # CRACMM2 | CRACMM1_aq | RACM2: used as a whole
exp.addons = ["NOAA_BB"]                             # with the MCM: furans, furfurals, phenolics
```

To start from a template instead of an empty experiment:

```python
exp = juice.new_experiment("my_run", template="voc_photooxidation_low_nox")
```

## 6. From the window to code

When the Run step prepares an experiment, the run folder gets `<name>_run.py` besides the ENZ, EASY and
JSON files. The field window does the same. The file follows the window's steps and writes what the
experiment sets: every field that differs from the chamber's defaults, with what it means. Shortened:

```python
def build():
    # STEP 1 — EXPERIMENT: which chamber, when, and where the conditions come from
    exp = juice.new_experiment('saphir_isoprene', chamber='SAPHIR')
    c = exp.conditions
    c.date = '2024-07-09'    # YYYY-MM-DD, UTC
    c.start = '06:00'        # HH:MM UTC
    c.temperature_mode = 'const'    # const = temperature_K | enz = the T column of the file/archive
    ...
    # STEP 2 — CHEMISTRY: which mechanism, which compounds
    exp.base_mechanism = 'MCM_v3.3.1'
    exp.precursors = ['C5H8']
    ...
    # STEP 3 — INJECTIONS     exp.injections = [juice.Injection(...), ...]
    # STEP 4 — CHAMBER        (only what differs from the profile, e.g. ch.hono_const = '2.6*0.7*1.5')
    return exp

def main():
    exp = build()
    # STEP 5 — RUN
    solver = 'fast'          # fast | precise | twostep | None (the settings file's)
    kv_pd = True             # rates and production/destruction for Process Analysis
    files = juice.write(exp, OUT, TOML)                     # 1  input files (edit them here if you must)
    json_path = juice.prepare(exp, OUT, solver=solver, toml=TOML, kv_pd=kv_pd, files=files)
    cs_path = juice.solve(json_path, TOML)                  # 3  the solver
    r = juice.results(cs_path, exp)                         # 4  results as a table
```

- `python <name>_run.py` rebuilds the experiment, writes its files into `OUT`, runs the solver and prints
  a summary.
- `build()` returns exactly the experiment the window held, field by field. The test suite checks this
  for every template and example and through both windows.
- A field the file does not mention has the chamber's default (the date and the time window are always
  written). Section 2's table and `juice.OPTIONS` list the fields; the `_explained.py` examples show them all.
- Running the file again in its own folder overwrites that run. Change `OUT` to keep both.

[`examples/chamber/saphir_isoprene_explained.py`](../examples/chamber/saphir_isoprene_explained.py) is the same structure written by hand,
with longer notes: where a conditions file goes and how to hold NO2 or O3 to it, injections from a liquid
volume or a flow, the solver choice, and an off-road section (loops, your own chemistry, editing the input
files, another chamber).

A saved experiment can also be loaded directly:

```python
exp = juice.load(r"D:\runs\my_day\my_day.experiment.json")   # or a .field.json
exp.name = "my_day_no_hono"
exp.chamber.hono = False
juice.run(exp, r"D:\runs\my_day_no_hono")
```

## 7. Going further

Code makes it easy to do things the window does one at a time.

**Loops.** For example the same day at three NO levels (`python examples/chamber/saphir_isoprene_explained.py --scan`):

```python
for no in (1, 5, 20):
    e = isoprene_day(f"isoprene_NO{no}", no)
    r = juice.run(e, os.path.join(OUT, e.name), kv_pd=False).results()
    print(no, r.ppb("O3")[-1])
```

For large sweeps with parallel workers and summary tables, use the study runner
(`python -m juice.core.study_runner <name>.study.json`).

**Changing a rate constant.** Give the reaction exactly as the mechanism writes it, the new expression, and a
comment. `"0"` removes the reaction:

```python
exp.reaction_adjustments = [("OH + NO2 --> HNO3", "1.1*KMT08", "test: OH + NO2 10 % faster")]
```

**Your own chemistry.** Put an EASY include file in a folder, name the folder and the file, and add the
species it needs at the start:

```python
exp.include_dirs = [r"D:\my chemistry"]
exp.extra_includes = ["my_extra_reactions.incl"]
```

**Verbatim EASY.** Statements in `exp.overrides` are written into the EASY file as they are, after the
mechanism and the rate-constant changes. A later definition replaces an earlier one.

**Holding species to measurements.** `c.constrain_species = ["NO2", "O3"]` holds them to their columns in the
input table (`c.source = "file"` or `"archive"`).

**A lamp chamber.** Use `chamber="generic_lamps"` or a chamber profile of your own (made with the chamber
editor). Set `c.lights_on = [("09:00", "13:00")]` (the same field as `roof_open`). The lamps (jNO2 at full
power, spectrum) are set in the profile, not in the experiment, and `c.j_source` does not change them. See
[`examples/chamber/uv_lamp_chamber_explained.py`](../examples/chamber/uv_lamp_chamber_explained.py): the same steps as
the base case, with the light, the OH source, the walls and a profile of your own written out.

**Field cases.** A site, a plume or a Lagrangian air parcel:

```python
case = juice.new_field_case("my_site", run_type="constrained", table=r"D:\data\site.csv")
case.chemistry = {"precursors": ["C5H8", "TOLUENE"]}
case.generic_lifetime_h = 6.0
juice.run(case, r"D:\runs\my_site")
```

The input table is described in the main tutorial (Sect. 3.2).
[`examples/field/field_site_explained.py`](../examples/field/field_site_explained.py) writes the shipped `site_day` case out in the field
window's steps, with the table's columns, the dilution choices and the skill table. The field window, too, writes `<name>_run.py`
next to every case it prepares, and `juice.load("<case>.field.json")` reads a saved or shipped case
(`field_spec.list_examples()` lists the three examples).

## 8. Functions

| function | what it does |
|---|---|
| `new_experiment(name, chamber="SAPHIR", template="")` | a chamber experiment with the chamber's defaults, or a template |
| `new_field_case(name, run_type, table)` | a field case of one kind of run |
| `load(path)`, `save(spec, path)` | `.experiment.json` / `.field.json` |
| `check(spec)` | the review panel's checks: `[(level, message)]` |
| `write(spec, out_dir)` | the input files; refuses an experiment with errors |
| `prepare(spec, out_dir, solver=None, kv_pd=True)` | + the solver input `<name>.json` |
| `solve(json_path)` | runs the solver; returns the results file `<name>.cs` |
| `results(cs_path, spec=None)` | the results as a table (`Results`) |
| `skill(case, cs_path)` | field case: the model against the `O_` columns; writes `<name>_skill.csv` |
| `run(spec, out_dir="", solver=None, kv_pd=True)` | write, prepare, solve (`out_dir` '' = `<output_root>/<name>`); returns a `Run` with `.results()` and `.summary()` |
| `Results.plot(species, path="")` | one panel per species, light periods shaded; returns the figure |
| `write_script(spec, path, solver=None)` | the experiment as a runnable `<name>_run.py` |
| `build_from_script(path)` | `build()` of such a file, without running it |
| `chambers()`, `mechanisms()`, `solvers()`, `templates()`, `OPTIONS` | what can be chosen |
| `paths()` | the folders of the settings file (`include_dir`, `output_root`) |

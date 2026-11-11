# JUICE — Jülich Utility for Integrating Chemical Equations

Box-model simulations of simulation-chamber experiments and field observations with the Master Chemical
Mechanism (MCM) and other mechanisms. JUICE replaces the IDL/EASY + FACSIMILE workflow of the SAPHIR chamber
(Forschungszentrum Jülich) with Python, and reproduces FACSIMILE on the SAPHIR reference experiments.

- **Describe the experiment, not the files.** Chamber, conditions, chemistry, injections: from that one
  description JUICE writes the ENZ and EASY input files (they run in FACSIMILE unchanged) and solves them.
- **Any chamber, sun or lamps.** Chamber profiles for SAPHIR and about 18 other chambers; photolysis under
  UV lamps from the lamp spectrum.
- **Field observations.** A constrained site, a plume, or a Lagrangian air parcel, from one input table.
- **More than the MCM.** CRACMM 2, CRACMM 1, RACM2 and a biomass-burning add-on from the F0AM library, and
  GECKO-A generator output up to a million reactions.
- **Results and studies.** Process analysis, OH reactivity, comparison with measurements, figure lists;
  sweeps and isopleths over many runs.

## Install

```powershell
git clone https://github.com/GGkatzelis/JUICE.git
cd JUICE
python -m venv .venv
.venv\Scripts\activate
pip install -e .                                   # JUICE and its dependencies (Python >= 3.11)
python -m juice.core.experiment_spec --selftest    # quick check: "2/2 golden specs pass"
python -m juice.tools.make_shortcut                # the JUICE icon on the Desktop (Windows)
```

macOS / Linux: `source .venv/bin/activate`. Windows is the tested platform for the windows; the solver and
the code interface are plain Python.

## Use it

**With the windows:** double-click the JUICE icon (or `JUICE.bat`). It asks whether you model a chamber
experiment or field observations and walks you through five steps. [docs/tutorial.md](docs/tutorial.md) goes
through a chamber experiment, the same under UV lamps, a field case and a study in an afternoon.

**From Python:** the same in a script ([examples/chamber/saphir_isoprene.py](examples/chamber/saphir_isoprene.py)):

```python
import juice

exp = juice.new_experiment("saphir_isoprene", chamber="SAPHIR")
c = exp.conditions
c.date, c.start, c.end = "2024-07-09", "06:00", "16:00"
c.roof_open = [("08:00", "15:00")]
c.j_source = "clear_sky"                                        # J values from the sun's position
exp.precursors = ["C5H8"]                                       # isoprene: the MCM subset is cut for you
exp.injections = [juice.Injection(name="NO", products="NO", total_ppb="5", start="07:30"),
                  juice.Injection(name="isoprene", products="C5H8", total_ppb="10", start="08:30")]

run = juice.run(exp)                                            # input files + solver, in runs/saphir_isoprene
run.results().plot(["O3", "NO2", "C5H8", "OH"])
```

More in [examples/](examples) (each one short and step by step) and [docs/juice_from_code.md](docs/juice_from_code.md).
Every run the windows prepare also gets a `<name>_run.py` that rebuilds it in code.

## What is where

| folder | what it holds |
|---|---|
| [juice/](juice) | the code: `core/` (experiment → input files → solver), `chemistry/`, `importers/`, `gui/`, `tools/` |
| [examples/](examples) | chamber and field examples as scripts, and three synthetic field cases |
| [docs/](docs) | tutorial, user guide, JUICE from code, architecture, chamber survey, lamp photolysis |
| [chemistry/](chemistry) | `mechanisms/` (MCM v3.3.1 and the others, add-ons) and `photolysis/` (cross sections, J parametrisations) |
| [setups/](setups) | what an experiment starts from: `chambers/` (one profile per chamber), `templates/`, `mixtures/` |
| [validation/](validation) | the FACSIMILE reference scenes and the UV-lamp reference |
| [tests/](tests) | the test suite: `pip install -e .[dev]`, then `python -m pytest` |
| `runs/` | your experiments and results (created on the first run, not in git) |

`settings.toml` holds the validated solver settings; your own paths go in `settings.local.toml` (not in git).

## Validation

`python -m juice.tools.run_benchmark` runs every SAPHIR reference scene through JUICE and compares every
column with FACSIMILE. With the default solver profile, every setting agrees and every species above 1 ppt
is within 0.6 % in the six SAPHIR experiments (O3, NOx, OH and HO2 within 0.25 %). The tables, the solver
profiles (fast, precise, twostep) and the known limitations are in the
[user guide](docs/user_guide.md#validation) and in [docs/architecture.md](docs/architecture.md).

## Chemistry from others

JUICE ships chemistry that others developed. Please cite it when you use it.

| chemistry | source | please cite | terms |
|---|---|---|---|
| MCM v3.3.1, complete, as downloaded (checksums in its manifest) | [mcm.york.ac.uk](https://mcm.york.ac.uk/MCM/) | Jenkin et al. (1997); Saunders et al. (2003); Jenkin et al. (2003); Bloss et al. (2005); Jenkin et al. (2015), as the [MCM citation page](https://mcm.york.ac.uk/MCM/about/cite) asks | the MCM website states no licence; it asks for the citation wording on that page |
| CRACMM 2, CRACMM 1 (aqueous), RACM2, NOAA biomass burning (add-on) | converted from the [F0AM](https://github.com/AirChem/F0AM) chemistry library (`Chem/`) | Skipper et al. (2024); Place et al. (2023); Goliff et al. (2013); Decker et al. (2019) and Coggon et al.; and F0AM: Wolfe et al., GMD 9, 3309 (2016) | F0AM: GPL-3.0 |
| GECKO-A ethanol (an imported example; larger ones are imported locally) | generator output of [GECKO-A](https://gitlab.in2p3.fr/) | Aumont et al. (2005); Valorso et al. (2011) | GPL-3.0-or-later (Zenodo) / AGPL-3.0 (GitLab); JUICE reads the output and ships none of its code |
| cross sections × quantum yields for lamp photolysis | NCAR [TUV-x](https://github.com/NCAR/tuv-x) (JPL, IUPAC evaluations) | as listed per J in `chemistry/photolysis/sigma_phi_298K.sources.md` | Apache-2.0 |

Each converted mechanism keeps a manifest with its source, version, citation and the conversion notes
(`chemistry/mechanisms/*/*/manifest.json`).

## Citation and collaboration

If you use JUICE in a publication, please cite it. GitHub's "Cite this repository" button gives the
reference from [CITATION.cff](CITATION.cff). Please also contact us for the corresponding methods reference.

JUICE is new research software. If you plan to use it for a publication, we strongly encourage you to
contact us early. We are happy to help check that the setup and the results are correct. Where we
contribute substantially to the analysis, co-authorship can be discussed in line with standard
authorship guidelines. This is a request, not a condition of the license.

## License

Copyright (C) 2026 Georgios Gkatzelis and Franz Rohrer, Forschungszentrum Jülich GmbH.

JUICE is free software: you can redistribute it and/or modify it under the terms of the
[GNU General Public License](LICENSE) as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version. It is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Contact

- Georgios Gkatzelis, Forschungszentrum Jülich GmbH — [g.gkatzelis@fz-juelich.de](mailto:g.gkatzelis@fz-juelich.de)
- Franz Rohrer, Forschungszentrum Jülich GmbH — [f.rohrer@fz-juelich.de](mailto:f.rohrer@fz-juelich.de)

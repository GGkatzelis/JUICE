# Examples

Every example comes twice:

| file | what it is |
|---|---|
| `<name>.py` | the short version: only what the experiment sets, then run and a figure. Start here. |
| `<name>_explained.py` | the same experiment step by step, in the order of the window, with what every setting means, what else can be chosen, and an "off-road" section (loops, your own chemistry, editing the input files, your own chamber). |

Both write byte-identical input files (`tests/test_juice_api.py` checks it).

| folder | example | |
|---|---|---|
| [chamber/](chamber/) | `saphir_isoprene` | isoprene photo-oxidation in SAPHIR, open roof, 5 ppb NO |
| | `uv_lamp_chamber` | ethanol + OH in a small chamber with UV-A blacklights, H2O2 as the OH source |
| [field/](field/) | `field_site` | a constrained site: measured species held, the rest computed, HCHO compared |
| | `site_day/`, `plume_ageing/`, `la_parcel/` | the three synthetic field cases the field window offers under "Load an example" ([field/README.md](field/README.md)) |

Run one with the Python of JUICE's environment, from the JUICE folder:

```
python examples/chamber/saphir_isoprene.py
```

The run goes to `runs/<name>/` (or wherever `settings.local.toml` points). Once `pip install -e .` has been
done, the scripts also run from any other folder; copy one there and change it for your own experiment.
[docs/juice_from_code.md](../docs/juice_from_code.md) is the reference for everything `import juice` offers.

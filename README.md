# Minimal Cantera QSSA NH2 reconstruction

This directory is a small, standalone reproduction of one case from the
parent project: `cases/NH3_KAUST_DEF_NP_1bar/sim_props.yaml`.

It deliberately contains only three steps:

1. solve one NH3/H2 counterflow flame with Cantera;
2. reconstruct `NH2` twice with a fixed-temperature QSSA integration;
3. compare the flame result with QSSA baseline and QSSA baseline + OH.

There is no PyCSP, CSP/CEM calculation, NSP workflow, Fortran executable,
parallel scheduler, or batch/multi-mechanism logic.

## Requirements

Use Python 3.10 or newer, create an environment, and install the dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

The NUIG 39-species mechanism used by this example is included at
`mechanism/NUIG_39s306r.yaml`.  Its authorship and recommended citation are in
[`mechanism/README.md`](mechanism/README.md).

## Run

From this directory, run:

```bash
make run
```

The generated files are:

```text
results/simulation.csv
results/qssa_baseline.csv
results/qssa_baseline_plus_OH.csv
results/nh2_comparison.csv
results/nh2_error_summary.csv
results/nh2_comparison.png
```

To use a different mechanism:

```bash
make run MECHANISM=/path/to/mechanism.yaml
```

The equivalent individual commands are:

```bash
python3 src/simulate_counterflow.py --mechanism /path/to/mechanism.yaml
python3 src/reconstruct_nh2.py --mechanism /path/to/mechanism.yaml
python3 src/plot_nh2_comparison.py
```

For a quick installation check before the full reconstruction, add
`--max-rows 3` to the reconstruction command.  On the development machine the
serial full reconstruction takes roughly one minute; `--max-rows` is not a
scientific result and should not be used for the final comparison.

## QSSA convention

At every flame-grid point, the species in `fixed_species` are prescribed from
the Cantera flame.  All other species are initialized to zero, then integrated
at the local temperature and pressure while the prescribed species are held
fixed.  The resulting `X_NH2_qssa` is the reconstructed NH2 mole fraction.

The two supplied variants are:

| Variant | Prescribed species |
| --- | --- |
| `baseline` | NH3, H2, O2, H2O, N2 |
| `baseline_plus_OH` | baseline species + OH |

The scripts are intentionally serial and transparent.  This makes the
numerical definition easy to inspect; it is not intended as a replacement for
the parent project's high-throughput reconstruction workflow.

## License

The source code and supplied results are released under the [MIT License](LICENSE).
The NUIG kinetic mechanism retains its separate authorship and citation; see
[`mechanism/README.md`](mechanism/README.md).

## Download QR code

<img src="assets/github_qr.png" alt="QR code for the GitHub repository" width="280">

The QR code encodes: <https://github.com/zhuxu317/QSSARecon>

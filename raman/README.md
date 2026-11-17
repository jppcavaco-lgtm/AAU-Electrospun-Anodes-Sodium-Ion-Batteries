# Raman figures for electrospun carbon anodes

Supplementary analysis package for Joao Pacheco Pereira Cavaco's MSc thesis,
*Electrospun electrodes for Sodium-Ion Batteries*, FCUL, Universidade de Lisboa.

The package contains 15 measured Raman spectra, the saved five-band fits used for
the thesis figures, and portable copies of the plotting and fitting scripts.
The final thesis discusses the spectra qualitatively. The fits illustrate
overlapping contributions. Their parameters and covariance-based errors are
exploratory numerical outputs.

## Folder contents

- scripts/raman_integracao.py - baseline correction, observed D/G maxima and five-band refitting.
- scripts/raman_picos.py - separate baseline-corrected plots from one results folder.
- scripts/apresentacao_raman.py - paired baseline/fit figures, including residuals.
- data/ - copies of the original two-column TXT spectra, grouped by material.
- samples.json - input manifest linking each spectrum ID to a relative data path.
- saved_fits/ - the saved numerical arrays and fitting metadata used for the thesis.
- presentation_selection.json - the spectrum selected for each main-text figure.
- provenance.json - SHA-256 hashes of the source scripts, spectra and saved arrays.
- Raman_Scripts_Documentation.txt - detailed processing, function and interpretation notes.
- requirements.txt - package versions recorded for the original analysis.

## Installation

Use Python 3.11 or newer. The original analysis used Python 3.11.4.
Open a terminal inside this raman folder, then run:

    python -m pip install -r requirements.txt

An isolated Python virtual environment can be used. Run all commands below with
the Python interpreter in which these dependencies are installed.

## Recreate the paired figures used in the thesis

    python scripts/apresentacao_raman.py --output results/thesis_figures

This reads the supplied saved fits. It creates 15 PNG/SVG figure pairs, covering
all main-text and additional spectra. Both panels always show the same spectrum.
The selection JSON identifies main-text versus appendix figures.
The output directory must be new. Omitting --output generates a timestamped folder.

## Recalculate the fits from the measured TXT spectra

    python scripts/raman_integracao.py --output results/refit

This recalculates all 15 fits with the archived numerical settings and random seed.
It writes JSON, CSV, NPZ and grouped fit plots. To recalculate one material:

    python scripts/raman_integracao.py --group ffes_ca_cb10 --output results/refit_cb10

Group IDs are listed in samples.json. A new Python/SciPy environment may produce
slightly different fitted values for overlapping bands.

## Plot the recalculated data

Separate baseline-corrected plots:

    python scripts/raman_picos.py --results results/refit

Paired baseline-corrected/five-band figures:

    python scripts/apresentacao_raman.py --results results/refit --output results/refit_figures

For a one-group refit, replace results/refit with results/refit_cb10.
The presentation script retains the thesis spectrum selection in
presentation_selection.json instead of selecting a new spectrum from a new fit.

## Scientific context

Four Lorentzian components (D4, D1, G, D2) and one Gaussian component (D3) follow
the model combination examined by Sadezky et al. for soot and related carbonaceous
materials. Numerical limits, baseline windows and the selection rule were choices
made in this analysis. Several saved fits have parameters at their imposed bounds.
An excellent fit or high R-squared does not by itself validate the band assignments.

Recorded centre/width/height standard errors and ratio errors describe local
least-squares behaviour under model assumptions. They do not provide a full
measurement uncertainty. Daily calibration could not be confirmed. Measured
spectra at different locations are not documented independent fabrication batches.
The data do not support universal ranking of graphitization, carbon-phase
fractions, porosity or conductivity through fitted D/G ratios.

See Raman_Scripts_Documentation.txt for the full explanation and references.

## Changes from the workspace scripts

The numerical model and saved arrays were retained. File discovery now uses
samples.json with relative paths instead of a personal Windows/OneDrive directory.
The paired-figure script accepts input/output arguments and uses the archived
main-text selection. Original input files and the earlier scripts remain in the
author's workspace. The old three-band documentation has been replaced in this
package by documentation of the final five-band qualitative workflow.

## Adding this package to the supplementary GitHub repository

Place this entire raman folder at the repository root, alongside the existing
supplementary G-code. Preserve its subfolders. The root README can link to
raman/README.md. Generated results and virtual environments are ignored by
.gitignore. The measured data and saved fits are included so that another reader
can reproduce the figures without the author's local directory structure.

The package was prepared locally. Repository visibility and upload were not
verified. No publication licence is assigned by this package.

## References

Sadezky et al. (2005), Carbon 43, 1731-1742.
https://doi.org/10.1016/j.carbon.2005.02.018

Ferrari and Robertson (2000), Physical Review B 61, 14095-14107.
https://doi.org/10.1103/PhysRevB.61.14095

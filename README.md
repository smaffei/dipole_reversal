# Repository for dipole reversal study.

### Contents
- Models (all spherical harmonic models used)
  - IGRF
  - IMM
  - CHAOS5
  - CALS10K1B
- instantaneous_flows (original study by Sam/Phil)
  - original code (original fortran code files provided by Phil)
  - modified code (modified fortran code to be called by python)
  - dipole_bound.py (module containing required functions)
  - igrf.ipynb (IGRF analysis)
  - imm.iypnb (IMM analysis)
  - IMM_roc.txt (data file from analysis output)
  - IMM_roc_ex.txt (data file from analysis output)
  
-----

### Instantaneous flows

Flows are optimised for the instantaneous Rate Of Change (ROC) of the axial dipole (g10). ROC is calculated for unconstrained and for purely toroidal flows given a magnetic field and an RMS flow velocity.

#### Using the code

To compile the fortran code (requires fftw3):
1. First navigate to `instantaneous_flows/modified_code` and change the path to the fftw3 library within the Makefile to your own location, then run the Makefile with `make`.

To setup the python environment:
1. To create the conda environment for the python code, run `conda env create -f environment.yml` from within `instantaneous_flows`
2. Activate the conda environment with `source activate dipole_reversal`
3. Finally add the kernel to those available to jupyter notebooks with `python -m ipykernel install --user --name dipole_reversal --display-name "Python (dipole_reversal)"`
4. The jupyter notebooks can now be run with `jupyter notebook` (you may need to select the kernel you have just installed)

-----

# IMP (Intrinsic membrane properties)

MATLAB script to calculate intrinsic membrane properties from whole cell patch clamp (current clamp) data.

## Data format and organisation
Traces can be in .abf (Axoclamp/Clampex) or .dat (HEKA Patchmaster) formats. HEKA .dat files require exporting each protocol into .mat file beforehand using Patchmaster.

## Required functions
* import_abf (reads .abf file)
* abfload (F. Collman) **third party**
* import_heka (reads .mat file)
* rmp (calculates resting membrane potential)
* iv (calculates Rin, tau, capacitance and single spike parameters, if depolarization reaches threshold)
* Ezyfit toolbox (F. Moisy; http://www.fast.u-psud.fr/ezyfit/) **third party**
* sag
* f_Icurve (calculates firing parameters. Single spike parameters also calculated if I/V protocol contained no spikes

## Typical protocols for optimal analysis
All protocols are optional and parameters will be set to NaN in the results if the relative protocol is not present
* rmp: 20s-60s current clamp trace no current injection
* I/V: increasing current injections (dI: 5 or 10 pA, 0.4-1s length) starting from -50pA. At least 7-8 sweeps without spikes
* sag: long (>=1s) and strong hyperpolarizing current injection (-100/-150 pA or enough current to bring the cell between -100 and -120 mV, depending on cell types)
* f/I curve: progressively depolarizing current injections (typically
* from +20 pA, dI: 50 pA)

## Saved results
All values are exported in an excel file with the recording name (extracted from I/V protocol). All missing values are set as NaN.

## Use
Run IMP.m pipeline, which runs all the required functions sequentially.

Marco Bocchio, last update 17/4/2019

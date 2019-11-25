# IMP
Script to calculate Intrinsic Membrane Properties from whole cell patch clamp (current clamp) data.


Traces can be in .abf (Axoclamp/Clampex) or .dat (HEKA Patchmaster)
formats. HEKA .dat files require exporting each protocol into .mat file beforehand using
Patchmaster.

# REQUIRED FUNCTIONS
   - import_abf (reads .abf file)
   - import_heka (reads .mat file)
   - importHEKA (Malcolm Lidierth)
   - abfload (F. Collman)
   - rmp (calculates resting membrane potential)
   - iv (calculates Rin, tau, capacitance and single spike parameters, if
       depolarization reaches threshold)
   - Ezyfit toolbox (F. Moisy; http://www.fast.u-psud.fr/ezyfit/)
   - sag
   - f_Icurve (calculates firing parameters. Single spike parameters also
       calculated if I/V protocol contained no spikes


# TYPICAL PROTOCOLS FOR OPTIMAL ANALYSIS
 All protocols are optional and parameters will be set to NaN in the results
 if the relative protocol is not present

   - rmp: 20s-60s current clamp trace no current injection
   - I/V: increasing current injections (dI: 5 or 10 pA, 0.4-1s length) starting from
       -50pA. At least 7-8 sweeps without spikes
   - sag: long (>=1s) and strong hyperpolarizing current injection
       (-100/-150 pA or enough current to bring the cell between -100 and
       -120 mV, depending on cell types)
   - f/I curve: progressively depolarizing current injections (typically
       from +20 pA, dI: 50 pA)
# SAVED RESULTS
 All values are exported in an excel file with the recording name
 (extracted from I/V protocol). All missing values are set as NaN


Marco Bocchio, updated 17/4/2019

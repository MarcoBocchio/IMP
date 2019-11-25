% IMP - Script to calculate Intrinsic Membrane Properties from whole
% cell patch clamp (current clamp) data.

%%
% Traces can be in .abf (Axoclamp/Clampex) or .dat (HEKA Patchmaster)
% formats. HEKA .dat files require exporting each protocol into .mat file beforehand using
% Patchmaster.

%%
% >>> REQUIRED FUNCTIONS >>>
%   - import_abf (reads .abf file)
%   - import_heka (reads .mat file)
%   - importHEKA (Malcolm Lidierth)
%   - abfload (F. Collman)
%   - rmp (calculates resting membrane potential)
%   - iv (calculates Rin, tau, capacitance and single spike parameters, if
%       depolarization reaches threshold)
%   - Ezyfit toolbox (F. Moisy; http://www.fast.u-psud.fr/ezyfit/)
%   - sag
%   - f_Icurve (calculates firing parameters. Single spike parameters also
%       calculated if I/V protocol contained no spikes

%%
% >>> TYPICAL PROTOCOLS FOR OPTIMAL ANALYSIS >>>
% All protocols are optional and parameters will be set to NaN in the results
% if the relative protocol is not present
%
%   - rmp: 20s-60s current clamp trace no current injection
%   - I/V: increasing current injections (dI: 5 or 10 pA, 0.4-1s length) starting from
%       -50pA. At least 7-8 sweeps without spikes
%   - sag: long (>=1s) and strong hyperpolarizing current injection
%       (-100/-150 pA or enough current to bring the cell between -100 and
%       -120 mV, depending on cell types)
%   - f/I curve: progressively depolarizing current injections (typically
%       from +20 pA, dI: 50 pA)

%%
% >>> SAVED RESULTS >>>
% All values are exported in an excel file with the recording name
% (extracted from I/V protocol). All missing values are set as NaN

%% 
% Marco Bocchio, updated 17/4/2019

%% 1. Pre-allocation

close all

fileFormat = input('Which file format would you like to load? Type 1 for ABF (default), type 2 for HEKA (.dat):','s');
if fileFormat ~= '2';
    fileFormat = 1;
end

rmpValue=NaN;
inputRes=NaN;
membraneTau=NaN;
membraneCapacitance=NaN;
spikeThreshold=NaN;
firstSpikeAmpl=NaN;
spikeHalfWidth=NaN;
fAHPampl=NaN;
sagRatio=NaN;
rebAmpl=NaN;
adaptIndex=NaN;
burstIndex=NaN;
firingCurve=NaN;
shortFiringCurve=NaN;
burstIndexCurve=NaN;
fileName = 'NaN_';
   
%% 2. Resting membrane potential

isThereInput = input('Is there a resting membrane potential protocol? Y/N [Y]:','s');
    if isThereInput=='N';
    else
        disp('select trace for resting membrane potential')
        [rmpValue] = rmp(fileFormat);
    end
    
    
%% 3. I/V    
isThereInput = input('Is there an I/V protocol? Y/N [Y]:','s');
if isThereInput=='N';
else
disp('select traces for I/V')
[inputRes, membraneTau, membraneCapacitance, spikeThreshold, firstSpikeAmpl, spikeHalfWidth, fAHPampl, fileName] = iv(fileFormat);
end

%% 4. Sag

isThereInput = input('Is there a sag protocol? Y/N [Y]:','s');
    if isThereInput=='N';
    sagRatio = 0;
    rebAmpl = 0;
    else
    disp('select trace for sag')
    [sagRatio rebAmpl] = sag(fileFormat);
    end

  
    singleSpikeParam = 0;
    if spikeThreshold == 0; %flag for no single spike parameters present
        singleSpikeParam = 1;
    end
  
%% 5. f/I curve

isThereInput = input('Is there a f/I protocol? Y/N [Y]:','s');
    if isThereInput=='N';
    else
         
    disp('select trace for f/I protocol')
        if singleSpikeParam == 1; %extraction of single spike parameters from protocol
        [firingCurve, shortFiringCurve, burstIndexCurve, burstIndex, adaptIndexCurve, adaptIndex, rheobase, spikeThreshold, firstSpikeAmpl, spikeHalfWidth, fAHPampl] = f_Icurve(fileFormat, singleSpikeParam);
        else 
        [firingCurve, shortFiringCurve, burstIndexCurve, ~, adaptIndexCurve, ~, rheobase] = f_Icurve(fileFormat, singleSpikeParam); %no extraction of single spike parameters from protocol (already measured in IV protocol)
        end
    end
   
   
%% 6. Saving results

toSave = [table(rmpValue), table(inputRes), table(membraneTau), table(membraneCapacitance), table(spikeThreshold), table(firstSpikeAmpl), table(spikeHalfWidth), table(fAHPampl), table(sagRatio), table(rebAmpl), table(adaptIndex), table(burstIndex), table(firingCurve), table(shortFiringCurve), table(burstIndexCurve)];
fileName = [fileName '.xls'];
writetable(toSave,fileName);

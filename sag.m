function [sagRatio rebAmpl] = sag(fileFormat);
%%
%OPERATION
% Calculates sag ratio and rebound amplitude for hyperpolarizing current
% injection
%
% >>> INPUT VARIABLES >>>
% NAME              TYPE, DEFAULT      DESCRIPTION
% fileFormat        scalar, 1          file format to be loaded: 1 for ABF, 2 for
%                                      HEKA (.mat)
%%
% >>> PARAMETERS >>>
% NAME              DEFAULT            DESCRIPTION
% sagSweep          1                  sweep number to use for sag
% currInjDuration   1000               length of current injection (in ms)
% base              131 (ABF)          baseline in ms
%                   100 (HEKA)                                             
% firstCurrAmpl     -50                amplitude of first current injection
% currInjDuration   500                length of current injection (in ms)
% SpikeThreshold    0                  Threshold to detect spikes (value in Y
% minSearch         300                time to search for minimum point after baseline (in ms)
% sagSearch         50                 time to average to measure sag amplitude at the end of the current injection (in ms)
% rebSearch         200                time to search for maximum point for rebound (in ms)
% plotting          1                  Plotting on/off


%% 
% Marco Bocchio, updated 17/4/2019


%% Loading files

%ABF file (default)
if nargin == 0;
    fileFormat = 1;
end

if fileFormat == 1;
    [Ytrace,time, Fs, si_ms] = import_abf;
else
    
%HEKA file
   [YtraceCh1,YtraceCh2,time, Fs, si_ms] = import_heka;
   Ytrace = YtraceCh1; %use trace from channel 2
end

%% Parameters
plotting = 1; %select 1 yes and 0 for no
sagSweep = 1; %sweep to use to measure sag
%si_ms = si.si./1000; %sampling interval converted to ms (from µs)
SpikeThreshold = 0; %spike threshold (trace value, not amplitude)
currInjDuration = 1000; %current injection duration in ms
sweepTime = (0:length(Ytrace)-1).*si_ms; %time of each sweep (in ms)
%Fs = si.si *0.1; %sampling rate (in kHz), obtained from sampling interval in abfload
base = 131; % baseline time (in ms)
currInjEnd = base + currInjDuration; %time of end of current injection (in ms)
minSearch = 300; %time to search for minimum point after baseline (in ms)
sagSearch = 50; %time to average to measure sag amplitude at the end of the current injection (in ms)
rebSearch = 200; %time to search for maximum point for rebound (in ms)

%% Sag and rebound depolarization
minHyperpSearchTrace = Ytrace(base*Fs:(minSearch*Fs),sagSweep);  % part of trace to search for min hyperpolarization point
sagSearchTrace = Ytrace(round(currInjEnd*Fs)-(round(sagSearch/2)*Fs):round(currInjEnd*Fs)+(round(sagSearch/2)*Fs),sagSweep); % part of trace to search for sag (end of hyperp step)
rebSearchTrace = Ytrace(round(currInjEnd*Fs):round(currInjEnd*Fs+(rebSearch*Fs)),sagSweep); %part of trace to search for rebound depolarization

sagTraceMean = mean(sagSearchTrace); % mean sag value (membrane value, not amplitude)
baseMean = mean(Ytrace(1:base*Fs)); % mean baseline
minHyperpValue = min(minHyperpSearchTrace); % mean min hyperpolarization value ((membrane value, not amplitude)
maxRebValue = max(rebSearchTrace); % mean max rebound depolariz value (membrane value, not amplitude)
minHyperpLoc = find(minHyperpSearchTrace==minHyperpValue,1); % location of min hyperpolarization value
maxRebLoc = find(rebSearchTrace==maxRebValue,1);  % location of max rebound depolarization value

minHyperpMean = mean(minHyperpSearchTrace(minHyperpLoc-1:minHyperpLoc+1)); % mean of 3 data points around the min hyperpolariz value
maxRebMean = mean(rebSearchTrace(maxRebLoc-1:maxRebLoc+1)); % mean of 3 data points around the max rebound depolariz value

hyperpAmplStart = baseMean - minHyperpMean; % hyperpolarization amplitude at the beginning of the current injection (HCN channels are closed)
hyperpAmplEnd = baseMean - sagTraceMean; % hyperpolarization amplitude at the end of the current injection (HCN channels are open)
sagRatio = hyperpAmplEnd / hyperpAmplStart; % sag ratio: values closer to zero imply stronger sag
rebAmpl = maxRebMean - baseMean; %rebound amplitude

minHyperpLocFullSweep = minHyperpLoc + (base*Fs);
maxRebLocFullSweep = maxRebLoc + (currInjEnd*Fs);

if maxRebMean > -20;
   warning('WARNING: rebound spike!')
   rebAmpl = 0;
end

if plotting == 1; %if plotting is set on 'yes'

    %% Plotting
    %close all;
    figure;

    plot(Ytrace(:,sagSweep));
    hold;
    scatter(minHyperpLocFullSweep,minHyperpMean,100,'r');
    scatter(maxRebLocFullSweep,maxRebMean,100,'r');
    title('Sag and rebound depolarization')
    ylabel('Membrane potential (mV)')
    xlabel('Data points')
    
else
    
end

%% Printing results
sagRatioDisplay = ['Sag ratio = ', num2str(sagRatio)];
disp (sagRatioDisplay);
rebAmplDisplay = ['Rebound amplitude = ', num2str(rebAmpl), ' mV'];
disp (rebAmplDisplay);

end


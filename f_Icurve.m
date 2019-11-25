function [firingRate, burstFiringRate, burstIndexCurve, burstIndex, adaptIndexCurve, adaptIndex, rheobase, spikeThreshold, firstSpikeAmpl, spikeHalfWidth, fAHPampl] = f_Icurve(fileFormat, singleSpikeParam);

%%
% OPERATION
% Calculates membrane and single spike parameters from I/V protocol
%
%%
% >>> INPUT VARIABLES >>>
% NAME                  TYPE, DEFAULT      DESCRIPTION
% fileFormat            scalar, 1          file format to be loaded: 1 for ABF, 2 for
%                                      HEKA (.mat)
% singleSpikeParam      scalar, 0          flag to calculate or not parameters from first spike 

%%
% >>> PARAMETERS >>>
% NAME                  DEFAULT             DESCRIPTION
% base                  100                 baseline in ms
% firstCurrAmpl         20                  amplitude of first current injection
% currInjDuration       500                 length of current injection (in ms)
% currInjDelta          50                  increasing current amplitude
% minSpikeDetInterval   2                   minimum time between spikes for
%                                           detection (in ms)
% burstIndexWindow      50                  time to calculate burst index (in ms)

% firstSpikeWindow      30                  time window to quantify and plot the first spike in the IV (in ms)
% SpikeThreshold_det    20                  Threshold to detect spikes (value in Y
%                                           trace, not amplitude)
% minSpikesFiringParam  15                  minimum number of spikes to
%                                           calculate firing parameters
% plotting              1                   Plotting on/off

%% 
% Marco Bocchio, updated 20/5/2019


%% Default input if no input is provided
if ~exist('singleSpikeParam', 'var')
    singleSpikeParam = 0; %do not calculate first spike parameters if no input is provided
end


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
%si_ms = si.si./1000; %sampling interval converted to ms (from µs)
%Fs = si.si *0.1; %sampling rate (in kHz), obtained from sampling interval in abfload
spikeThreshold_det = -30; %spike threshold (trace value, not amplitude)
firstCurrAmpl = 20; %current injected in first sweep (in pA)
currInjDuration = 500;
currInjDelta = 50; %delta used for current steps (in pA)
%time = (0:length(Ytrace)-1).*si_ms; %time of each sweep (in ms)
base = 100; %time of baseline in ms
basePoints = base .* Fs; %baseline in points
minSpikeDetInterval = 2; %minimum time between spikes for detection
minSpikeDetPnts = minSpikeDetInterval * Fs; %minimum number of points between spikes for detection
sweepN = size(Ytrace,2); %number of sweeps
burstIndexWindow = 50; %time to calculate burst index (usually 50 ms)
%binSize = 50; % bin size in ms
minSpikesFiringParam = 15; % minimum number of spikes to calculate firing parameters
firstSpikeWindow = 10; %time window that is taken to quantify and plot the first spike in the IV (in ms)
fAHPwindow = 5; %window to measure fAHP


% Time of the last spike
%lastSpikeLocs = max(spikePeaksLocs); %finds location of the last spike in each sweep
%lastSpikeTimes = (lastSpikeLocs ./ Fs) - baseTime; %time of the last spike in each sweep
%lastSpikeTimesWhenNoSpike = find(lastSpikeTimes,-baseTime);

%% Preallocation for some variables that may not be assigned
spikeThreshold=0;
firstSpikeAmpl=0;
spikeHalfWidth=0;
fAHPampl=0;

%% Spike detection
spikePeaks = zeros(200,size(Ytrace,2)); %preallocation
spikePeaksLocs = zeros(200,size(Ytrace,2)); %preallocation
lastSpikeLocs = zeros(size(Ytrace,2),1); %preallocation
spikeN = zeros(1,size(Ytrace,2)); %preallocation 
lastPartCurrInjTime = (base+400):(base+500);
lastPartCurrInjPnts = lastPartCurrInjTime .* Fs;
lastPartCurrInjTrace = Ytrace (lastPartCurrInjPnts,:);
lastPartCurrInjMean = mean(lastPartCurrInjTrace,1);
sweepsWithSpikes = zeros(1,size(Ytrace,2)); %preallocation

    for i = 1 : size(Ytrace,2);
        if lastPartCurrInjMean(i) > spikeThreshold_det - 2;
            spikeThreshold_det = lastPartCurrInjMean(i) + 5;
        end
       [spikePeaks_i,spikePeaksLocs_i] = findpeaks(Ytrace(1:(base+currInjDuration)*Fs,i),'MinPeakHeight',spikeThreshold_det, 'MinPeakDistance', minSpikeDetPnts); %finds spike peaks for each sweep
       noSpikes = isempty(spikePeaks_i); %checks if the sweep has no spikes
       if noSpikes == 0  %if there are spikes
          sweepsWithSpikes(i) = 1;     %value of vector corresponding to sweep is one
          spikePeaks(1:length(spikePeaks_i),i) = spikePeaks_i; %allocation of spike peaks values to matrix (each sweep is a column)
          spikePeaksLocs(1:length(spikePeaksLocs_i),i) = spikePeaksLocs_i;%allocation of spike peak locations values to matrix (each sweep is a column)
          spikeN(i) = length(spikePeaks_i);
          lastSpikeLoc_i = max(spikePeaksLocs_i);
       else
           spikeN(spikePeaks_i) = 0;
          lastSpikeLoc_i = NaN;
       end
       lastSpikeLocs(i) = lastSpikeLoc_i;
    end
    
firingRate = spikeN ./ (currInjDuration/1000); 
    
    
lastCurrInj = (((size(Ytrace,2)-1)*currInjDelta)+firstCurrAmpl);
current = [firstCurrAmpl:currInjDelta:lastCurrInj];

spikeTimes = spikePeaksLocs ./ Fs;

%% First sweep with spikes
%only used if calculation of first spike parameters is selected
    firstSweepWithSpikes = find(sweepsWithSpikes,1); %first sweep with spikes
    firstSpikePeakLoc = spikePeaksLocs(1,firstSweepWithSpikes); %location of the first spike in the first sweep with spikes
    firstSpikePeak = spikePeaks(1,firstSweepWithSpikes);
    
%% Spike times for sweeps to plot
spikePeaksLocs1 = spikePeaksLocs(:,1); %spike locations for 1st sweep
spikePeaksLocs1(spikePeaksLocs1==0) = []; %remove zeros
spikePeaks1 = spikePeaks(:,1); %spike peaks for 1st sweep
spikePeaks1(spikePeaks1==0) = []; %remove zeros
spikePeaksLocsLast = spikePeaksLocs(:,end); %spike locations for last sweep
spikePeaksLocsLast(spikePeaksLocsLast==0) = []; %remove zeros
spikePeaksLast = spikePeaks(:,end); %spike peaks for last sweep
spikePeaksLast(spikePeaksLast==0) = []; %remove zeros


%% Burst index and adaptation index

firingParamSweep = find(spikeN>=minSpikesFiringParam,1);

%Burst index
initialSpikesN = zeros(1,sweepN);
for i = 1:sweepN;
    initialSpikes_i = find(0 < spikeTimes(:,i) & spikeTimes(:,i)<(base+burstIndexWindow));
    initialSpikesN_i = length(initialSpikes_i);
    initialSpikesN(i) = initialSpikesN_i;
end

burstIndexCurve = initialSpikesN ./ spikeN;
burstFiringRate = initialSpikesN ./ (burstIndexWindow./1000);

spikeTimes(spikeTimes==0)=NaN;

%Adaptation index
isi = diff(spikeTimes);
first3isi = isi(1:3,:);

 %find location of the last spike and last 3 ISIs
lastSpikesLoc=zeros(sweepN,1);
for i=1:sweepN
    isi_tmp = isi(:,i);
    k = find(isnan(isi_tmp)); %find last NaN for each sweep (columns)
    lastSpikesLoc(i) = k(1);
end
lastSpikesLoc=lastSpikesLoc-1;

 %isolate last 3 ISIs
last3isi = zeros(3,sweepN);
for i=1:sweepN
    if lastSpikesLoc(i)==0; %skip sweeps with no spikes
        last3isi(1:3,i) = 0;
        i=i+1;
    elseif spikeN(i) < 7;
        last3isi(1:3,i) = 0;
    else
       last3isi(1:3,i) = isi(lastSpikesLoc(i)-2:lastSpikesLoc(i),i);
    end
end

adaptIndexCurve = mean(first3isi,1)./mean(last3isi,1);
rheobase=find(firingRate,1);

adaptIndex=adaptIndexCurve(firingParamSweep);
burstIndex=burstIndexCurve(firingParamSweep);


    
%% Analysis of the first spike
% Calculates spike threshold, half-width, amplitude and after-hyperpolarization

if singleSpikeParam == 1; %if flag to analyse first spike parameters is selected

    % Threshold
    firstSpikeTrace = Ytrace(firstSpikePeakLoc-(firstSpikeWindow/2*Fs):firstSpikePeakLoc+(firstSpikeWindow/2*Fs),firstSweepWithSpikes); %location of the first spike in the first sweep with spikes
    timeFirstSpikeTrace = ((0:(length(firstSpikeTrace))-1)./Fs)'; %time vector for first spike waveform
    firstSpikeTraceMinus1 = firstSpikeTrace(1:end-1); %spike waveform minus last point to plot against derivative
    diffFirstSpikeTrace = diff(firstSpikeTrace); %derivative of first spike
    spikeThresholdLoc = find(diffFirstSpikeTrace>2,1); %location of spike threshold (first point above threshold in the derivative). Default: 0.5
    spikeThresholdTime = spikeThresholdLoc/Fs; %time of spike threshold
    spikeThreshold = firstSpikeTraceMinus1(spikeThresholdLoc); %value of spike threshold

    % Amplitude
    firstSpikeAmpl = (firstSpikePeak-spikeThreshold); %amplitude of first spike

    % fAHP
    [~,tempPeakLoc]=findpeaks(firstSpikeTrace,'SORTSTR','descend','NPEAKS',1); %find relative time of spike in the trace
    fAHPcut = -firstSpikeTrace(tempPeakLoc+1:tempPeakLoc+fAHPwindow*Fs); %isolate fragment for fAHP (from spike peak)
    [tempfAHPloc] = find(diff(-fAHPcut)>-1,1); %finds fAHP point based on first derivative
    fAHPloc = tempfAHPloc+tempPeakLoc; % location of fAHP in trace with only first spike
    fAHPtime = fAHPloc / Fs; %time of data point to measure fAHP
    fAHPvalue = firstSpikeTrace(fAHPloc); 
    fAHPampl = spikeThreshold-(fAHPvalue); %fAHP amplitude

    % Half-width
    firstSpikeHalfAmpl=firstSpikeAmpl/2; %first spike half amplitude
    
    firstSpikeTracePeakLoc = find(firstSpikeTrace==firstSpikePeak); %peak location in the spike waveform vector
    firstSpikeTraceUntilPeak = firstSpikeTrace(1:firstSpikeTracePeakLoc); %trace containing the first spike only until the spike peak (to calculate half amplitude in the ascending phase only)

        %Ascending phase
        [firstSpikeHalfAmpl1Delta firstSpikeTraceHalfAmpl1Loc] = min(abs(firstSpikeTraceUntilPeak-(spikeThreshold+firstSpikeHalfAmpl))); %point of the waveform corresponding to half amplitude
        firstSpikeTraceHalfAmpl1 = firstSpikeTraceUntilPeak(firstSpikeTraceHalfAmpl1Loc); %trace value corresponding to half amplitude in the ascending phase of the spike
        firstSpikeTraceHalfAmpl1Time = firstSpikeTraceHalfAmpl1Loc / Fs; %time of half amplitude data point (ascending phase)
        
        %Descending phase
        [firstSpikeHalfAmpl2Delta firstSpikeDescHalfAmpl2Loc] = min(abs(firstSpikeTrace(firstSpikeTracePeakLoc:end)-firstSpikeTraceHalfAmpl1)); %point of the waveform corresponding to half amplitude
        firstSpikeDescHalfAmpl2Loc = firstSpikeDescHalfAmpl2Loc-1; %take data point before in the trace as usually it is more accurate
        firstSpikeTraceHalfAmpl2 = firstSpikeTrace(firstSpikeTracePeakLoc+firstSpikeDescHalfAmpl2Loc); %trace value corresponding to half amplitude
        firstSpikeTraceHalfAmpl2Loc = (firstSpikeDescHalfAmpl2Loc + firstSpikeTracePeakLoc); %trace value corresponding to half amplitude in the descending phase of the spike
        firstSpikeTraceHalfAmpl2Time = (firstSpikeDescHalfAmpl2Loc + firstSpikeTracePeakLoc) / Fs; %time of half amplitude data point (descending phase)
        spikeHalfWidth = firstSpikeTraceHalfAmpl2Time - firstSpikeTraceHalfAmpl1Time; %spike half-width in ms
end

if plotting == 1; %if plotting is set on 'yes'

%% Plotting
    %close all;
    figure('Position', [100, 100, 1049, 895]);

    %first sweep
    subplot(3,3,1);
    plot(time,Ytrace(:,1));
    hold;
    xlim([0 max(time)]);
    scatter(time(spikePeaksLocs1),spikePeaks1,50,'r');
    ylabel('Membrane potential (mV)');
    xlabel('Time (ms)');
    title('Spike detection 1st sweep');
 
    %last sweep
    subplot(3,3,2);
    plot(time,Ytrace(:,end));
    hold;
    xlim([0 max(time)]);
    scatter(time(spikePeaksLocsLast),spikePeaksLast,50,'r');
    ylabel('Membrane potential (mV)');
    xlabel('Time (ms)');
    title('Spike detection last sweep');

    %raster
    subplot(3,3,3);
    hold;
    xlim([0 max(time)]);
    ylim([0 sweepN]);
    for i = 1:sweepN;
    subplot(3,3,3);
    plot (spikeTimes(:,i), i,'*b');
    end
    ylabel('Sweep n');
    xlabel('Time');
    title('Raster plot');

    %f/I curve
    subplot(3,3,4);
    plot (current,firingRate,'-ob');
    title('f/I curve')
    ylabel('Firing rate (Spikes/s)')
    xlabel('Injected current (pA)')
    xlim([0 max(current)]);

    %f/I curve first 50 ms
    subplot(3,3,5);
    plot (current,burstFiringRate,'-or');
    title('f/I curve for first 50 ms')
    ylabel('Firing rate (Spikes/s)')
    xlabel('Injected current (pA)')
    xlim([0 max(current)]);

    %burst index
    subplot(3,3,6);
    plot (current,burstIndexCurve,'-or');
    title('Burst index')
    ylabel('Fraction of spikes in the first 50 ms')
    xlabel('Injected current (pA)')
    xlim([0 max(current)]);
    
     %adaptation index
    subplot(3,3,7);
    plot (current,adaptIndexCurve,'-or');
    title('Adaptation index')
    ylabel('Adaptation index')
    xlabel('Injected current (pA)')
    xlim([0 max(current)]);

    if singleSpikeParam == 1;
    %First spike
        %Half-width, threshold and fAHP
        subplot(3,3,8);
        plot(timeFirstSpikeTrace,firstSpikeTrace)
        xlim([min(timeFirstSpikeTrace) max(timeFirstSpikeTrace)]);
        hold on;
        scatter(spikeThresholdTime,spikeThreshold,100,'r');
        scatter(firstSpikeTraceHalfAmpl1Time,firstSpikeTraceHalfAmpl1,100,'r');
        scatter(firstSpikeTraceHalfAmpl2Time,firstSpikeTraceHalfAmpl2,100,'r');
        scatter(fAHPtime,fAHPvalue,100,'r');
        title('First spike')
        ylabel('Membrane potential (mV)')
        xlabel('time (ms)')
        hold off;

        %dV/V
        subplot(3,3,9);
        plot(firstSpikeTraceMinus1,diffFirstSpikeTrace);
        title('First spike dV/V')
        ylabel('Membrane potential (mV)')
        xlabel('dV')
        %xlim([diffFirstSpikeTrace(1) max(diffFirstSpikeTrace)]);
        %ylim([firstSpikeTrace(1) max(firstSpikeTrace)]);
      
    end
end
        
end

  

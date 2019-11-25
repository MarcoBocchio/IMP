function [inputRes, membraneTau, membraneCapacitance, spikeThreshold, firstSpikeAmpl, spikeHalfWidth, fAHPAmpl, fileName] = iv_old(fileFormat);
%%
%OPERATION
% Calculates membrane and single spike parameters from I/V protocol
%
% >>> INPUT VARIABLES >>>
% NAME              TYPE, DEFAULT      DESCRIPTION
% fileFormat        scalar, 1          file format to be loaded: 1 for ABF, 2 for
%                                      HEKA (.mat)
%%
% >>> PARAMETERS >>>
% NAME              DEFAULT            DESCRIPTION
% base              100                baseline in ms
% firstCurrAmpl     -50                amplitude of first current injection
% currInjDuration   500                length of current injection (in ms)
% currInjDelta      10                 increasing current amplitude
% voltageRange      10                 window that is averaged to measure membrane potential at the end of each current injection (in ms)
% voltageEndTime    50                 time before the end of the current
%                                      step in which membrane potential is measured (in ms)
% tauSweepN         3                  sweep with -30pA current injection to calculate tau
% tauStartCut       2                  time to cut off from beginning of tau trace to
%                                      ensure a good fit (in ms)
% firstSpikeWindow  30                 time window to quantify and plot the first spike in the IV (in ms)
% SpikeThreshold    -20                Threshold to detect spikes (value in Y
%                                      trace, not amplitude)
% plotting          1                  Plotting on/off

%%
% >>> ADDITIONAL FUNCTIONS REQUIRED >>>
% - Ezyfit

%%
% Marco Bocchio, updated 17/4/2019

%% Loading files

%ABF file (default)
if nargin == 0;
    fileFormat = 1;
end

if fileFormat == 1;
    [Ytrace,time, Fs, si_ms, fileName] = import_abf;
else
    
%HEKA file
   [YtraceCh1,YtraceCh2,time, Fs, si_ms, fileName] = import_heka;
   Ytrace = YtraceCh1; %use trace from channel 2
end

%% Parameters
plotting = 1; %select 1 yes and 0 for no
conversion = 0; % if set on 1, the trace will be converted to mV (multiplied by 1000)
%si_ms = si.si./1000; %sampling interval converted to ms (from µs)
%Fs = si.si *0.1; %sampling rate (in kHz), obtained from sampling interval in abfload
base = 100; % baseline time (in ms)
firstCurrAmpl = -50;
currInjDelta = 10;  % delta by which current step amplitude is increased
currInjDuration = 500; %length of depolarizing/hyperpolarizing current injections (in ms)
voltageRange = 10;  %range that is averaged to measure membrane potential at the end of each current injection (in ms)
voltageEndTime = 50; %ms before the end of the current step in which membrane potential is measured (usually set to 50 ms)
voltageEndPnts = voltageEndTime*Fs;
SpikeThreshold = -20; % spike threshold (trace value, not amplitude)
tauSweepN = 3; %sweep with -30pA current injection to calculate tau
tauSweepLength = 40; %length of the sweep to use to fit exponential for tau (in ms)
tauStartCut = 2;
sweepTime = (0:length(Ytrace)-1).*si_ms; %time of each sweep (in ms)
firstSpikeWindow = 30; %time window that is taken to quantify and plot the first spike in the IV (in ms)


%% Conversion to mV (if data are in V)
if conversion == 0
else
    Ytrace = Ytrace*1000;
end


%% IV for Input resistance
sweepsWithSpikes = zeros(1,size(Ytrace,2)); %preallocation

    %% Spike detection
    spikePeaks = zeros(200,size(Ytrace,2)); %preallocation
    spikePeaksLocs = zeros(200,size(Ytrace,2)); %preallocation
    
    for i = 1 : size(Ytrace,2);
       [spikePeaks_i,spikePeaksLocs_i] = findpeaks(Ytrace(:,i),'MinPeakHeight',SpikeThreshold); %spike peak detection in every sweep
       noSpikes = isempty(spikePeaks_i); %checks if the sweep has no spikes
       if noSpikes == 0;  %if there are spikes
          sweepsWithSpikes(i) = 1;     %value of vector corresponding to sweep is one
          spikePeaks(1:length(spikePeaks_i),i) = spikePeaks_i; %allocation of spike peaks values to matrix (each sweep is a column)
          spikePeaksLocs(1:length(spikePeaksLocs_i),i) = spikePeaksLocs_i;%allocation of spike peak locations values to matrix (each sweep is a column)
       else
           
       end
            
    end

       firstSweepWithSpikes = find(sweepsWithSpikes,1); %first sweep with spikes
       firstSpikePeakLoc = spikePeaksLocs(1,firstSweepWithSpikes); %location of the first spike in the first sweep with spikes
       firstSpikePeak = spikePeaks(1,firstSweepWithSpikes);
       
           
 %% Matrices
if sweepsWithSpikes == 0;
    inputfileIV = Ytrace;
else
    inputfileIV = Ytrace(:,1:firstSweepWithSpikes-1);
end
    
basePnts = base * Fs;
currPnts = currInjDuration * Fs;
avgPnts = voltageRange * Fs;
avgStartPnts = currPnts - voltageEndPnts; %point in the trace where measurement of membrane potential starts
[dataPnts,sweeps] = size(inputfileIV);
baseV = mean(inputfileIV(1:basePnts,:)); %baseline membrane potentials
V = mean(inputfileIV(avgStartPnts:avgStartPnts+avgPnts,:)); %membrane potentials at the end of each current step
deltaV = V - baseV; %vector with delta membrane potentials
current = (firstCurrAmpl:currInjDelta:sweeps*currInjDelta+firstCurrAmpl-currInjDelta); %vector with current steps
deltaV(find(current==0)) = 0;
inputRes = deltaV/current*1000; %input resistance in M? (slope of the I/V fit with intercept forced through zero)

    
%% tau

%only parameters (calculation of tau done with Ezyfit function in plotting
%section
tauSweep = inputfileIV(:,tauSweepN);%create vector with sweep to calculate tau (usually -20pA)
startCut=(tauStartCut+base)*Fs; %empirically works
endCut=(base+tauStartCut+tauSweepLength)*Fs;
tauSweepCut = tauSweep(startCut:endCut); %sweep to calculate tau limited to current injection
timeTauSweepCut = [1:length(tauSweepCut)]'./Fs; %time vector

if sweepsWithSpikes == 0;
    warning('WARNING: NO SPIKES DETECTED');
    spikeThreshold = 0;
    firstSpikeAmpl = 0;
    spikeHalfWidth = 0;
    fAHPAmpl = 0;
else

%% Analysis of the first spike
% Calculates spike threshold, half-width, amplitude and after-hyperpolarization

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
    firstSpikeTraceNorm=firstSpikeTrace-mean(firstSpikeTrace(1:100));
    [fAHPvalue,fAHPloc] = findpeaks(-firstSpikeTraceNorm,'MinPeakHeight',10, 'MinPeakDistance', minSpikeDetPnts); %finds fAHP from normalised negative trace
    fAHPvalue = firstSpikeTrace(fAHPloc); %real fAHP value (from original trace)
    fAHPampl = spikeThreshold-fAHPvalue; %fAHP amplitude
    fAHPtime = fAHPloc / Fs; %time of data point to measure fAHP


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

   if plotting == 1;     
        
%% Plotting
    %close all
    figure('Position', [100, 100, 1024, 640]);


    %% IV
    subplot(3,2,2)
    plot(sweepTime,inputfileIV)
    xlim([0 max(sweepTime)])
    title('Sweeps used for IV')
    ylabel('Membrane potential (mV)')
    xlabel('time (ms)')


    subplot(3,2,3); 
    scatter (current,deltaV)
    title('IV')
    ylabel('Depolariz/Hyperpolariz amplitude (mV)')
    xlabel('Current (pA)')

    %% tau
    subplot(3,2,4);
    plot(timeTauSweepCut, tauSweepCut,'k')
    title('Membrane tau fit')
    ylabel('Membrane potential (mV)')
    xlabel('time (ms)')
    f=showfit('f(t) = a+exp(-t/tau1)+exp(-t/tau2)','fitlinewidth',2,'fitcolor','red');
    membraneTau=max(f.m);

    %% Membrane capacitance
    membraneCapacitance = membraneTau/inputRes*1000; %membrane capacitance in pF
                  
            
    if sweepsWithSpikes == 0;
    else

%% Spikes
            subplot(3,2,4);
            plot(sweepTime,Ytrace(:,firstSweepWithSpikes));
            hold;
            scatter(sweepTime(firstSpikePeakLoc),firstSpikePeak,50,'r')
            xlim([0 max(sweepTime)]);
            title('First spike detection')

            %% First spike
            subplot(3,2,5);
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


            subplot(3,2,6);
            plot(firstSpikeTraceMinus1,diffFirstSpikeTrace);
            title('First spike dV/V')
            ylabel('Membrane potential (mV)')
            xlabel('dV')
            %xlim([diffFirstSpikeTrace(1) max(diffFirstSpikeTrace)]);
            %ylim([firstSpikeTrace(1) max(firstSpikeTrace)]);
    end
            
   else
       
   end

%% Printing results
inputResDisplay = ['Rin = ', num2str(inputRes), ' MOhm'];
disp (inputResDisplay);
membraneTauDisplay = ['Membrane tau = ', num2str(membraneTau), ' ms'];
disp (membraneTauDisplay);
membraneCapacitanceDisplay = ['Membrane capacitance = ', num2str(membraneCapacitance), ' pF'];
disp (membraneCapacitanceDisplay);
spikeThresholdDisplay = ['Spike threshold = ', num2str(spikeThreshold), ' mV'];
disp (spikeThresholdDisplay);
spikeAmplDisplay = ['Spike amplitude = ', num2str(firstSpikeAmpl), ' mV'];
disp (spikeAmplDisplay);
spikeHalfWidthDisplay = ['Spike half-width = ', num2str(spikeHalfWidth), ' ms'];
disp (spikeHalfWidthDisplay);
fAHPAmplDisplay = ['After-hyperpolarization amplitude = ', num2str(fAHPAmpl), ' mV'];
disp (fAHPAmplDisplay);
        
end


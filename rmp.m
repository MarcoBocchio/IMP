function [rmpValue] = rmp(fileFormat);
%%
%OPERATION
% Calculates average resting membrane potential from a continuous trace.
% Spike detection is performed and if firing rate is higher than 2Hz, the
% measurement is discarded

%%
% >>> INPUT VARIABLES >>>
% NAME        TYPE, DEFAULT      DESCRIPTION
% fileFormat  scalar, 1          file format to be loaded: 1 for ABF, 2 for
%                                HEKA (.mat)
%%
% >>> PARAMETERS >>>
% NAME              DEFAULT      DESCRIPTION
% SpikeThreshold    -10          Threshold to detect spikes (value in Y
%                                trace, not amplitude)
% plotting          1            Trace is plotted

%% Loading files

%ABF file (default)
if nargin == 0;
    fileFormat = 1;
end

if fileFormat == 1;
    [Ytrace,time, Fs] = import_abf;
else
    
%HEKA file
   [Ytrace1,~,time, Fs] = import_heka;
   Ytrace = Ytrace1; %use trace from channel 2
end

%% Parameters
plotting = 1; %select 1 yes and 0 for no
SpikeThreshold = -10; %membrane potential value, not amplitude


%% Spike detection
[spikePeaks,spikePeaksLocs] = findpeaks(Ytrace,'MinPeakHeight',SpikeThreshold);
traceDuration = length(Ytrace)/Fs; % trace duration in ms
noSpikes = isempty(spikePeaks); %checks if the trace has no spikes
       if noSpikes == 0  %if there are spikes
        spikeN = length(spikePeaks); %calculate number of spikes in the trace
        firingRate = spikeN / traceDuration * 1000; % and firing rate: spikes/s (Hz)
            if firingRate > 2 %if firing rate is > 2 spikes/s
                warning('Cell is spontaneously active and fires > 2 Hz. Reject resting membrane potential measurement')
            end
       else
           
       end


%% Resting membrane potential
rmpValue = mean(Ytrace);

%% Printing results
rmpDisplay = ['Resting membrane potential = ', num2str(rmpValue)];
disp (rmpDisplay);

if plotting == 1; %if plotting is set on 'yes'
    
    %% Plotting
    %close all;
    figure;
    plot (time,Ytrace)
    title('Resting membrane potential')
    ylabel('Membrane potential (mV)')
    xlabel('Time (s)')

else
    
end

end
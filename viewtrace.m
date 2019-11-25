function viewtrace(fileFormat, varargin); 

%% 
%OPERATION
% Plots traces from ABF files or HEKA files (exported in .mat format)

%%
% >>> INPUT VARIABLES >>>
% NAME        TYPE, DEFAULT      DESCRIPTION
% fileFormat  scalar, 1          file format to be loaded: 1 for ABF, 2 for
%                                HEKA (.mat)
% sweepN      scalar/array       number of sweep(s) to be plotted

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


if length(varargin) == 0;
    sweepN = size(Ytrace,2);
else
sweepN = cell2mat(varargin); % creates vector with number of sweeps  to be displayed
end


%% Plotting sweeps
close all;
figure;

if length(varargin) == 0;
    plot(time,Ytrace);
else
plot(time,Ytrace(:,sweepN));
end
xlabel('Time (ms)');
xlim([0 max(time)]);

end
function [varargout] = spl(varargin)
% [Leq, LAeq] = spl(filename)
% [Leq, LAeq] = spl(data, fs)
% [Leq, LAeq] = spl(data, t, fs)
%
% [Leq, LAeq, LCeq, LAmax, LCpeak] = spl(...)
%
% [Leq, LAeq, LCeq, LAmax, LCpeak, L, LA, LC, t] = spl(...) where L, LA and
% LC are the linear, A-weighted, and C-weighetd levels depending on the
% timefilter ('fast','slow') and t is the time vector for the levels.
%
% [Leq, LAeq, LCeq, LAmax, LCpeak, Loct, LAoct, Fc] = spl(...,'bandwidth', bw)
% performs octave-band analysis according to bandwidth bw = '1/3 octave' or
% '1 octave'. The sound presure levels are the result after the octave-band
% filtering.
%
% [output] = spl(...,'outputformat','struct') the output will be in a sruct
%
%
% [...] = spl(...,'var',value)
% var : 'channel','calib','p0', 'timefilter', 'bandwidth', 'filterorder'
%


if nargin<1 || isempty(varargin{1})
    [file,path] = uigetfile('*.wav');
    filename = fullfile(path,file);
    
    [data,fs] = audioread(filename);
    t = [];
    
    narg = 2;
elseif isstruct(varargin{1})
    t = [];
    fnames = fieldnames(varargin{1});
    for k = 1:length(fnames)
        fn = fnames{k};
        switch validatestring(fn,{'data','fs','time','t'})
            case 'data'
                data = varargin{1}.(fn);
            case {'time','t'}
                t = varargin{1}.(fn);
            case 'fs'
                fs = varargin{1}.(fn);
        end
    end
    narg = 2;
else
    data = varargin{1};
    if length(varargin{2})>1
        narginchk(3,nargin)
        t = varargin{2};
        fs = varargin{3};
        narg = 4;
    else
        narginchk(2,nargin)
        fs = varargin{2};
        t = [];
        narg = 3;
    end
end

valid_inputs = {...
    'channel', 'calib', 'p0', 'timefilter', ...
    'timelimits', 'limits', ...
    'bandwidth', 'filterorder', ...
    'outputformat'};

% Defaults values
outputformat = 'variables';
channel = 1:size(data,2);
calib = 1;
timelimits = [];
p0 = 2e-5;
timefilter = 'fast';
bandwidth = 'full';  % Bandwidth
filterorder = 6;  % Filter order
if isempty(t)
    t = (0:length(data)-1)'./fs;
end

assignin('base','argin',varargin);
args = {};
if nargin>=narg
    if isstruct(varargin{narg})
        args = namedargs2cell(varargin{narg});
    else
        args = varargin(narg:end);
    end
end

for n = 1:2:length(args)
    switch validatestring(args{n},valid_inputs)
        case 'channel'
            channel = args{n+1};
        case 'calib'
            calib = args{n+1};
        case 'p0'
            p0 = args{n+1};
        case 'timefilter'
            timefilter = validatestring(args{n+1},{'fast','slow'});
        case 'bandwidth'
            if isnumeric(args{n+1})
                bandwidth = args{n+1};
            else
                switch validatestring(args{n+1},{'full','1/3 octave','1 octave','octave'})
                    case 'full'
                        bandwidth = 'full';
                    case {'1/3 octave'}
                        bandwidth = '1/3 octave';
                    case {'1 octave','octave'}
                        bandwidth = '1 octave';
                end
            end
        case 'filterorder'
            filterorder = args{n+1};
        case 'timelimits'
            timelimits = args{n+1};
        case 'limits'
            limits = args{n+1};
            timelimits = t(end)./100.*limits;
        case 'outputformat'
            outputformat = validatestring(args{n+1},{'variables','struct'});
            if strcmp(outputformat,'struct')
                nargoutchk(0,1)
            end
    end
end

if ~isempty(timelimits)
    timelimits(1) = max(timelimits(1),t(1));
    timelimits(2) = min(timelimits(2),t(end));
    i = t>=timelimits(1) & t<=timelimits(2);
else
    i = true(length(data),1);
end

y = calib.*data(i,channel);
y = y-mean(y);
t = t(i);

% A and C-weighting filters
AweightFilt = weightingFilter('A-weighting',fs);
CweightFilt = weightingFilter('C-weighting',fs);

if isnumeric(bandwidth)
    Fc = bandwidth;
    if ~isfinite(Fc(1))
        bandwidth = 'lowpass';
        [B,A] = butter(filterorder, Fc(2)/(fs/2), 'low');
    elseif ~isfinite(Fc(2))
        bandwidth = 'highpass';
        [B,A] = butter(filterorder, Fc(1)/(fs/2), 'high');
    else
        bandwidth = 'bandpass';
        [B,A] = butter(filterorder, Fc./(fs/2));
    end
    y = filter(B,A,y);
end

switch timefilter
    case 'fast'
        frame_duration = 125e-3;
    case 'slow'
        frame_duration = 1;
end
frame_length = fs * frame_duration; % Number of samples per frame
[frames, num_frames] = filter_time(y, frame_length);
t_filter = (t(1):frame_duration:t(end))';
t_frames = t_filter(1:num_frames);

switch bandwidth
    case {'full','lowpass','highpass','bandpass'}
        [Lp, LA, LC] = deal(zeros(num_frames, size(y,2)));
        
        % Leq
        Leq = 10*log10(var(y)./p0^2);
        Lp(:,:) = 10*log10(var(frames)./p0^2);
        
        % A-weighted
        yA = AweightFilt(y);
        LAeq = 10*log10(var(yA)./p0^2);       
        LA(:,:) = 10*log10(var(filter_time(yA, frame_length))./p0^2);
        LAmax = max(LA);
        
        
        % C-weighted
        yC = CweightFilt(y);
        LCeq = 10*log10(var(yC)./p0^2);
        LCpeak = 10*log10(max(yC.^2)./p0^2);
        LC(:,:) = 10*log10(var(filter_time(yC, frame_length))./p0^2);
                
        output = cell2struct(...
            {Leq, LAeq, LCeq, LAmax, LCpeak, Lp, LA, LC, t_frames}, ...
            {'Leq', 'LAeq', 'LCeq', 'LAmax', 'LCpeak', 'L', 'LA', 'LC', 't'},...
            2);
        
        if strcmp(outputformat,'variables')
            varargout = struct2cell(output);
        else
            varargout = {output};
        end
        
        
    otherwise        
        % octave
        octFiltBank = octaveFilterBank(bandwidth,fs,'FilterOrder',filterorder);
        Fc = octFiltBank.getCenterFrequencies;
        
        [Leq, LAeq, LCeq] = deal(zeros(1,size(y,2)));
        [Loct, LAoct] = deal(zeros(length(Fc), size(y,2)));
        LA = zeros(num_frames, size(y,2));
                
        yo = octFiltBank(y);
        Loct(:,:) = 10*log10(var(yo)./p0^2);
        Leq(:,:) = 10*log10(sum(var(yo)./p0^2));
        
        yAo = octFiltBank(AweightFilt(y));
        LAoct(:,:) = 10*log10(var(yAo)./p0^2);
        LAeq(:,:) = 10*log10(sum(var(yAo))./p0^2);
        
        [B,A] = butter(filterorder,octFiltBank.FrequencyRange/(fs/2));
        yf = filter(B,A,y);
        
        LA(:,:) = 10*log10(var(filter_time(AweightFilt(yf), frame_length))./p0^2);
        LAmax = max(LA);
        
        yCo = octFiltBank(CweightFilt(y));
        LCeq(:,:) = 10*log10(sum(var(yCo))./p0^2);
        LCpeak = 10*log10(max(CweightFilt(yf).^2)./p0^2);
        
        output = cell2struct(...
            {Leq, LAeq, LCeq, LAmax, LCpeak, Loct, LAoct, Fc(:)}, ...
            {'Leq', 'LAeq', 'LCeq', 'LAmax', 'LCpeak', 'Loct', 'LAoct', 'Fc'},...
            2);
        
        if strcmp(outputformat,'variables')
            varargout = struct2cell(output);
        else
            varargout = {output};
        end
        
end

function [frames, num_frames] = filter_time(y, frame_length)
% Split the signal into frames
num_frames = floor(size(y,1) / frame_length);
if size(y,2)>1
    frames = zeros(frame_length, num_frames, size(y,2));
else
    frames = zeros(frame_length, num_frames);
end

for i = 1:num_frames
    start_sample = (i - 1) * frame_length + 1;
    end_sample = min(i * frame_length, size(y,1));
    if size(y,2)>1
        frames(1:(end_sample - start_sample + 1),i,:) = y(start_sample:end_sample,:);
    else
        frames(1:(end_sample - start_sample + 1),i) = y(start_sample:end_sample);
    end
end
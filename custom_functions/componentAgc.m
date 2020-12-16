function c = componentAgc(c,comp,varargin)

% This function applies automatic gain control (AGC) to each trace. This 
% process, commonly used in seismic reflection processing, applies a
% variable scale to each trace such that the amplitude along the entire
% trace is roughly uniform. By minimizing amplitude variations along the
% trace, well-correlated but low-amplitude phases become more visible. A
% time window may be specified to control how tightly this scaling is
% applied. Use a longer window for lower frequecy signals. The function
% returns a correlation object.
%
% EXAMPLES:
%  c=agc(c)    				apply agc using the default time window (0.5 s)
%  c=agc(c,0.8)            	apply agc using a window of 0.8 s

% Author: Michael West, Geophysical Institute, Univ. of Alaska Fairbanks
% $Date$
% $Revision$


if (length(varargin) >= 2)
    error('Too many inputs');
end


if length(varargin)==1
	agcwin = varargin{1};
else
	agcwin = 0.5;
end


ncc = length(comp);
cstation = fieldnames(c.(comp{1}));

for s = 1:1:length(cstation)
    if ~isa(c.(comp{1}).(cstation{s}).corr,'correlation')
        error('First input must be a correlation object');
    end
    
    %get traces for same record across components
    wav = get(c.(comp{1}).(cstation{s}).corr,'waveform');
    nw = length(wav);
    wav = repmat(waveform(),nw, ncc);
    
    for nc = 1:1:ncc
        wav(:,nc) = get(c.(comp{nc}).(cstation{s}).corr,'waveform');
    end
    
    % LOOP THROUGH TRACES APPLYING GAIN
    for tracenum = 1:1:nw        
        agcsamp = round( agcwin * get(wav(1,1),'Fs') );

        %w = get(wav(tracenum,:),'DATA');

        for nc = 1:1:ncc
            w(:,nc) = get(wav(tracenum,nc),'DATA');
            % correct for average amplitude on each channel
            w(:,nc) = w(:,nc) ./ nanmean(abs(w(:,nc)));
        end
        
        
        scale = zeros( size(w,1)-2*agcsamp , 1 );
        for index = -1*agcsamp:agcsamp
            % contains abs of all channels
            win = abs( w(agcsamp+index+1:agcsamp+index+length(scale),:));
            % average over channels
            scale = scale + nanmean(win,2);
        end
        
        %for index=-1*agcsamp:agcsamp
        %    scale=scale + abs( w(agcsamp+index+1:agcsamp+index+length(scale)) );
        %end;
        
        %average across window
        scale = scale/nanmean(abs(scale));
        scale = [ones(agcsamp,1)*scale(1) ; scale ; ones(agcsamp,1)*...
            scale(end)];
        for nc = 1:1:ncc
            w(:,nc) = w(:,nc)./scale;
            wav(tracenum,nc) = set(wav(tracenum,nc),'DATA', w(:,nc));
        end
    end
    
    for nc = 1:1:ncc
        c.(comp{nc}).(cstation{s}).corr = set(c.(comp{nc}).(cstation{s}...
            ).corr,'waveform', wav(:,nc));
    end
end

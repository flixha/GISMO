function c2 = removeEmptyTraces(c2)

% removeEmptyTraces removes empty traces from a network correlation object.
%   Such traces may appear during processing if the raw trace contained
%   e.g. a NaN value.
% Input:
%   c2        : network-correlation object
% Output:
%   c2        : network-correlation object without empty traces

% remove empty / zero traces
comp = fieldnames(c2);
cstation = fieldnames(c2.(comp{1}));

for k=1:1:length(cstation)
    nw = length(get(c2.Z.(cstation{k}).corr,'waveforms'));
    % for all comps except Three-comp-obj
    wav = repmat(waveform(),nw, length(comp)-1);
    for nc = 1:1:length(comp)-1
        wav(:,nc) = get(c2.(comp{nc}).(cstation{k}).corr,'waveforms');
    end

    % Check whether any of the waveform of one event is all zero
    keepWavIdx = true(nw,1);
    for j=1:1:nw
        if any(( nanmax(wav(j,:)) == 0 & nanmin(wav(j,:)) == 0 ) |...
                (isnan(max(wav(j,:))) & isnan(min(wav(j,:))) ) )
            keepWavIdx(j) = false;
        end
    end

    % only keep the ones where nothing is zero / NaN
    for nc = 1:1:length(comp)
        if isa(c2.(comp{nc}).(cstation{k}),'threecomp')
            tc = c2.(comp{nc}).(cstation{k});
            c2.(comp{nc}).(cstation{k}) = tc(keepWavIdx);
        else
            c2.(comp{nc}).(cstation{k}).corr = subset(c2.(comp{nc}).(...
                cstation{k}).corr, keepWavIdx);
            c2.(comp{nc}).(cstation{k}).cat = subset(c2.(comp{nc}).(...
                cstation{k}).cat, keepWavIdx);
        end
    end
end
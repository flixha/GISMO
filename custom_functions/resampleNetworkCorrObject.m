function cOut = resampleNetworkCorrObject(c, targetSamplingRate)

% resampleNetworkCorrObject resamples all waveforms in a network-
%   correlation object to a selected target sampling rate.
% Input:
%   c   : a network correlation object 
%   targetSamplingRate: sampling ratue to resample all traces to
% Output:
%   c2  : a network correlation object with resampled traces

    comp = fieldnames(c)';
    cstation = fieldnames(c.(comp{1}))';
    
    cOut = c;
    
    for k=1:1:numel(comp)
        for s=1:1:numel(cstation)
            %check that we're not resampling the three component objects
            if isa(c.(comp{k}).(cstation{s}),'threecomp')
                continue
            end
            if ~isa(c.(comp{k}).(cstation{s}).corr,'correlation')
                continue
            end
            
            wavs = get(c.(comp{k}).(cstation{s}).corr, 'waveform');

            
            %resample waveform to target Sampling rate
            crunchFactor = get(wavs,'freq') ./ targetSamplingRate;
            crunchFactor = round(crunchFactor, 4);

            % We can calculate the crunchFactor very accurately and 
            % efficiently with symbolic toolbox:
            if license('test', 'symbolic_toolbox')
                if length(sym(crunchFactor)) > 20
                    crunchFactor = round(crunchFactor, 2);
                end
                [Q, P] = numden(sym(crunchFactor));
            else
                % workaround so that we don't need symbolic toolbox
                crunchFactor = round(crunchFactor, 9, "significant");
                n_decimals = length(char(string( ...
                    crunchFactor - double(int64(crunchFactor))))) - 1;
                Q = int32(crunchFactor * 10 ^ n_decimals);
                P = int32(1 * 10 ^ n_decimals);
            end

            Q = double(Q);
            P = double(P);
            if all(Q==Q(1)) && all(P==P(1))
                Q = Q(1); P = P(1);
                D = double(wavs);
                
                
                
                ResampleD = resample(D,P,Q);
            end

            for w=1:1:numel(wavs) % can be safely parallelized with "parfor"
                wav = wavs(w);
                % put back into waveform, but don't forget to 
                % update the frequency               
                wav = set(wav,'data',ResampleD(:,w), 'Freq',...
                    targetSamplingRate);
                wavs(w) = wav;
            end
            cOut.(comp{k}).(cstation{s}).corr = set(...
                    cOut.(comp{k}).(cstation{s}).corr, 'waveform', wavs);
        end
    end

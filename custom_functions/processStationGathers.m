function c2 = processStationGathers(c, cstation, comp, minFrequency,...
    maxFrequency)

    c2 = c;

    for nc=1:1:length(comp)
        corr_matrices.(comp{nc}) = struct();
    end
    
    for k=1:1:length(cstation)
        for nc=1:1:length(comp)
            if isfield(c2.(comp{nc}), cstation{k})
                c2.(comp{nc}).(cstation{k}).corr = demean(...
                    c2.(comp{nc}).(cstation{k}).corr);
                c2.(comp{nc}).(cstation{k}).corr = detrend(...
                    c2.(comp{nc}).(cstation{k}).corr);
                
                % Fill NaN-value gaps
                waves = get(c2.(comp{nc}).(cstation{k}).corr,'waveform');
                waves = fillgaps(waves,0);
                c2.(comp{nc}).(cstation{k}).corr =...
                    set(c2.(comp{nc}).(cstation{k}).corr,'waveform',waves);
                
                c2.(comp{nc}).(cstation{k}).corr = crop(...
                    c2.(comp{nc}).(cstation{k}).corr,-4,30);
                c2.(comp{nc}).(cstation{k}).corr = taper(...
                    c2.(comp{nc}).(cstation{k}).corr);
                c2.(comp{nc}).(cstation{k}).corr = butter(...
                    c2.(comp{nc}).(cstation{k}).corr,[minFrequency 15]);
            end
        end
    end
    
    c2 = requestMinimumSNR(c2, 1.5);
    
    % ALIGNMENT
    for k=1:1:length(cstation)
        for nc=1:1:length(comp)
            if isfield(c2.(comp{nc}),cstation{k})   
                %trigger adjustments don't seem to make the whole traces align
                %better. if it was possible to define a max shift then it might
                %help. Source wavelets are just too different?

                %ALIGNMENT
                % only when we have Z-components
                if any(contains(comp,'Z'))
                    %only calculate the lags once; when running through
                    %Z-component loop
                    if strcmp(comp{nc},'Z')
                        %1 get lags from the hilbert transform of first
                        %arrivals
                        hil = hilbert(c2.Z.(cstation{k}).corr);
                        hil = xcorr(hil,[-0.2 0.5]);
                        hil = linkage(hil);
                        hilLags_Z = get(hil,'lag');

                        % calculate lags from all components ; but weight
                        % Z-component strongest
                        
                        hilLags = hilLags_Z;
                        
%                         hilLags = (hilLags_Z*2+hilLags_N+hilLags_E)/4;
                        %hil = adjusttrig(hil,'MMM', 0.2);
                    end
                    % 2. Correct all components with the same adjustments from
                    % hilbert transform
                    c2.(comp{nc}).(cstation{k}).corr = set(...
                        c2.(comp{nc}).(cstation{k}).corr, 'lag', hilLags);
                    c2.(comp{nc}).(cstation{k}).corr = adjusttrig(...
                        c2.(comp{nc}).(cstation{k}).corr, 'MMM', 0.2);
                end
                
                % Only filter again if the maximum frequency is lower
                % than during the first bandpass filtering
                if maxFrequency < 15
                    c2.(comp{nc}).(cstation{k}).corr = butter(...
                       c2.(comp{nc}).(cstation{k}).corr,...
                       [minFrequency maxFrequency]);
                else
                    disp('Will not filter again, maximum frequency too high')
                end
                c2.(comp{nc}).(cstation{k}).corr = xcorr(...
                    c2.(comp{nc}).(cstation{k}).corr,[-0.2 14]);
                
                % Fill linkage 
                if length(get(c2.(comp{nc}).(cstation{k}).corr,'trig'))>1
                    c2.(comp{nc}).(cstation{k}).corr = linkage(...
                        c2.(comp{nc}).(cstation{k}).corr);
                end
            end
        end
    end

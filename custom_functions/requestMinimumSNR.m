function c2 = requestMinimumSNR(c, minimumSNR)

% takes a network-object and removes all 3-c traces where no component has
% a signal-to-noise ratio below minimumSNR. 

comp = fieldnames(c);
cstation = fieldnames(c.Z)';

c2 = c;

number_removed = 0;
number_kept = 0;

for k=1:1:length(cstation)  
    if isfield(c.(comp{1}),cstation{k})
        if length(get(c.(comp{1}).(cstation{k}).corr,'waveform')) > 0
            noiseWZ = get(crop(c.(comp{1}).(cstation{k}).corr,-2.5,-0.5),...
                'waveform');
            signalWZ = get(crop(c.(comp{1}).(cstation{k}).corr,0,2),...
                'waveform');
            SNR_Z = mean(abs(signalWZ))./mean(abs(noiseWZ));

            noiseWN = get(crop(c.(comp{2}).(cstation{k}).corr,-2.5,-0.5),...
                'waveform');
            signalWN = get(crop(c.(comp{2}).(cstation{k}).corr,0,2),...
                'waveform');
            SNR_N = mean(abs(signalWN))./mean(abs(noiseWN));

            noiseWE = get(crop(c.(comp{3}).(cstation{k}).corr,-2.5,-0.5),...
                'waveform');
            signalWE = get(crop(c.(comp{3}).(cstation{k}).corr,0,2),...
                'waveform');
            SNR_E = mean(abs(signalWE))./mean(abs(noiseWE));

            goodSNR = find(SNR_Z > minimumSNR | SNR_N > minimumSNR |...
                SNR_E > minimumSNR);
        else
            goodSNR = [];
        end

        for nc=1:1:length(comp)
            c2.(comp{nc}).(cstation{k}).corr = subset(...
                c.(comp{nc}).(cstation{k}).corr, goodSNR);
            c2.(comp{nc}).(cstation{k}).cat = subset(...
                c.(comp{nc}).(cstation{k}).cat, goodSNR);
        end

        number_removed = number_removed + length(SNR_Z) - length(goodSNR);
        number_kept = number_kept + length(goodSNR);

        disp([num2str(round(k/length(cstation)*100)), ' % done, kept: ',...
            num2str(number_kept),', removed: ',num2str(number_removed)])
    end
end
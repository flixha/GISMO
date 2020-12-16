function [ax, SonRp_vs_PonZp_ord_smo, SonRp_vs_SonTp_ord_smo,...
    SonTp_vs_PonZp_ord_smo] = plotWavComparisonHistogram(c, sta, varargin)

% plotWavComparisonHistogram measures amplitude ratios and plots these in a
%   histogram to put alongside the waveform plots.

    % measure amplitude ratios / energy ratios:
    % P vs S
    % P vs 2-3 seconds after P    
    % P vs. time between P and S
    % first 2-3 seconds after P vs time until S
    % S vs. 2-3 seconds times 1.7 after S
    % S vs. time after S
    % first 2-3 seconds times 1.7 after S vs after that

    SECPERDAY = 24*60*60;
    prePhas = 0.25;

    p = inputParser;

    defaultMethod = 'DFST';
    validMethods = {'LAT', 'LATITUDE', 'LON', 'LONGITUDE', 'X',...
        'DISTFROMTRENCH', 'DISTANCEFROMTRENCH', 'Y', 'DISTALONGSTRIKE',...
        'Z', 'DEPTH', 'DISTFROMSLABTOP', 'DISTANCEFROMSLABTOP', 'SLABTOP',...
        'DFST', 'S-P','SMP','SMINUSP'};
    checkMethod = @(x) any(validatestring(x,validMethods));
    addParameter(p,'method', defaultMethod, checkMethod)

    defaultArrivals = struct();
    addParameter(p,'arrivals', defaultArrivals, @isstruct)
    
    defaultRunMedianWidth = 4;
    addParameter(p,'RunMedianWidth', defaultRunMedianWidth, @isnumeric)

    defaultWindowLength = 1.0;
    addParameter(p,'windowLength', defaultWindowLength, @isnumeric)
    
    defaultScaling = 'linear';
    validScaling = {'linear', 'lin', 'log', 'logarithmic'};
    checkScaling = @(x) any(validatestring(x,validScaling));
    addParameter(p,'scaling', defaultScaling, checkScaling);
    
    addRequired(p,'c',@isstruct);
    addRequired(p,'sta',@ischar);

    parse(p, c, sta, varargin{:})
    res = p.Results;
    
    method = res.method;
    arrivals = res.arrivals;
    winLen = res.windowLength;
    scaling = res.scaling;
    RunMedianWidth = res.RunMedianWidth;

    comp = 'Zp';
    

    isStackedCatalog = false;
    % Get trace Order , values etc.
    [traceOrder, sortedValues, titlename, ylabelname, stdValues,...
        invertplot] = getTraceOrderAndPlotMethodProps(c, method, comp,...
        sta, isStackedCatalog);
    nLines = numel(traceOrder);

    % Get arrival information
    if isempty(fields(arrivals))
        load('arrivals_0.02degGrid.mat','arrivals');
    end
    
    % initialize some variables for the whole gather
    pickOffset = zeros(nLines,1);
    SmPdifference = zeros(nLines,1);
    SmP = zeros(nLines,1);
    ratio_SonRp_vs_SonTp = zeros(nLines,1);
    
    MeanEnergy_Rp = nan(nLines,1);
    MeanEnergy_Tp = nan(nLines,1);
    Rpwav = repmat(waveform(),nLines,1);
    Tpwav = repmat(waveform(),nLines,1);
    
    MeanSenergy_Rp = nan(nLines,1);
    MeanSenergy_Tp = nan(nLines,1);
    S_Rpwav = repmat(waveform(),nLines,1);
    S_Tpwav = repmat(waveform(),nLines,1);
    
    MeanS_AmpRatio = nan(nLines,1);
    
    
    triggers = get(c.(comp).(sta).corr,'trig');
    
    q = 0;
    for j=traceOrder'
        %q = q+1;
        evID = c.(comp).(sta).cat.table.EventID(j);
        
        if arrivals(evID).EventID ~= evID
            warning('event IDs do not match')
            continue
        end
                           
        % get the right line for plotting
        [arrivalLine, arrIdx] = find(traceOrder==j);
        if invertplot
            arrivalLine = length(traceOrder) - arrivalLine + 1;
        end

        % check that theere are arrivals for that event
        if isempty(arrivals(evID).arrivals)
            continue
        end
        % check that there are arrivals for that station 
        if isempty(arrivals(evID).arrivals.(sta).time)
            continue
        end
        
        narrivals = numel(arrivals(evID).arrivals.(sta).time);

        % Two ways to find the first arrival: from the FM3D paths
        % or from the picks (triggers)
        firstArrival = getFirstArrivalTimeFromFM3D(c, arrivals,...
            sta, comp, traceOrder, evID, arrIdx);

        firstArrivalFromPicks = (triggers(j) - c.(comp).(sta...
            ).cat.table.otime(j)) * SECPERDAY;
        pickOffset(j) = firstArrival - firstArrivalFromPicks;

        
        
        % get_S_arrivalTime
        directIdx = find(arrivals(evID).arrivals.(sta).nConv == 0 &...
            arrivals(evID).arrivals.(sta).nReflections == 0);
        if length(directIdx) < 1; continue; end
        Sphase = subsetArrival(arrivals(evID).arrivals.(sta),...
            directIdx(end));
        % convert S-traveltime to relative time from P
        SmPtimeCalc = Sphase.time - firstArrivalFromPicks;
        
        % try to find picked S minus P-time
        SmPtimePicked = (c.(comp).(sta).cat.table.arrivals{j}.s_time -...
            c.(comp).(sta).cat.table.arrivals{j}.p_time) * SECPERDAY;
        
        if (contains(Sphase.phase,'P') || SmPtimeCalc < 0 || isempty(...
                Sphase)) && (SmPtimePicked < 0 || isnan(SmPtimePicked))
            continue
        end
        
        % Try to check calculated vs. picked difference
        SmPdiff = SmPtimeCalc - SmPtimePicked;
        if ~isempty(SmPdiff)
            SmPdifference(j) = SmPdiff;
        end
        
        if abs(SmPdiff) < winLen / 3
            SmP(j) = SmPtimePicked;
        else
            SmP(j) = SmPtimeCalc;
        end
    end
    
    % here we still have P-triggers
    ZpC = taper(crop(c.Zp.(sta).corr, [-0.5 2*max(SmP)]));
    ZpC = hilbert(demean(detrend(ZpC)));
    RpC = taper(crop(c.Rp.(sta).corr, [-0.5 2*max(SmP)]));
    RpC = hilbert(demean(detrend(RpC)));
    TpC = taper(crop(c.Tp.(sta).corr, [-0.5 2*max(SmP)]));
    TpC = hilbert(demean(detrend(TpC)));
    
    % P-wavelet
    P_ZpC = taper(crop(c.Zp.(sta).corr, [-winLen*prePhas, winLen*(1-prePhas)]));
    P_ZpC = hilbert(demean(detrend(P_ZpC)));
    
    S_trigger = get(c.Zp.(sta).corr,'trig') + SmP./SECPERDAY;
    c.Zp.(sta).corr = set(c.Zp.(sta).corr,'trig', S_trigger);
    c.Rp.(sta).corr = set(c.Rp.(sta).corr,'trig', S_trigger);
    c.Tp.(sta).corr = set(c.Tp.(sta).corr,'trig', S_trigger);
        
    % SV-wavelet
    S_RpC = taper(crop(c.Rp.(sta).corr, [-winLen*prePhas, winLen*(1-prePhas)]));
    S_RpC = hilbert(demean(detrend(S_RpC)));
    % SH-wavelet
    S_TpC = taper(crop(c.Tp.(sta).corr, [-winLen*prePhas, winLen*(1-prePhas)]));
    S_TpC = hilbert(demean(detrend(S_TpC)));
    
    
    ratio_SonRp_vs_SonTp = median(get(S_RpC,'waveform')) ./...
        median(get(S_TpC,'waveform'));
    
    ratio_SonRp_vs_PonZp = median(get(S_RpC,'waveform')) ./...
        median(get(P_ZpC,'waveform'));
    
    ratio_SonTp_vs_PonZp = median(get(S_TpC,'waveform')) ./...
        median(get(P_ZpC,'waveform'));
    
    
    % total S-energy vs total P-energy, normalized by direct P- and direct
    % S-arrivals

    Spost_RpC = taper(crop(c.Rp.(sta).corr, [-0.5, 2*max(SmP)]));
    Spost_RpC = hilbert(demean(detrend(Spost_RpC)));
    % SH-wavelet
    Stot_TpC = taper(crop(c.Tp.(sta).corr, [-0.5, 2*max(SmP)]));
    Stot_TpC = hilbert(demean(detrend(Stot_TpC)));
    
    Spre_RpC = taper(crop(c.Rp.(sta).corr, [-1*max(SmP) - 0.5, -0.5]));
    Spre_RpC = hilbert(demean(detrend(Spre_RpC)));
    
    % integrating the envelope means that the last value = total count, is 
    % the max
    Sdir_energy = max(integrate(get(S_RpC,'waveform')));
    Spost_energy = max(integrate(get(Spost_RpC, 'waveform')));
    Spre_energy = max(integrate(get(Spre_RpC, 'waveform')));
    
    ratio_refl_S_vs_conv_S =  (Spost_energy ./ Sdir_energy) ./...
        (Spre_energy ./ Sdir_energy);
    
    % smooth Out the Transverse vs. Radial-ratio of S
    SonRp_vs_SonTp_ordered = ratio_SonRp_vs_SonTp(traceOrder);
    SonRp_vs_PonZp_ordered = ratio_SonRp_vs_PonZp(traceOrder);
    SonTp_vs_PonZp_ordered = ratio_SonTp_vs_PonZp(traceOrder);
    nwin = RunMedianWidth;
    SonRp_vs_SonTp_ord_smo = SonRp_vs_SonTp_ordered;
    SonRp_vs_PonZp_ord_smo = SonRp_vs_PonZp_ordered;
    SonTp_vs_PonZp_ord_smo = SonTp_vs_PonZp_ordered;
    
    ratio_refl_S_vs_conv_S_ord = ratio_refl_S_vs_conv_S(traceOrder);
    ratio_refl_S_vs_conv_S_ord_smo = ratio_refl_S_vs_conv_S_ord;
    
    for j=1:1:nLines
        lower = j-nwin; if lower<1; lower = 1;end
        upper = j+nwin; if upper>nLines; upper = nLines;end
        SonRp_vs_SonTp_ord_smo(j) = nanmedian(...
            SonRp_vs_SonTp_ordered(lower:upper));
        SonRp_vs_PonZp_ord_smo(j) = nanmedian(...
            SonRp_vs_PonZp_ordered(lower:upper));
        SonTp_vs_PonZp_ord_smo(j) = nanmedian(...
            SonTp_vs_PonZp_ordered(lower:upper));
        
        ratio_refl_S_vs_conv_S_ord_smo(j) = nanmedian(...
            ratio_refl_S_vs_conv_S_ord(lower:upper));
    end
            
    
    figure
    set(gca,'Visible','off')
    set(gcf,'Position',[900 50 140 1000]);
    
    blue = [0    0.4470    0.7410];
    red =  [0.83    0.3 0];
    
    % Plot SV vs P ratio
    ax(1) = axes('XAxisLocation','bottom', 'YAxisLocation','left',...
        'Color','w');
    hold(ax(1),'on')
    plot(ax(1),SonRp_vs_PonZp_ord_smo, 1:1:nLines,'-','color',blue)
    ylim(ax(1),[0 nLines])
    
    
    if contains(scaling,'log')    
        upXLim = prctile(abs(log10(SonRp_vs_PonZp_ord_smo(nwin:...
            end-nwin))),98);
        if upXLim < 1
            upXLim = 1;
        end
        xlim(ax(1),[10^-upXLim 10^upXLim]);
        ax(1).XScale = 'log';
    else
        upX2Lim = prctile(SonRp_vs_PonZp_ord_smo(nwin:end-nwin), 99);
        xlim(ax(1),[0 upX2Lim])
    end
    ax(1).XGrid = 'on';
    ax(1).XMinorGrid = 'off';
    ax(1).XColor = blue;
    ax(1).YColor = 'k';
    ax(1).YTick = [];
    ax(1).XLabel.String = 'SV / P - Amp.';
    
    % Plot SH vs SV ratio
    ax(2) = axes('Position', ax(1).Position, 'XAxisLocation','top',...
        'YAxisLocation','right', 'Color','none');
    hold(ax(2),'on')
    plot(ax(2),SonRp_vs_SonTp_ord_smo, 1:1:nLines,'-','color',red)
    ylim(ax(2),[0 nLines])
    
    if contains(scaling,'log')     
        upXLim = prctile(abs(log10(SonRp_vs_SonTp_ord_smo(nwin:...
            end-nwin))),98);
        if upXLim < 1
            upXLim = 1;
        end
        xlim(ax(2),[10^-upXLim 10^upXLim]);
        ax(2).XScale = 'log';
    else
        xlim(ax(2),[0 prctile(SonRp_vs_SonTp_ord_smo(nwin:...
            end-nwin),99)]);
    end
    %ax(2).XGrid = 'on';
    
    ax(2).XColor = red;
    ax(2).YColor = 'k';
    ax(2).YTick = [];
    ax(2).XLabel.String = 'SV / SH - Amp.';
    
     


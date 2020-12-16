
%% This is the basic script to analyze seismograms for secondary arrivals.
%  The script creates a number of plots on which potential (chosen)
%  secondary arrivals are marked with a line and a label. Here the input
%  comes from previously saved mat-files that contain earthquake
%  hypocenters, picks, and waveforms readily in one GISMO-catalog object.
%  Alternatively, the input can be read from a wide range of databases with
%  the available GISMO functions. Note that some additional parameters
%  (x,y,z-coordinates in a local coordinate system) are added to the GISMO
%  catalog here.

%set parameters
reloadFM3Darrivals = false;
targetSamplingRateForProcessing = 100;
targetSamplingRateAfterFilter = 30;

% Choose stations to process + correlate waveforms for
cstation = {'S009','S010','S011','PE05','PE07','S012',...
    'PE02','S013','S014','S015','S016'};
corcomp = {'Z','E','N'};

minFrequency = 1.5;
maxFrequency = 10;

catalogFile = 'data/greekevents_468_2020_12_16.mat'; % Full Catalog
% Catalog file with waveforms needs to be split so each file is less than
% 100 MB
catalogFile_part1 = 'data/greekevents_468_part1.mat'; % Part 1
catalogFile_part2 = 'data/greekevents_468_part1.mat'; % Part 2
catalogFile_part3 = 'data/greekevents_468_part1.mat'; % Part 3
arrivalsFile = 'data/arrivals_0.02degGrid.mat';


%% Load in functions and data files
gismopath = pwd;
startup_GISMO(gismopath)
addpath([gismopath, '/custom_functions'])
addpath([gismopath, '/data'])

% load events and waveforms from file
% load(catalogFile, 'greekevents'); % If I have one file for all events
% Need to split the catalog in 3 to fit on GITHUB within 100 MB limit
load(catalogFile, 'gr1');
load(catalogFile, 'gr2');
load(catalogFile, 'gr3');
% Combine the files to one catalog again
greekevents = Catalog();
greekevents.table = [gr1.table; gr2.table; gr3.table];
greekevents.waveforms = [gr1.waveforms; gr2.waveforms; gr3.waveforms];


% Load in computed arrival times from an FastMarching FM3D run
% Load from previously saved mat file:
load(arrivalsFile, 'arrivals');
% Or load from original FM3D output files with function
% arrivals = loadFM3Darrivals(fm3Dpath, [1:1:greekevents.numberOfEvents]);

%% Select a set of earthquakes from the full catalog.
%  Here: select the active clusters of earthquakes in the mantle wedge, on
%  the interface, and in the subducting slab.

% Southern Tripoli cluster:
% Make area as wide as possible, but still rather vertical incoming rays.
maxY = 25;
minY = -35;
minZ = 20;
maxZ = 200;

% Get earthquakes within a region around the GRT cross section image
% Interface fault cluster
I_IFCluster = find(...
    greekevents.table.y > -15  & greekevents.table.y < -4.9 &...
    greekevents.table.x > 83.5 & greekevents.table.x < 93.5 &...
    greekevents.table.depth > 54 & greekevents.table.depth < 60 );
% Within-MW-cluster
I_MWCluster = find(...
    greekevents.table.y > -20 & greekevents.table.y < 0 &...
    greekevents.table.x > 69 & greekevents.table.x < 83 &...
    greekevents.table.depth < 56);
% Intraslabclsuter
I_ISCluster = find(...
    greekevents.table.y > -10 & greekevents.table.y < 0 &...
    greekevents.table.x > 65 & greekevents.table.x < 113 );
I_ISCluster = I_ISCluster(~ismember(I_ISCluster,I_IFCluster) &...
    ~ismember(I_ISCluster,I_MWCluster));

eventI = [I_MWCluster; I_IFCluster; I_ISCluster]';

% Other events below the array
I_ArrayCluster = find(...
    greekevents.table.x > 65 & greekevents.table.x < 125 &...
    greekevents.table.y > minY & greekevents.table.y < maxY &...
    greekevents.table.z > minZ & greekevents.table.z < maxZ);
I_ArrayCluster = I_ArrayCluster(~ismember(I_ArrayCluster,eventI));

eventI = [I_MWCluster; I_IFCluster; I_ISCluster; I_ArrayCluster]';
eventI = sort(eventI);
nSelectedEvents = length(eventI);


%% Collect all waveforms for the relevant events into a correlation object
% c is an object that contains GISMO correlation objects for the whole
% network, accessible as: c.(component).(stationname)
c = buildCorrelationCatalog(greekevents, eventI, cstation, corcomp,...
    targetSamplingRateForProcessing);

% Process the waveforms per station and channel to boost secondary arrival
% visibility
c2 = processStationGathers(c, cstation, corcomp,...
    minFrequency, maxFrequency);

%% Create threecomponent-objects, rotate and apply polarizationfilter

% First, throw out waveforms with bad Signal-to-noise-ratio
minimumSNR = 2.5;
c3 = requestMinimumSNR(c2, 2.5);

% Do polarization-filtering. Works well with these values:
% n, J, K, dt, width: (zrt, 1, 1, 2, 0.02, 0.6)
% the width is crucial! chose e.g. 0.5 s width for 1.5 - 10 Hz filter
windowLength = 10/maxFrequency * 0.5;

dt = 1/targetSamplingRateForProcessing;
width = 0.5;
c4 = threeComponentProcessing(c3, cstation, 0.5, 1, 2, dt, width,...
    true, false);
c4 = removeEmptyTraces(c4);

c5 = componentAgc(c4, {'Zp','Rp','Tp'}, 2);

%% Plot waveforms at a station sorted by some parameter of the event
%  (e.g., depth, lat, x-)
%  make nice printable plots for correlation objects

plotComp = {'Zp'};
printFigure = true;
plotEnvelope = true;
baseScale = 1.2;
xlims = [-2 16];
fileNameAddition0 = 'Whist_';

if plotEnvelope
    fileNameAddition0 = [fileNameAddition0, 'Envelope_'];
    %resample to reduce filesize
    cOut = resampleNetworkCorrObject(c5, targetSamplingRateAfterFilter);
end

% Plot the three-component figures for different stages of the processing: 
% (set e.g., p = [1,5] for raw plots and fully labeled and processed
% envelopes:
% p = 1: plot raw ZRT channels,
% p = 2: plot against properly scaled y-axis
% p = 3: plot with polarization-filtered ZRT channels,
% p = 4: plot polarization-filtered ZRT channels plus theoretical arrivals
% p = 5: plot polarization-filtered ZRT channels, plus arrivals and labels
for p=[1, 5]
    [plotArrivals, labelArrivals, plotComp, fileNameAddition0,...
        plotEnv, plotType, linewidth, scale, doPlotHistograms] =...
        getThreeCompPlotOptions(p, fileNameAddition0, plotEnvelope,...
        baseScale);
    if plotEnv
        scale = baseScale * 1.65;
    else
        scale = baseScale;
    end

    for j=1:1:length(cstation)
        clear ax;
        if doPlotHistograms
            [axx, SonRp_vs_PonZp, SonRp_vs_SonTp, SonTp_vs_PonZp] =...
                plotWavComparisonHistogram(c4, cstation{j},...
                'method','DFST','arrivals',arrivals,'windowLength',...
                1,'scaling','lin', 'RunMedianWidth', 4);
            if mean(SonRp_vs_PonZp) < 1
                plotPMS = true;
                avoidPtoS = false;
            end
        end


        for k=1:1:numel(plotComp)
            %try
            plotWavesAtStationSortedBy(cOut, plotComp{k},...
                cstation{j}, 'DFST', plotType,...
                'scale', scale, 'linewidth', linewidth,...
                'markarrivaltimes', plotArrivals, 'arrivals',...
                arrivals,'maxconversions',1,'maxreflections', 1,...
                'maxConvAndRefl', 2, 'maxPhaseLength', 5,...
                'labelarrivals',labelArrivals, 'plotEnvelope',...
                plotEnv, 'avoidPtoS', false, 'plotPMS',true,...
                'colorful',true)
            %catch
            %    continue
            %end

            oldPos=get(gcf,'Position');
            set(gcf,'Position',[1+(j-1)*230, oldPos(2), 350,...
                oldPos(4)])
            ax(k) = gca();
        end

        if doPlotHistograms
            ax = [ax, axx];
        end

        newax = formatThreeComponentWaveformFigures(ax,...
            cstation{j}, plotComp, plotType,...
            plotArrivals, plotEnvelope, labelArrivals,...
            fileNameAddition0, printFigure, xlims);

        if length(cstation) > 1
            close all;
        end
    end
end

set(gcf,'PaperPositionMode','Auto','PaperSize', [pos(3), pos(4)])
set(gcf,'renderer','painters')


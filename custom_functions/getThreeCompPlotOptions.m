function [plotArrivals, labelArrivals, plotComp, fileNameAddition0,...
        plotEnv, plotType, linewidth, scale, doPlotHistograms] =...
        getThreeCompPlotOptions(p, fileNameAddition0, plotEnvelope, scale)

% getThreeCompPlotOptions returns the default plot options for a set of 5
%   different choices to plot seismograms for secondary arrival analysis.
% Input:
%   p           : choice of 1-5:
%     p = 1: plot raw ZRT channels,
%     p = 2: plot against properly scaled y-axis
%     p = 3: plot with polarization-filtered ZRT channels,
%     p = 4: plot polarization-filtered ZRT channels plus theoretical arrivals
%     p = 5: plot polarization-filtered ZRT channels, plus arrivals and labels
%   fileNameAddition0: 
%   plotEnvelope: boolean, whether to plot the envelopes of the seismograms
%   scale       : scaling factor for all seismograms
% Output:
%   plotArrivals        : bool
%   labelArrivals       : bool
%   plotComp            : a set of {'Z','R','T'} or {'Zp','Rp','Tp'}
%   fileNameAddition0   : str
%   plotEnv             : bool
%   plotType            : one of
%   linewidth           : float
%   scale               : float
%   doPlotHistograms    : bool

    switch p
        case 1
            plotArrivals = false;
            labelArrivals = false;
            plotComp = {'Z','R','T'};
            fileNameAddition0 = [fileNameAddition0, ''];
            plotEnv = false; doPlotHistograms = false;
            plotType = 'wig'; linewidth = 0.5;
        case 2
            plotArrivals = true; labelArrivals = true;
            plotComp = {'Zp','Rp','Tp'};
            fileNameAddition0 = [fileNameAddition0, ''];
            plotEnv = false; doPlotHistograms = false;
            plotType = 'wigbyy'; linewidth = 0.5; scale = 0;
        case 3
            plotArrivals = false; labelArrivals = false;
            plotComp = {'Zp','Rp','Tp'};
            fileNameAddition0 = [fileNameAddition0, ''];
            plotEnv = plotEnvelope; doPlotHistograms = false;
            plotType = 'bwig'; linewidth = 0.1;
        case 4
            plotArrivals = true; labelArrivals = false;
            plotComp = {'Zp','Rp','Tp'};
            fileNameAddition0 = [fileNameAddition0, ''];
            plotEnv = plotEnvelope; doPlotHistograms = true;
            plotType = 'bwig'; linewidth = 0.1;
        case 5
            plotArrivals = true;
            labelArrivals = true;
            plotComp = {'Zp','Rp','Tp'};
            fileNameAddition0 = [fileNameAddition0, ''];
            plotEnv = plotEnvelope; doPlotHistograms = true;
            plotType = 'bwig'; linewidth = 0.1;
    end
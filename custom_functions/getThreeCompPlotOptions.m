function [plotArrivals, labelArrivals, plotComp, fileNameAddition0,...
        plotEnv, plotType, linewidth, scale, doPlotHistograms] =...
        getThreeCompPlotOptions(p, fileNameAddition0, plotEnvelope, scale)

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
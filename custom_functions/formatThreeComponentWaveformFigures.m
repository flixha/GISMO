function newax = formatThreeComponentWaveformFigures(ax, stationORevent,...
    plotComp, plotType, plotArrivals, plotEnvelope, labelArrivals,...
    fileNameAddition0, printFigure)

    imagePaperA4 = false;
    for k=1:1:numel(ax)
        box(ax(k),'on');
        imagePaperA4 = true;
        %shorten the title down to station + channel
        oldTitle = strsplit(ax(k).Title.String);
        ax(k).Title.String = oldTitle{1};
        
        % try to reorder children of plot such that text is always on top.
        children = ax(k).Children;
        textBoxHandles = children(isgraphics(children,'Text'));
        nonTextBoxHandles = children(~isgraphics(children,'Text'));
        ax(k).Children = [textBoxHandles; nonTextBoxHandles];
        ax(k).SortMethod = 'ChildOrder';
        drawnow;
        
        hAll = findall(ax(k));
        for idx = 1 : length(hAll)
          try
              % do not change the colorful picks lines
              if hAll(idx).Color == [0 0 0]
                set(hAll(idx),'LineWidth', 0.1);
              end
          catch
            % never mind...
          end
          try
           set(hAll(idx),'fontsize',7);
          catch
           % never mind...
          end
          
          try
              stc = hAll(idx).String;
              % change text label layout only if: 1. station gather, 2. not
              % the y-r or xlabel or title
              if ~contains(stationORevent,'event') && ~contains(...
                      stationORevent,'200') && (contains(stc,'P_') ||...
                      contains(stc,'S_')) || strcmpi(stc,'P') ||...
                      strcmpi(stc,'S')
                  set(hAll(idx),'BackgroundColor','w');
                  set(hAll(idx),'EdgeColor','k');
                  set(hAll(idx),'Margin',1); set(hAll(idx),'Margin',0);
                  uistack(hAll(idx),'top')
              end
          catch
            % never mind
          end
          
          % if it's a wiggle plot, then all wiggles should be made black
          % instead of blue
          try
             if strcmpi(plotType,'WIG') && hAll(idx).Color == [0 0 1] 
               set(hAll(idx),'Color','k');
             end
          catch
           % never mind...
          end
          
          % Remove lines from an envelope plot to reduce file size
          try 
              if strcmpi(plotType,'BWIG') && plotEnvelope && hAll(...
                      idx).Color == [0 0 0] && isgraphics(hAll(idx),'line')
                 delete(hAll(idx))
              end
          catch
              %nothing
          end
        end
    end
    for k=1:1:numel(ax)
        ax(k).Title.FontSize = 9;
        ax(k).XLabel.FontSize = 8;
        % make sure the title is in the middle
        ax(k).Title.Position(1) = sum(ax(k).XLim)/2;
        
        if k==1
            ax(k).YLabel.FontSize = 8;
        else
            ax(k).YTickLabel = {};
            ax(k).YLabel.String = '';
        end
    end
    drawnow;
    
%     if numel(ax) > numel(plotComp)
        

    np = numel(plotComp);
    newax = ax;
    if np > 1 && np >= numel(ax)
        fig = figure;
        for k=1:1:np
            newax(k) = copyobj(ax(k),fig);           
            % with y label only kept on leftmost axis:
            newax(k).Position(1) = (1-0.085)/np * (k-1)+0.08;
            newax(k).Position(3) = (1-0.085)/np - 0.01;
        end
    elseif np < numel(ax)
        % If I want to include a histogram on the right
        fig = figure;
        for k=1:1:numel(ax)
            newax(k) = copyobj(ax(k),fig);           
            % with y label only kept on leftmost axis:
            if k <= np
                newax(k).Position(1) = (1 - 0.145)/np * (k-1)+0.055;
                newax(k).Position(3) = (1-0.15)/np - 0.01;
            else
                newax(k).Position(1) = (1 - 0.145)/np * (np)+0.055;
                newax(k).Position(3) = 0.085;
                newax(k).Position(2) = newax(1).Position(2);
                newax(k).Position(4) = newax(1).Position(4);
                box(newax(k),'off');
            end
        end
    else
        fig = gcf;
    end
    
    % change font sizes
    for k=1:1:numel(ax)
        newax(k).Title.FontSize = 9;
        newax(k).XLabel.FontSize = 8;
        newax(k).YLabel.FontSize = 8;
        newax(k).YLabel.FontSize = 8;
    end

    set(fig, 'Color', 'w');
    set(fig,'Units','centimeters');
    set(fig,'PaperUnits','centimeters');
    pos = get(fig,'Position');
    
    %make a a4-page plot only for the station gathers
    if imagePaperA4
        set(fig,'Position',[pos(1), pos(2), 21.0, 29.7]);
        set(fig,'PaperPositionMode','Auto','PaperSize', [21.0, 29.7] )
    else
        for j=1:1:numel(newax)
            newax(j).Position(2) = 0.05;
            newax(j).Position(4) = newax(j).Position(4)-0.09;
            ybounds = get(newax(j),'YLim');
            newax(j).Title.Position(2) = ybounds(1) - diff(ybounds)*0.11;
        end
        set(gcf,'PaperPositionMode','Auto','PaperSize',[pos(3),pos(4)])
    end

    % chose the file name by station/event, components, processing steps
    % etc.
    compName = char(plotComp);
    compName = reshape(compName', 1, numel(compName));
    drawnow;

    if contains(compName,'p')
        fileNameAddition1 = 'PolFiltered_';
    else
        fileNameAddition1 = 'Filt_';
    end
    
    if plotArrivals && labelArrivals
        fileNameAddition2 = 'wArrivalz_Labels_';
    elseif plotArrivals
        fileNameAddition2 = 'wArrivals_';
    else
        fileNameAddition2 = '';
    end

    if printFigure
        figName = ['Figures/', fileNameAddition0, fileNameAddition1,...
            fileNameAddition2, '_', plotType, '_', stationORevent, '_',...
            compName];
        print(fig, [figName, '.pdf'], '-dpdf', '-painters', '-r600')
    
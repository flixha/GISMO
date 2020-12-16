function wiggleplot(c,scale,ord,norm,Y,reverseBlackWhite)

% Private method. See ../plot for details.

% Author: Michael West, Geophysical Institute, Univ. of Alaska Fairbanks
% $Date$
% $Revision$


% PREP PLOT
figure('Color','w','Position',[50 50 850 1100]);
box on;
% hold on;

% LOOP THROUGH WAVEFORMS
wstartrel = 86400 *( get(c.W(ord),'START_MATLAB') - c.trig(ord));% relative start time (trigger is at zero)
freq = get(c.W(ord),'Fs');
lengths = get(c.W(ord),'data_length');
tr = nan(max(lengths),numel(ord)); %pre allocate with nan to not plot
abs_max =  max(abs(c.W(ord)));
for count = 1:numel(ord)
    % time vector per trace
    tr(1:lengths(count),count) = ...
        wstartrel(count) + [ 0:lengths(count)-1]'/freq(count);   
end;


% Y=[1:numel(ord)];
if size(Y,2) > size(Y,1)
    Y = Y';
end
% scale is negative because it is reversed below
if norm==0
    % GET MEAN TRACE AMPLITUDE FOR SCALING BELOW (when norm = 0)
    maxlist = max(abs(c.W(ord)));
    normval = mean(maxlist);
    d =  double(c.W(ord) .*( -scale ./ normval)+ Y,'nan'); % do not normalize trace amplitudes
else
    abs_max(abs_max==0) = 1; % ignore zero traces
    
    d = double(c.W(ord) .* (-scale ./ abs_max)+ Y,'nan'); % normalize trace amplitudes
end

% amount of whitespace from trace mean before coloring begins
whiteAreaRatio = 0.00;
%wb = (Y' - scale) .* whiteAreRatio;
wb = scale  .* whiteAreaRatio;

%red fill (positive)
% Y: trace offset from 0; wb: no-color area; sign... is 1 if trace above threshold; d-wb is :
redmat = Y' + wb  + (sign( d - Y' - wb) + 1)/2 .* (d -Y' - wb);
%blue fill (negative)
bluemat = Y' - wb + (sign(-(d - Y' + wb)) + 1)/2 .* (d -Y' + wb);


hold on

if reverseBlackWhite
    fill(tr, redmat, 'k','LineStyle','none','EdgeColor',[1 1 1]);
else
    fill(tr, bluemat, 'k','LineStyle','none','EdgeColor',[1 1 1]);
end
% fill(tr, redmat, 'w','LineStyle','none','EdgeColor',[1 1 1]);

% plot(tr,d,'k-','LineWidth',0.2);

% fill(tr, bluemat, 'b','LineStyle','none','EdgeColor',[1 1 1]);
% fill(tr, redmat, 'r','LineStyle','none','EdgeColor',[1 1 1]);


% adjust figure
%axis([tmin tmax 0 length(ord)+1]);
axis([min(tr(:)) max(tr(:)) 0 length(ord)+1]);
set(gca,'YDir','reverse',...
    'YTick',1:length(ord),...
    'YTickLabel',datestr(c.trig(ord),'yyyy-mm-dd HH:MM:SS'),...
    'FontSize',8);

xlabel('Relative Time,(s)','FontSize',8);



% replace dates with station names if stations are different
if ~check(c,'STA')
    sta  = get(c,'STA');
    chan = get(c,'CHAN');
    
    for i=1:get(c,'TRACES')
        labels(i) = strcat( sta(i) , '_' , chan(i) );
    end
    set( gca , 'YTick' , [1:1:get(c,'TRACES')] );
    set( gca , 'YTickLabel' , labels );
end




%PRINT OUT FIGURE
set(gcf, 'paperorientation', 'portrait');
set(gcf, 'paperposition', [.25 .25 8 10.5] );
%print(gcf, '-depsc2', 'FIG_alignwfm.ps')
%!ps2pdf FIG_alignwfm.ps
%!convert FIG_alignwfm.ps FIG_alignwfm.gif
%!rm FIG_alignwfm.ps


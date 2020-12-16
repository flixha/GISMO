function TC = polarizationfilter(TC, varargin)
% POLARIZATIONFILTER calculates particle motion vectors for threecomp object
% TC = POLARIZATIONFILTER(TC) calculates particle motion parameters for
% threecomp object TC. The results are stored as threecomp properties:
%     TC.rectilinearity
%     TC.planarity
%     TC.energy
%     TC.azimuth
%     TC.inclination
%     TC.weightFunction
%
% TC = POLARIZATIONFILTER(TC, n, J, K, DT, WIDTH) applies a polarization filter
% based on a gain function defined from calculated particle motions. The
% Characteristics of the gain function are controlled by parameters n, J
% and K (see Montalbetti et al., 1970). 
% DT is the time step through the traces. WIDTH is the width of the time
% window. If the parameters are not included, the function estimates
% appropriate values for the data.
%
% POLARIZATIONFILTER can accept traces that have been rotated in the horizontal
% plane (i.e. type ZRT and Z21) so long as the orientation field is filled
% in. This is accomplished internally by first rotating these traces to
% type ZNE and then carrying out the particle motion analysis. In other
% words, the particle motion coefficients are relative to a fixed
% geographic reference frame and are independent of the orientation of the
% input traces. See THREECOMP(DESCRIBE) for a complete description of the
% particle motion fields.
%
% For an overview of this methodology see:
%    Montalbetti, J. F. and Kanasewich, E. R. (1970), Enhancement 
%    of Teleseismic Body Phases with a Polarization Filter. Geophys. 
%    J. Royal Astro. Soc., 21, 119???129. 
%    doi: 10.1111/j.1365-246X.1970.tb01771.x
%
% see also threecomp/describe

% Author: Felix Halpaap, Department of Earth Science, University of Bergen,
% after the particlemotion-script of Michael West,Geophysical Institute, 
% Univ. of Alaska Fairbanks
% $Date$
% $Revision$



if isa(TC,'threecomp')
    originalTraces = TC.traces;
else
    error('Threecomp:polarizationfilter:mustBeThreecompObject',...
        'First argument must be a threecomp object');
end
    
    
% SET UP INPUTS
%n
if length(varargin) >= 1
    n = varargin{1};
else
    n = 0.5;
end

% J
if length(varargin) >= 2
    J = varargin{2};
else
    J = 1;
end

% K
if length(varargin) >= 3
    K = varargin{3};
else
    K = 2;
end

% l
if length(varargin) >= 4
    dt = varargin{4};
else
    dt = 86400 * get(TC(1).traces(1),'DURATION' ) / 200; %100
end

% w
if length(varargin) >= 5
    width = varargin{5};
else
    width = 86400 * get(TC(1).traces(1),'DURATION' ) / 10; % 10
    disp(['n: ' num2str(n) ' J: ' num2str(J) ' K: ' num2str(K),...
        'Time step: ' num2str(dt,'%4.3f') ' Window width: ',...
        num2str(width,'%4.3f') ]);
end




% CHECK ORIENTATIONS
orientation = get(TC,'ORIENTATION');
if isempty(orientation) || any(any(isnan(orientation)))
    error('Threecomp:polarizationfilter:requiresOrientation',...
        'Orientations must be provided for all channels');
end
if any(orientation(:,2)~=0)
    error('Threecomp:polarizationfilter:verticalCompMustPointUp',...
        'vertical component must have a vertical orientation of 0.');
end



% TODO: ROTATE TEMPORARY THREECOMP TO ZNE, IF NEEDED


% STEP THROUGH EACH OBJECT
%disp('applying polarization filter ->      ');
qMax = numel(TC);
for q = 1:qMax
    %fprintf('\b\b\b\b\b%4.0f%%',q/qMax*100);
    [Z, R, T] = do_one(TC(q), n, J, K, dt, width);
    TC = set(TC,'waveform',[Z, R, T]);
end
%fprintf('\n');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process a single object

function [Z, R, T, rec, plan, energy, az, incl] = do_one(TC, n, J, K, dt, width)

% TEMPORARILY ADJUST ORIENTATION IF NEEDED
chn1 = get(TC.traces(1),'ChannelTag');
chn2 = get(TC.traces(2),'ChannelTag');
chn3 = get(TC.traces(3),'ChannelTag');
if (strcmp(chn2.channel(end),'N') && strcmp(chn3.channel(end),'E'))
    if isempty(TC.backAzimuth)
        error(['Orientations are ZNE, need to rotate to ZRT but no ',...
            'backazimuth provided'])
    else
        disp(['temporarily rotating traces to Z-R-T for polarization ',...
            'filtering ...']);
        TC = rotate(TC);
    end
end
w = TC.traces;
freq = get(w(1), 'Freq');
alldata = double(w);
%ellipsoidMatrix = cov(alldata);
Tstart = get(w(1), 'Start_Matlab');
Tend = get(w(1), 'End_Matlab');
Tc = Tstart : dt/86400 : Tend;
T1 = Tc - 0.5*width/86400;
T2 = Tc + 0.5*width/86400;
f = (T1 < Tstart);
T1(f) = Tstart;
f = (T2 > Tend);
T2(f) = Tend;
T1i = round((T1-Tstart)*86400*freq + 1);
T2i = round((T2-Tstart)*86400*freq);
% test = [(T2'-T1')*86400 T2i' T1i'];
% Tc is the center of the time ranges bounded by T1 and T2
% T1i and T2i are the sample indexes corresponding to these time ranges
% Is this formulation correct to the sample?


steps = numel(T1);
pm.Rec = zeros(1,steps); RL = pm.Rec; Dir = zeros(3,steps);
% pm.RL = pm.Rec;
% pm.Dir = zeros(3,steps);

% fprintf('Calculating covariance matrix : %05d of %05d',0, steps);
% backspacing = repmat('\b',1,15);
% formstr = [backspacing,' %05d of %05d'];
 
for q=1:steps
    %fprintf(formstr,n, totcount);
    snippet = alldata( T1i(q) : T2i(q) ,:); %grab a chunk
    %snippet = snippet .* repmat(hanning(size(snippet,1)),1,3);
    covmatrix = cov(snippet); %covariance matrix
    [V, D] = eig(covmatrix); %columns of V are eigenvectors, D are eigenvalues
    lambda = diag(D); %lambda are eigenvlues
    [lambda_sorted, I] = sort(lambda,1,'descend'); %lambda are sorted eigenvlues, I is index
    EprincAx = V(:,I(1)); % eigenvector of the principal axis with respect to Z,R,T
    if numel(lambda) < 3
        warning('less than three eigenvalues!');
    end   
   
    % RECTILINEARITY
    pm.Rec(q) = 1 - lambda_sorted(2) / lambda_sorted(1);   
    %pm.Rec(n) = 1 - (lambda(2)+lambda(3))/(2*lambda(1)); % alternate definition
    
    % Rectilinearity gain function for Z,R,T-components
    RL(q) =  (1 - (lambda_sorted(2) / lambda_sorted(1)).^n).^J;
    % Directivity gain function for Z,R,T-components
%     Dir(:,n) = lambda.^K; 
    Dir(:,q) = EprincAx.^K;
end

% smooth the gain functions within windows that are half the width of the
% original window
halfwidth = width/2;
nhalfwidth = halfwidth/dt;
% If nhalfwidth cannot be divided by two, then just make the halfwidth by 1
% sample wider
if mod(nhalfwidth,2)>0
    nhalfwidth = nhalfwidth +1;
end

steps = numel(RL);
lowerstep = (1:1:steps) - nhalfwidth/2;
upperstep = (1:1:steps) + nhalfwidth/2;
f = (lowerstep < 1);
lowerstep(f) = 1;
f = (upperstep > steps);
upperstep(f) = steps;
% smooth each sample with running mean
for q=1:steps
    RLsnippet = RL(1, lowerstep(q) : upperstep(q)); %grab a chunk
    RL(q) = mean(RLsnippet);
    Dsnippet = Dir(:, lowerstep(q) : upperstep(q)); %grab a chunk
    Dir(:,q) = mean(Dsnippet,2);
end

% if dt is not the same as sampling interval, then interpolate the gain
% functions
if dt ~= 1/freq
    newdt = (1/freq)/86400;
    newTc = Tc(1):newdt:Tc(end);
    RL = interp1(Tc, RL, newTc, 'spline');
    newDir = zeros(size(RL));
    for p=1:1:3
        newDir(p,:) = interp1(Tc, Dir(p,:), newTc, 'spline');
    end
    Dir = newDir;
end


% shorten traces by the last sample to work around fence post problem
RL = RL(1:end-1);
Dir = Dir(:,1:end-1);

% Control gain on Z, R and T-components
% R = R.* pm.F.* pm.Dir;
gaineddata = alldata .* RL'.* Dir';

Z = w(1);
Z = set(Z,'data',gaineddata(:,1),'units','');
Z = addhistory(set(Z,'history',{}),'Gain adjusted by polarization filter');

R = w(2);
R = set(R,'data',gaineddata(:,2),'units','');
R = addhistory(set(R,'history',{}),'Gain adjusted by polarization filter');

T = w(3);
T = set(T,'data',gaineddata(:,3),'units','');
T = addhistory(set(T,'history',{}),'Gain adjusted by polarization filter');









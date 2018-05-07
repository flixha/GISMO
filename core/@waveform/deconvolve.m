function w = deconvolve(w, d, varargin)
    %DECONVOLVE a waveform ("wavelet") from the waveform object(s).
   
    %Deconvolving a waveform of length n with a wavelet of length m will
    %return a waveform of length n+m-1
    
    % AUTHOR of deconvolution functions: Kathrin Spieker, Unversity of 
    %        Bergen
    % IMPLEMENTED into GISMO by Felix Halpaap, University of Bergen
    %
    % Spieker et al. (2018, XXX) provide a review of all implemented 
    % deconvolution methods. They discuss for which applications each method
    % is particularly suitable.
    % $Date$ March 2018
    % $Revision$
    
    p = inputParser;
    % Available Methods 
    defaultMethod = 'spectraldivision';
    validMethods = {'iterative','spectraldivision','timedomain'};
    checkMethod = @(x) any(validatestring(x,validMethods));  
    addOptional(p,'method',defaultMethod,checkMethod)
    
    % Parameter-value pairs
    defaultIterations = 100;
    addParameter(p,'iterations',defaultIterations,@isnumeric)
    defaultGaussianWidth = 2.5;
    addParameter(p,'gaussianWidth',defaultGaussianWidth,@isnumeric)
    defaultMinimalErrorChange = 0.001;
    addParameter(p,'minimalErrorChange',defaultMinimalErrorChange,@isnumeric)
    
    defaultRegularization = 'wat';
    validRegularization = {'con','wat','fqd'};
    checkRegularization = @(x) any(validatestring(x,validRegularization));
    addParameter(p,'regularization',defaultRegularization,checkRegularization)
    
    defaultTShift = 0;
    addParameter(p,'tshift',defaultTShift,@isnumeric)
    defaultBackTransformShift = 0;
    addParameter(p,'backtransformshift',defaultBackTransformShift,@isnumeric)
    
    % Parse function input
    parse(p, varargin{:})
    % Display what options were chosen
    if ~isempty(fieldnames(p.Unmatched))
       disp('Extra inputs:')
       disp(p.Unmatched)
    end
    if ~isempty(p.UsingDefaults)
       disp('Using defaults: ')
       disp(p.UsingDefaults)
    end
    
       
    %Check the number of waveforms vs. wavelets
    Nwavelets = numel(d);
    Nmax = numel(w);
    if Nwavelets ~= 1
        if Nwavelets ~= Nmax
            msg = ['Need only one deconvolution wavelets or the same ',...
                'number as waveforms'];
            error(msg)
        end
    end
    
    
    for j = 1:1:Nmax
        if isempty(w(j)) || sum(w(j).data)==0, continue, end
        %chose appropriate wavelet
        if Nwavelets > 1
            if isempty(d(j)), continue, end
            wavelet = d(j).data;
        else
            wavelet = d.data;
        end
        
        % zero-pad inputs if the length of waveform and wavelet differ
        n_samples_w = length(w(j).data);
        if n_samples_w > length(wavelet)
            wavelet(end:n_samples_w) = zeros;
        end

        % sampling rate
        ndt = 1/w(j).Fs;
        
        switch lower(p.Results.method)
          case 'iterative'
            [rf,rms] = iterdecon(w(j).data', wavelet', p.Results.iterations,...
                p.Results.minimalErrorChange, p.Results.backtransformshift,...
                ndt, p.Results.gaussianWidth);
            w(j).data = rf';
          case 'spectraldivision'
            [rf,qrf] = spectraldivision(w(j).data, wavelet, ndt, ...
                p.Results.tshift, p.Results.regularization);
          case 'timedomain'
            [rf,qrf,eps] = timedeco_noise(w(j).data', wavelet', ndt,...
                p.Results.tshift);
            w(j).data = rf;
        end
        
        
    end
   
end


function [vrf,rms]=iterdecon(v,u,iter,mini,shift,ndt,gauss)
    %Function LIGORRIA is the iterative time domain deconvolution described by
    %Ligorria & Ammon (1999)
    %Input: 
    %       v:      trace which should be deconvolved (numerator)
    %       u:      direct wave for deconvolution (denominator)
    %       iter:   maximum iterations
    %       t2:     end time (s)
    %       mini:   minimal error change stops iteration
    %       shift:  time shift for back transformation
    %       width:  filter width (generally 2.5 which is equal to a pulse width 
    %               of 1s)
    %Output:
    %       rf:     receiver function
    %
    % The following functions are needed: correlat.m, gfilter.m, convol.m,
    % sshift.m, gaussian.m

    %number of points in fft
    N = length(v);
    nf = 2^nextpow2(N);
    %filter
    %[gauss]=gaussian(ndt,nf,width);

    %root-mean-square vector
    rms = zeros(iter,1); 
    %vector spikes
    s = zeros(1,nf);
    u(N+1:nf) = zeros;
    v(N+1:nf) = zeros;
    %residual numerator
    res = v;
    %variables for iteration
    it = 0; 
    er = mini+100*sum(v.^2); %error change
    rms_old = 1;
    maxi = nf; 
    %begin iteration loop
    while (abs(er) > mini && it < iter)
        it = it+1; 
        %cross correlation of numerator u (here res) and denominator w
        C = correlat(res,u,nf);
        C = C/sum(u.^2);
        %amplitude
        [~,dum] = max(abs(C(1:maxi)));
        %ampli=C(dum)/ndt;
        ampli = C(dum);
        %deconvolution
        s(dum) = s(dum)+ampli; %create spike
        s_new = gfilter(s,gauss,nf,ndt);
        %convolve spike and denominator wavelet w to get u resp. res
        res_new = convol(s_new,u,nf,ndt); 
        res = v-res_new;
        %adjust error
        rms_new = sum(res.^2);
        rms(it) = rms_new; %root-mean-square vector
        er = (rms_old-rms_new)*100; %change in error
        rms_old = rms_new;
    end
    [num2str(abs(er)), ' error after ', num2str(it), ' iterations']
    rms = sqrt(rms(it)/numel(res));
    res_new = gfilter(s,gauss,nf,ndt);
    %shift
    res_new = sshift(res_new,nf,ndt,shift);
    %receiver function
    vrf = res_new(1:N);
end



function [lrf,qrf]=spectraldivision(w,u,ndt,tshift,regul)
    % Function spectraldivision(w,u,ndt,tshift,regul) is a standard spectral
    % division frequency domain deconvolution. The regularization can be chosen
    % by the variable "regul", this can be 'con', 'wat', or 'fqd' for constant
    % damping factor, waterlevel, or frequency-dependent damping, respectively.
    % "ndt" is the sampling rate in seconds and "tshift" is the time before
    % onset of the direct wave.
    % w=lcomp;
    % u=qcomp;
    %
    % output is the receiver function "qrf" and the direct wave deconvolved by
    % itself
    N = length(w);
    % pre-event noise needed for regularisation parameter
    wn = zeros(1,N);
    wn(1:(tshift-3)/ndt) = w(1:(tshift-3)/ndt);
    %number of points in fft
    nf = 2^nextpow2(N);
    nft = nf/2+1;
    %fourier transform
    uf = fft(u,nf);
    wf = fft(w,nf);
    wnf = fft(wn,nf);
    %denominator and regularisation
    den = wf.*conj(wf);
    noise = wnf.*conj(wnf);
    den0 = den;
    % which regularization do you want?
    freqdep = strcmp(lower(regul),'fqd');
    const = strcmp(lower(regul),'con');
    water = strcmp(lower(regul),'wat');
    if (freqdep==0) && (const==0) && (water==0)
       error(['Regularization not defined (your input: regul=',regul ...
           ').' ...
           ' Use either "fqd" for frequency-dependent' ...
           'or "con" for constant value regularization' ...
           'or "wat" for water-level.'])
    end
    % constant damping factor regularization
    if const
       eps = max(real(noise));
       den = den+eps; 
    end
    % frequency-dependent damping 
    if freqdep
       den = den+noise;
    end
    % waterlevel regularization
    if water
        eps = (max(real(noise)));
        den(real(den)<eps) = eps;
    end
    % numerator
    num = uf.*conj(wf);
    numl = wf.*conj(wf);
    % deconvolution
    rfq = num./den;
    rfl = numl./den;
    if freqdep==1
         N2=floor(numel(rfq)/2)+1;
        for i=1:N2
            fac=cos(pi/2*(i-1)/N2)^2;
            rfq(:,i)=rfq(:,i)*fac;
        end
    end
    % back transformation
    w=(0:nf/2)*2*pi/(nf*ndt);
    w=[w,-fliplr(w(2:end-1))];
    rfq=rfq.*exp(-sqrt(-1)*w*tshift);
    rfl=rfl.*exp(-sqrt(-1)*w*tshift);

    qrft=ifft(rfq,nf);
    qrf=real(qrft);
    qrf=qrf(1:N);
    lrft=ifft(rfl,nf);
    lrf=real(lrft(1:N));
    lrf=lrf(1:N);
end


function [lrf,qrf,eps]=timedeco_noise(u,w,ndt,tshift)
    % Implementation of the time deconvolution described by Gurrola 1995
    % "Simultaneous time-domain deconvolution with application  to  the
    % computation of receiver functions".
    % u: numerator, w: denominator, eps: weighting/regularisaton parameter, 
    % ndt: sampling interval, tshift: time before direct wave onset
    % The following function is needed: glev.m
    N = length(w);
    M = length(u);
    N2 = 2*N-1;
    %number of points in fft
    nf = 2^nextpow2(N);
    w = w';
    %noise; regularisation parameter eps is set to the maximum of the pre-event noise
    %% WHY substract 2 from the time in seconds? noise = w(1:(tshift-2)/ndt);
    noise = w(1:(tshift)/ndt-2);
    noiseM = convmtx(noise,length(noise));
    noised = noiseM'*noiseM;
    eps = max(max(noised));

    u01 = zeros(tshift/ndt,1);
    u02 = zeros(N-tshift/ndt-1,1);
    u2 = [u01; u'; u02];
    w2 = [u01; w ; u02];
    %creating convolution matrix
    CM = convmtx(w,N); 
    dumm1 = (eye(N)+1/eps*(CM'*CM));

    %Calculate rf
    dummq = 1/eps*(CM'*u2);
    qrff = glev(dumm1,dummq);
    qrf = qrff(1:N);
    dumml = 1/eps*(CM'*w2);
    lrff = glev(dumm1,dumml);
    lrf = lrff(1:N);
end


function [ab]=convol(a,b,nf,ndt)
    %convolution
    afr=fft(a,nf);
    bfr=fft(b,nf);

    afr=afr.*bfr*ndt;

    ab=real(ifft(afr,nf));
    return
end

function [b]=sshift(a,nf,ndt,shift)
    %convert to frequency domain
    afr=fft(a,nf);

    s=round(shift/ndt);
    p=2*pi*(1:nf).*s/nf;
    afr=afr.*(cos(p)-1i.*sin(p));
    %back to time-domain
    b=real(ifft(afr,nf))/cos(2*pi*s/nf);
    return
end


function x=glev(r,b)
    %GLEV   Levinson recursion.
    %----
    %USAGE: x = glev(r,b)
    %
    %       The Levinson recursion solves the Toeplitz equations
    %               R x = b
    %       where R=toeplitz(r) and b is an arbitrary vector.
    %       
    %  see also ATOG, ATOR, GTOA, GTOR, RTOA, RTOG
    %
    %---------------------------------------------------------------
    % copyright 1996, by M.H. Hayes.  For use with the book 
    % "Statistical Digital Signal Processing and Modeling"
    % (John Wiley & Sons, 1996).
    %---------------------------------------------------------------

    r=r(:);
    p=length(b);
    a=1;
    x=b(1)/r(1);
    epsilon=r(1);
    for j=2:p;
        g=r(2:j)'*flipud(a);
        gamma(j-1)=-g/epsilon;
        a=[a;0] + gamma(j-1)*[0;conj(flipud(a))];
        epsilon=epsilon*(1 - abs(gamma(j-1))^2);
        delta=r(2:j)'*flipud(x);
        q=(b(j)-delta)/epsilon;
        x=[x;0] + q*[conj(flipud(a))];
    end
end


function [x]=correlat(u,w,nf)
    x=real(ifft(fft(u,nf).*conj(fft(w,nf)),nf));
    return
end

function [gauss]=gaussian(ndt,nf,width)
    %function to create a Gaussian filter with sampling intervall (ndt),
    %number of points in frequency domain (nf) and width of the filter (width).
    %Output is a Gaussian filter of form G=exp(-w^2/(4*width^2))
    df=1/(ndt*nf);
    nft=0.5*nf+1;
    f=df*(0:1:nft-1);
    w=2*pi*f;

    gauss=zeros(1,nf);
    gauss(1:nf/2+1)=exp(-w.^2/(4.*width^2))/ndt;
    gauss(nf/2+2:end)=fliplr(gauss(2:nf/2));
    gauss(nf/2+1)=0.;
end

function [ab]=gfilter(a,b,nf,ndt)
    afr=fft(a,nf);
    afr=afr.*b*ndt;
    ab=real(ifft(afr,nf));
    return
end


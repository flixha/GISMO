function [vrf,rms]=iterdecon(u,v,iter,mini,shift,ndt,gauss)

    %Function LIGORRIA is the iterative time domain deconvolution described by
    %Ligorria & Ammon (1999)
    %Input: 
    %       u:      direct wave for deconvolution (denominator)
    %       v:      trace which should be deconvolved (nominator)
    %       iter:   maximum iterations
    %       t2:     end time (s)
    %       mini:   minimal error change stopps iteration
    %       shift:  time shift for back transformation
    %       width:  filter width (generally 2.5 which is equal to a pulse width 
    %               of 1s)
    %Output:
    %       rf:     receiver function
    %
    % The following functions are needed: correlat.m, gfilter.m, convol.m,
    % sshift.m, gaussian.m

    %number of points in fft
    N=length(u);
    nf=2^nextpow2(N);

    %filter
    %[gauss]=gaussian(ndt,nf,width);

    %root-mean-square vector
    rms=zeros(iter,1); 

    %vector spikes
    s=zeros(1,nf);

    %u=gfilter(u,gauss,nf,ndt);
    %v=gfilter(v,gauss,nf,ndt);
    %
    u(N+1:nf)=zeros;
    v(N+1:nf)=zeros;

    %residual numerator
    res=v;

    %variables for iteration
    it=0; 
    er=mini+100*sum(v.^2); %error change
    rms_old=1;
    maxi=nf; 

    %begin iteration loop
    while (abs(er) > mini && it < iter)
        it=it+1; 

        %cross correlation of nominator u (here res) and denominator w
        C=correlat(res,u,nf);
        C=C/sum(u.^2);

        %amplitude
        [~,dum]=max(abs(C(1:maxi)));
        %ampli=C(dum)/ndt;
        ampli=C(dum);

        %deconvolution
        s(dum)=s(dum)+ampli; %create spike
        s_new=gfilter(s,gauss,nf,ndt);
        %convolve spike and denominator wavelet w to get u resp. res
        res_new=convol(s_new,u,nf,ndt); 
        res=v-res_new;

        %adjust error
        rms_new=sum(res.^2);
        rms(it)=rms_new; %root-mean-square vector
        er=(rms_old-rms_new)*100; %change in error
        rms_old=rms_new;

    end

    rms=sqrt(rms(it)/numel(res));

    res_new=gfilter(s,gauss,nf,ndt);

    %shift
    res_new=sshift(res_new,nf,ndt,shift);

    %receiver function
    vrf=res_new(1:N);
end
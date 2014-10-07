function [Reg Reintcf Retcf Img Imintcf Imtcf]=bath_construction_TCY(J, F, KbT, t_max, dtg, dtt) 
%This function constructs the bath for simulation. The output includes the
%real parts and imaginary parts of lineshape function (g), the first integral 
%of the time correlation function (intcf), and time correlation function (tcf)

%Retcf is calcualted by FFT of input function, Redtcf.
%Reintcf is calcualted by numerically integrates Retcf.
%Img and Imintcf are calcualted by FFT.
%Both Reg and Img are calculated by propagating the time correlation function
%and its integral.

%INPUT
%Redtcf: for real part tcf, which is spectral density function J(w) * coth(w/(2KbT))
%Imdtcf: for imaginary part tcf, which is negative spectral density -J(w)
%Imdintcf: for imaginary part tcf integral, which is actually J(w)/w
%t_max: the maximum time calculated for the bath, depending on the coupling
%strength. For example, if couplingC <= 5, t_max must be larger than 15000
 
%Time Axis Assignment
%lineshape function
% dtg = 0.1; 
tg = 0:dtg:t_max;

%tcf and intcf
% dtt=0.1;
tt = 0:dtt:t_max;

%REAL PART
%==========================================================================

%Time Correlation Function Re. part
%==================================
%use ifft
%set the frequency domain vector, unit in fs-1
%Note that the 1st point in frequency domain must be avoided due to its
%singularity from coth().
dww=1/10/t_max*dtt;
%ww=dww*2*pi:(2*pi*dww):2*pi/dtt;
ww=(dww:dww:1/dtt)*2*pi; 

Redtcf=@(w) J(w).*coth(w/(2*KbT));
tcf_w0=2*F(0)*KbT; %prepare for the point at 0;

tcf_fft = 2*pi*(1/dtt)*real(ifft([tcf_w0 Redtcf(ww)]));
Retcf=tcf_fft(1:t_max/dtt+1);                         %select the time range

% Retcf=Retcf-Retcf(end);

% dtcf = (Retcf(2:end)-Retcf(1:end-1))/dtt;
% dtcf=[dtcf, 0];
% Retcf(find(abs(dtcf)<=1e-8))=0;

%TCF integral Re. part
%=====================
%integration of tcf with spline polynomial
ptcf=spline(tt, Retcf);
pintcf = ptcf;                      %make structure for intcf
coe = pintcf.coefs;                 %save the coefficient for integration
ndt   = pintcf.pieces;              %ndt is the pieces;
pintcf.order = 5;                   %assign the order after integration

%integration
pintcf.coefs = [coe(:,1)/4 coe(:,2)/3 coe(:,3)/2 coe(:,4) zeros(ndt, 1)];
coe = pintcf.coefs;
coefs_sum = [coe(:,1)*dtt^4 coe(:,2)*dtt^3 coe(:,3)*dtt^2 coe(:,4)*dtt];
pintcf.coefs(2:ndt,5)= cumsum(sum(coefs_sum(1:ndt-1,:),2),1);

Reintcf=ppval(tt, pintcf);

%Lineshape Function Re. part
%===========================
%propagation
gstep =ppval(tg, pintcf)*dtg+0.5*ppval(tg, ptcf)*dtg^2;
Reg = cumsum(gstep, 2); 

%IMAGINARY PART
%==========================================================================

%Time Correlation Function Im. part
%==================================
%use ifft
%reset the frequency domain axis, since there's no singularity to worry
%about anymore XD
Imdtcf =@(w) -J(w);

ww =(0:1/2/t_max:1/dtt)*2*pi;    

Imddg_ifft = 2*pi/dtt*ifft(Imdtcf(ww));
Imddg_tt = imag(Imddg_ifft);

Imtcf=Imddg_tt(1:t_max/dtt+1);
    
%TCF integral Im. part
%=====================
Imdintcf =@(w) F(w);

Imdg_ifft = 2*pi/dtt*ifft(Imdintcf(ww));
Imdg_tt = real(Imdg_ifft);

Imintcf=Imdg_tt(1:t_max/dtt+1)-real(Imdg_ifft(1))*ones(1,t_max/dtt+1);
 
%Lineshape Function Im. part
%======================================================================
gstep =spline(tt, Imintcf, tg)*dtg+0.5*spline(tt, Imtcf, tg)*dtg^2;
Img = cumsum(gstep, 2); 

end
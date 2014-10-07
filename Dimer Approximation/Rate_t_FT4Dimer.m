function [ pR ] = Rate_t_FT4Dimer(E, C, lambda, Reg, Img, dtg, dtt)
%INPUT:
%E: exciton energy vector, unit in fs-1
%Jcoe: coefficient matrix, dimension one represents the sites.
%Reg (lineshape real part, dt=0.1): each row contains the bath for each site.
%Redg(lineshape first derivative real part)
%Reddg(time correlation function realpart)
%Img (lineshape imaginary part, dt=0.1): each row contains the bath for each site.
%Imdg(lineshape first derivative imaginary part)
%Imddg(time correlation function imaginary part)

%OUTPUT:
%pR: a 3D structure of rate matrix predicted with modified Redfield theory, unit in fs-1
%   dimension 1 and 2 are elements indicating transfer processes
%   between exciton states; and dimension 3 are time, along which the
%   rate are changed.
%Rpd: the Markovian pure-dephasing rate matrix, dimension = 2.
%cRpd:  the Non-Markovian pure-dephasing rate coefficient matrix, dimension = 2.
%   The usage of cc is for time-dependent pure-dephasing rate Rpd(t).
%   Each element of Rpd(t)=cRpd*Redg(t).
%   All of the above matrice records EET rate from column sides to row
%   sides. e.g. R(1,2) means from exciton 2 to exciton 1.

%Rate integration:
%unit of R sould be fs-1
%======================================================================
tmax = (size(Reg,2)-1)*dtg; %**maxima time calculated in bath
%     dt = tmax/(size(Redg,2)-1);

tt = 0:dtt:tmax; %**time steps saved for each Redg, Imdg

nt=size(tt,2);

dR = zeros(1,nt);
g = Reg(1:dtt/dtg:end)+1i*Img(1:dtt/dtg:end);

dR = 2*real( exp(1i*(E(2)-E(1)-(C(1)+C(2))*lambda)*tt-(C(1)+C(2))*g) );

%Analytically integrate the rate equation==============================
pR = spline(tt, dR);             %make a structure for pR
clear dR

pR.order = 5;                  %assign the order after integration
coe = pR.coefs;                %save the coefficient for integration
ndt = pR.pieces;               %ndt is the pieces;

%integration
pR.coefs = [coe(:,1)/4 coe(:,2)/3 coe(:,3)/2 coe(:,4) zeros(size(coe,1), 1)];
coe = pR.coefs;
coefs_sum = [coe(:,1)*dtt^4 coe(:,2)*dtt^3 coe(:,3)*dtt^2 coe(:,4)*dtt];
pR.coefs = reshape(pR.coefs, 1, ndt, 5);
coefs_sum = reshape(coefs_sum, 1, ndt, 4);
pR.coefs(:, 2:ndt,5) = cumsum(sum(coefs_sum(:, 1:ndt-1, :), 3),2);
pR.coefs = reshape(pR.coefs, ndt, 5);

end
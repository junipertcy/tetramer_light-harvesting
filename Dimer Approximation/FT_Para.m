function [ pR ] = FT_Para( E, theta, t_max )
%0. Parameter Assignment
%==========================================================================
%Hamiltonian: either constructed here or load a Hamiltonian mat file

ne = 2;

%Unit transfer factor of cm-1 to fs-1
unitf=2*pi*3e-5;

%Bath Parameter
% unit in cm-1
%          1   2    3    4
%         wc   Lo   de  dis
bath = [100 35 0 0];
% wc: Ohmic cut-off frequency; for J_exp, it means lorentzian coupling
% Lo: Ohmic reorganization energy, coupling constant gg= Lo/wc; for J_exp
%     it means reorganization energy
% de: FWHM of site energy distribution
% dis: disorder number, 0 means no disorder

Kb = 0.6950*unitf;
T = 300;

%t_max=2e4; % time length of bath
dtg = 0.1;
dtt = 0.1;

%1.Bath Preparation
%==========================================================================
%Choose an exp or an Ohmic bath:
% Ohmic Bath parameters

wc = bath(1)*unitf;
Lo = bath(2)*unitf;

%spectral density
J=@(w) Lo/wc*w.*exp(-w/wc);
F=@(w) Lo/wc*exp(-w/wc);

%Bath Calculation
[ Reg, ~, ~, Img, ~, ~ ] = bath_construction(J, F, Kb*T, t_max, dtg, dtt);

%Excitons for unit in fs-1
W = E*unitf;
C = 1-0.5*(sin(theta)).^2;


%%
%3. Rate Calculation (preparing a rate matrix)
%==========================================================================
%Modified Redfield Rate
%Energy Transfer Rate and Pure-Dephasing Rate, output in fs-1
%tic;
[ pR ] = Rate_t_FT4Dimer(W, C, Lo, Reg, Img, dtg, dtt);
%toc;
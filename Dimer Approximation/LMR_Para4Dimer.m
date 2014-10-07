function [ pR ] = LMR_Para4Dimer( H, t_max )
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
[Reg Redg Reddg Img Imdg Imddg]=bath_construction(J, F, Kb*T, ...
    t_max, dtg, dtt);

L = -Imdg(:,end);

%%
%2. SYSTEM PREPARATION
%==========================================================================
%Exciton Basis, cm-1
[coe_mat Ee]=eig(H);
Ee=diag(Ee);                       %column vector of excitons
Jcoe = zeros(ne, ne, ne, ne, ne);  %overlap functions
for i = 1 : ne
    for j = 1 : ne
        for k = 1 : ne
            for l = 1 : ne
                Jcoe(:,i,j,k,l)=coe_mat(:,i).*coe_mat(:,j).*coe_mat(:,k).*coe_mat(:,l);
            end
        end
    end
end

%Excitons for unit in fs-1
W=Ee*unitf;

%%
%3. Rate Calculation (preparing a rate matrix)
%==========================================================================
%Modified Redfield Rate
%Energy Transfer Rate, output in fs-1
[pR c_PD R_PD] = Rate_t_LMR(W, Jcoe, Reg, Redg, Reddg, Img, Imdg, Imddg, dtg, dtt);
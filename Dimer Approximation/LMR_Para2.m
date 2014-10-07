function [ H , pR , c_PD , pRedg ] = LMR_Para2( r , theta )
%0. Parameter Assignment
%==========================================================================
%Hamiltonian: either constructed here or load a Hamiltonian mat file

% H is the final input, unit must be cm-1
%load('Tetramer_Hamiltonian', 'H02');
%H=H02;
%%
%Assign system size
%ne = size(H,1);
tic;
% in this section, Hamiltonian is defined in units of cm^-1
ne = 4; % dimension of H
H = zeros(ne);

H(1,1) = 1000;
H(2,2) =  950;
H(3,3) =  860;
H(4,4) =  830;

R = 40*10^(-10); % length between D1 and A2, in units of m 
mu = 7.75*3.34*10^(-30); % electric dipole, in units of C*m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caution: In DipoleCoupling(mu1,mu2,theta1,theta2,r), theta1,% 
%          theta2 are angles between vec(mu)_j and vec(r),    %
%          while theta here is angle between two donors'      %
%          electric dipole moments and horizontal line.       %                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nearest-neighbouring-sites coupling
H(1,2) = DipoleCoupling( mu , mu , pi/2-theta , pi/2-theta , r );
H(2,1) = H(1,2);
H(2,3) = DipoleCoupling( mu , mu , pi/2-theta , pi/2 , R-2*r );
H(3,2) = H(2,3);
H(3,4) = DipoleCoupling( mu , mu , pi/2 , pi/2 , r );
H(4,3) = H(3,4);

% next-nearest-neighbouring-sites coupling
H(1,3) = DipoleCoupling( mu , mu , pi/2-theta , pi/2 , R-r );
H(3,1) = H(1,3);
H(2,4) = H(1,3);
H(4,2) = H(2,4);

% end-sites coupling
H(1,4) = DipoleCoupling( mu , mu , pi/2-theta , pi/2 , R );
H(4,1) = H(1,4);

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

Kb=0.6950*unitf;
T = 300;


t_max=2e4; % time length of bath
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

% % %function to be integrated
% % Rtcf=@(w) J(w).*coth(w/(2*Kb*T));
% % Imdtcf =@(w) -J(w);
% % Imdintcf =@(w) Lo/wc*exp(-w/wc) ;

%Plot spectral density
% ww=0:2000;
% figure; plot(ww, J(ww*unitf)/unitf);

%Bath Calculation
[Reg Redg Reddg Img Imdg Imddg]=bath_construction(J, F, Kb*T, ...
    t_max, dtg, dtt);

L = -Imdg(:,end);

J12 = H(1,2);
e2p = (H(1,1)+H(2,2))/2-sqrt(((H(1,1)-H(2,2))/2)^2+J12^2);
theta12 = atan(2*J12/(H(1,1)-H(2,2)));

J34 = H(3,4);
e3p = (H(3,3)+H(4,4))/2+sqrt(((H(3,3)-H(4,4))/2)^2+J34^2);
theta34 = atan(2*J34/(H(3,3)-H(4,4)));

J31 = H(3,1);
J32 = H(3,2);
J41 = H(4,1);
J42 = H(4,2);
J32p = (J31*sin(theta12/2)-J32*cos(theta12/2))*cos(theta34/2)...
    +(J41*sin(theta12/2)-J42*cos(theta12/2))*sin(theta34/2);

E1p = (e2p+e3p)/2+sqrt(((e2p-e3p)/2)^2+J32p^2);
E2p = (e2p+e3p)/2-sqrt(((e2p-e3p)/2)^2+J32p^2);
theta = atan(2*J32p/(e2p-e3p));

np = 2;
a = zeros(ne,np,np);
a(1,1,1) = (cos(theta))^2*(sin(theta12/2))^2;
a(2,1,1) = (cos(theta))^2*(cos(theta12/2))^2;
a(3,1,1) = (sin(theta))^2*(cos(theta34/2))^2;
a(4,1,1) = (sin(theta))^2*(sin(theta34/2))^2;

a(1,2,2) = (sin(theta))^2*(sin(theta12/2))^2;
a(2,2,2) = (sin(theta))^2*(cos(theta12/2))^2;
a(3,2,2) = (cos(theta))^2*(cos(theta34/2))^2;
a(4,2,2) = (cos(theta))^2*(sin(theta34/2))^2;

a(1,1,2) = (1/2)*sin(2*theta)*(sin(theta12/2))^2;
a(2,1,2) = (1/2)*sin(2*theta)*(cos(theta12/2))^2;
a(3,1,2) = -(1/2)*sin(2*theta)*(cos(theta34/2))^2;
a(4,1,2) = -(1/2)*sin(2*theta)*(sin(theta34/2))^2;

a(1,2,1) = a(1,1,2);
a(2,2,1) = a(2,1,2);
a(3,2,1) = a(3,1,2);
a(4,2,1) = a(4,1,2);

Jcoe = zeros(ne, np, np, np, np);  %overlap functions
for i = 1 : np
    for j = 1 : np
        for k = 1 : np
            for l = 1 : np
                Jcoe(:,i,j,k,l)=a(:,i,j).*a(:,k,l);
            end
        end
    end
end

W = [ E1p E2p ]*unitf;

%%
%2. SYSTEM PREPARATION
%==========================================================================
%Exciton Basis, cm-1
% [coe_mat Ee]=eig(H);
% Ee=diag(Ee);                       %column vector of excitons
% Jcoe = zeros(ne, ne, ne, ne, ne);  %overlap functions
% for i = 1 : ne
%     for j = 1 : ne
%         for k = 1 : ne
%             for l = 1 : ne
%                 Jcoe(:,i,j,k,l)=coe_mat(:,i).*coe_mat(:,j).*coe_mat(:,k).*coe_mat(:,l);
%             end
%         end
%     end
% end

% %Excitons for unit in fs-1
% W=Ee*unitf;


%%
%3. Rate Calculation (preparing a rate matrix)
%==========================================================================
%Modified Redfield Rate
%Energy Transfer Rate and Pure-Dephasing Rate, output in fs-1
%tic;
[pR c_PD R_PD] = Rate_t_LMR(W, Jcoe, Reg, Redg, Reddg, Img, Imdg, Imddg, dtg, dtt);
%toc;

% R=ppval(pR, 10000); % Markovian Rate
% figure;bar3(R*1e3); % unit conversion to ps-1
% title('Exciton Transfer Rate');
% xlabel('Donor eigen state');
% ylabel('Acceptor eigen state');
% zlabel('Rate, ps^-^1');
%%
%NOTE: non-Markovian R_PD(t) = c_PD*Redg(t);
pRedg=spline(0:dtt:t_max, Redg);
%rPD =@(t) c_PD*ppval(pRedg,t);
toc

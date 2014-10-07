function R = DimerApproxEET4MP( r, theta )

tic;
R = zeros(3,1);

dt = 1; % time step for calculating rates
Tmax = 2e4; % maximum time for calculating rates
t = dt: dt: Tmax; % time sequence for calculating rates

[ H12, H34, E, MixAngle, J ] = Hamiltonian(r, theta);
% Prepare intermediate parameters for calculating rates
% H12: Hamiltonian of donor subsystem
% H34: Hamiltonian of acceptor subsystem
% E: eigen energies [E2, E3]
% MixAngle: mixing angle [theta12, theta34]
% J: effective coupling constants between |E2> and |E3>

% calculate R21 by CMRT
[ pR ] = LMR_Para4Dimer( H12, Tmax);
disR = ppval(pR,t); % dissipation rate
% R21 = mean(disR(1,2,round(0.5*Tmax):Tmax));
R(1) = mean(disR(1,2,round(0.5*Tmax):Tmax));

% calculate R43 by CMRT
[ pR ] = LMR_Para4Dimer( H34, Tmax);
disR = ppval(pR,t); % dissipation rate
% R43 = mean(disR(1,2,round(0.5*Tmax):Tmax));
R(2) = mean(disR(1,2,round(0.5*Tmax):Tmax));

% calculate R32 by generalized Forster
[ pR ] = FT_Para(E, MixAngle, Tmax);
disR = 2*(J^2)*ppval(pR,t); % dissipation rate
% R32 = mean(disR(round(0.5*Tmax):Tmax));
R(3) = mean(disR(round(0.5*Tmax):Tmax));

toc
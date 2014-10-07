tic;

RR = (10.8)*1e-10; % range of r
THETA = (0:0.01:0.5)*pi; % range of theta

dt = 1; % time step for calculating rates
Tmax = 2e4; % maximum time for calculating rates
t = dt: dt: Tmax; % time sequence for calculating rates

R21 = zeros(length(RR),length(THETA));
R43 = zeros(length(RR),length(THETA));
R32 = zeros(length(RR),length(THETA));

for j = 1: length(RR)
    for k = 1: length(THETA)
        
        r = RR(j);
        theta = THETA(k);
        
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
        R21(j,k) = mean(disR(1,2,round(0.5*Tmax):Tmax));
        
        % calculate R43 by CMRT
        [ pR ] = LMR_Para4Dimer( H34, Tmax);
        disR = ppval(pR,t); % dissipation rate
        R43(j,k) = mean(disR(1,2,round(0.5*Tmax):Tmax));
        
        % calculate R32 by generalized Forster
        [ pR ] = FT_Para(E, MixAngle, Tmax);
        disR = 2*(J^2)*ppval(pR,t); % dissipation rate
        R32(j,k) = mean(disR(round(0.5*Tmax):Tmax));
        
        k
    end;
end;

save('R21R43R32vsTheta.mat','R21','R43','R32');

toc

% calculate effective transfer rate
for j = 1: length(THETA)
    Reff(j)=R21(j)*R32(j)*R43(j)/(R21(j)*R32(j)+R32(j)*R43(j)+R21(j)*R43(j));
end;

plot(THETA,Reff)
tic;

RR = (6:0.2:14)*1e-10;
THETA = 0*pi;

dt = 1;
Tmax = 2e3;
t = dt: dt: Tmax;

for j = 1: length(RR)
    for k = 1: length(THETA)
        
        r = RR(j);
        theta = THETA(k);
        
        [ H12, H34, E, MixAngle, J ] = Hamiltonian(r, theta);
        
        [ pR ] = FT_Para(E, MixAngle, Tmax);
        disR = (J^2)*ppval(pR,t); % dissipation rate
        R32(j,k) = mean(disR(round(0.5*Tmax):Tmax));
        
    end;
end;

%plot(t,disR);
save('R32.mat','R32');

toc
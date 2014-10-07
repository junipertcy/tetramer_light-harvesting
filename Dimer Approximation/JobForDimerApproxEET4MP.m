r = (9: 0.2: 14)*1e-10;
theta = (0: 49)*pi/100;

nR = length(r);
nTheta = length(theta);

R = zeros(nR,nTheta,3);

matlabpool(7);

for k = 1 : nR
    parfor j = 1: nTheta
        
        R(k,j,:) = DimerApproxEET4MP(r(k), theta(j));
        
    end;
end;

matlabpool close;

R = E*1e3; % unit from fs^-1 to ps^-1
save('ThreeRvsRTheta.mat','R','r','theta'),
function [ H12, H34, Epsilon, MixAngle, Jeff ] = Hamiltonian( r , theta )

ne = 4; % dimension of H
H = zeros(ne);
H(1,1) = 13000;
H(2,2) = 12900;
H(3,3) = 12300;
H(4,4) = 12200;

R = 40*10^(-10); % length between D1 and A2, in units of m 
mu = 7.75*3.34*10^(-30); % electric dipole, in units of C*m

%Unit transfer factor of cm-1 to fs-1
unitf = 2*pi*3e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caution: In DipoleCoupling(mu1,mu2,theta1,theta2,r), theta1,% 
%          theta2 are angles between vec(mu)_j and vec(r),    %
%          while theta here is angle between two donors'      %
%          electric dipole moments and horizontal line.       %
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

H12 = H(1:2,1:2);
H34 = H(3:4,3:4);

E = eig(H12);
Epsilon(1) = E(1);
E = eig(H34);
Epsilon(2) = E(2);

MixAngle(1) = atan(2*H(1,2)/(H(1,1)-H(2,2)));
MixAngle(2) = atan(2*H(3,4)/(H(3,3)-H(4,4)));

J = (H(3,1)*sin(MixAngle(1)/2)-H(3,2)*cos(MixAngle(1)/2))*cos(MixAngle(2)/2)...
    +(H(4,1)*sin(MixAngle(1)/2)-H(4,2)*cos(MixAngle(1)/2))*sin(MixAngle(2)/2);
Jeff = J*unitf;
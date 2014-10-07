function [ J ] = DipoleCoupling( mu1 , mu2 , theta1 , theta2 , r )
% DipoleCoupling computes dipole coupling strength, in units of cm^-1
% input parameters:
% mu1: electric dipole moment vec(mu)_1
% mu2: electric dipole moment vec(mu)_2
% r: vec(r) displacement between two electic dipole moments
% theta1: angle between vec(mu)_1 and vec(r)
% theta2: angle between vec(mu)_2 and vec(r)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caution: theta1, theta2 are angles between vec(mu)_j        %
%          and vec(r)                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = 9*10^(9); % constant 1/(4*pi*epsilon_0)
hc = 6.626*3*10^(-24); % energy unit from J to cm^-1
J = C*(mu1*mu2*cos(theta1-theta2)-3*mu1*cos(theta1)*mu2*cos(theta2))/(r^3)/hc;
% electric diople coupling, in units of cm^-1
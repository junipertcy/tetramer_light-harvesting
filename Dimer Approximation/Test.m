RR = 11*1e-10;
THETA = 0.3*pi/100;
lambda = 35;
wc = 50;

ne = 4; % dimension of H
H = zeros(ne);
H(1,1) = 13000;
H(2,2) = 12900;
H(3,3) = 12300;
H(4,4) = 12200;

R = 40*10^(-10); % length between D1 and A2, in units of m
mu = 7.75*3.34*10^(-30); % electric dipole, in units of C*m
unitf = 2*pi*3e-2; % cm^-1 to ps^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caution: In DipoleCoupling(mu1,mu2,theta1,theta2,r), theta1,%
%          theta2 are angles between vec(mu)_j and vec(r),    %
%          while theta here is angle between two donors'      %
%          electric dipole moments and horizontal line.       %                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for j = 20 : length(RR)
%     for k = 1 : length(THETA)
        r = RR(1);
        theta = THETA(1);
        
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
        theta12 = atan(2*H12(1,2)/(H12(1,1)-H12(2,2)));
        
        H34 = H(3:4,3:4);
        theta34 = atan(2*H34(1,2)/(H34(1,1)-H34(2,2)));
        
        e12 = eig(H12);
        e34 = eig(H34);
        
        e2 = e12(1);
        e3 = e34(2);
        
        K2 = (2-0.5*(sin(theta12)^2+sin(theta34)^2))*pi*lambda/(2*wc);
        
        J = -H(3,2)*cos(theta12/2)*cos(theta34/2)
        % J = (H(3,1)*sin(theta12/2)-H(3,2)*cos(theta12/2))*cos(theta34/2)...
        %     +(H(4,1)*sin(theta12/2)-H(4,2)*cos(theta12/2))*sin(theta34/2)
        
        Rate(1,1) = 2*J^2*sin(K2)/(e2-e3);
        
%     end;
% end;

% imagesc(Rate);
% colorbar('fontsize', 14, 'ytick',[min(min(Rate)) max(max(Rate))])
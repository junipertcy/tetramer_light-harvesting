nR = 41;
nTheta = 50;

Reff=zeros(nR,nTheta);

for j=1:nR
    for k=1:nTheta
        R1 = R(1,j,k)*R(3,j,k)/(R(1,j,k)+R(3,j,k));
        R2 = R(2,j,k)*R(3,j,k)/(R(2,j,k)+R(3,j,k));
        Reff(j,k) = min(R1,R2);
    end;
end;
Reff(Reff<0)=1e-6;

Reff = Reff*1e3;

[C,h] = contourf(Reff);
colormap autumn;

set(gca,'YTick',[1 6 11 16 21 26 31 36 41],'XTick',[1 11 21 31 41 50]);
set(gca,'YTickLabel',{'6','7','8','9','10','11','12','13','14'},...
    'XTickLabel',{'0','0.1','0.2','0.3','0.4','0.5'},'fontsize', 14);

xlabel('\theta/\pi','fontsize',14);
ylabel('r (A)','fontsize',14);
colorbar('fontsize', 14);
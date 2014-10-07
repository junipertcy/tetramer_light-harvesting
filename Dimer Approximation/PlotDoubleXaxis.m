% PlotDoubleXaxis plot two sets of data for two different x axes

% plot R_eff vs r
hl1 = line(RR1,Reff1,'Color','k','LineStyle','--','LineWidth',2);
xlim([min(RR1) 13])
ylim([0 0.2])
ax1 = gca;
xlabel('r (A)','FontSize',16);
ylabel('R_{eff} (ps^{-1})','FontSize',16);
set(gca,'XTick',[6 7 8 9 10 11 12 13],...
    'YTick',[0 0.05 0.10 0.15 0.20],'FontSize',16)
set(gca,'XTickLabel',{'6','7','8','9','10','11','12','13'},...
    'YTickLabel',{'0','0.05','0.10','0.15','0.20'},'FontSize',16)

% plot R_eff vs theta
ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top',...
           'YAxisLocation','right','Color','none','XColor','r');
hl2 = line(THETA2,Reff2,'Color','r','LineStyle','-','LineWidth',2,'Parent',ax2);
xlabel('\theta/\pi','FontSize',16,'Position',[0.25 0.205]);
set(gca,'YTick',[0 0.05 0.10 0.15 0.20],'XTick',[0 0.1 0.2 0.3 0.4 0.5],...
    'FontSize',16)
set(gca,'YTickLabel',{'','','','',''},...
    'XTickLabel',{'0','0.1','0.2','0.3','0.4','0.5'},'FontSize',16)
xlim([min(THETA2) max(THETA2)])
ylim([0 0.2])
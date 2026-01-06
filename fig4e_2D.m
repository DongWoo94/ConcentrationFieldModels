clear all; close all

% define constants
C.D = 150 * 60; % [=] um2/min
C.RNAP = 1; %
C.RNase = 1; % [=] U/ml; 0.5 U/mL  1e6 uL = 1L; 1000L = 1m3; 1m3 = 1e18 um3; 1uL = 1e9 um3; 1mL = 1e12 um3
C.k_Dss = 0.036 * 60; % [=] (U/ml)-1 min-1 ; 20 (U/uL)-1 s-1  1uL = 1e9 um3; 1e9 uL = 1m3; 2.4e-2 s-1
C.L = 25; % um

lambda = sqrt(C.k_Dss*C.RNase/C.D); % um
C.Template = 25*[1/besselk(0,100*lambda) 1/besselk(0,200*lambda) 1/besselk(0,300*lambda)] * 1e-9 * 1e-15 * 0.02 * 60; % [=] mol/um3/min; 400nM initial; rejection rate 40%;  1mol/L = 1e-15 mol/um3    1L = 1e15 um3
% 150, 1000, 5500 nM/min [1 : 6.67 : 36.7]
1/besselk(0,lambda*100)
C.Template*1e24

% plotting
fig = figure();
newcolors = {'k','r','b'};
colororder(newcolors)
        axe = axes();
        fig.Units = "inches";
        axe.LineWidth = 1.5;
        axe.FontSize  = 16;
        axe.NextPlot  = "add";
        axe.Box       = "on";
        axe.XLim = [0 500]; %[-100, 400];
        % axe.YLim = [0 5]; %[-1, 6];
        axe.XLabel.String = "r / μm";
        axe.XAxis.FontWeight = 'bold';
        axe.XLabel.FontSize = 24;
        axe.XLabel.FontWeight = 'normal';
        axe.YLabel.String = "C / nM";
        axe.YAxis.FontWeight = 'bold';
        axe.YLabel.FontSize = 24;
        axe.YLabel.FontWeight = 'normal';
        axe.XMinorTick = "on";
        axe.YMinorTick = "on";
        axe.TickLength = [0.03 0.05];
        legend('FontSize',12, 'Location','northeast')
        legend('boxoff')
        axis square
        % xref = xline(C.L,":");
        % xref.LineWidth = 1.5;
        % xref.FontSize = 12;
        % xref.LabelOrientation = 'aligned';
        % xref.LabelHorizontalAlignment = 'left';
        % xref.LabelVerticalAlignment = 'bottom';
        hold on
        legend off
        set(axe,'yscale','log')

% Solve series ODE
xl = linspace(0.001,C.L,31);
xu = linspace(C.L, C.L+10000 ,10000);
xmesh = [xl xu];

K11 = besselk(1,lambda*C.L);
K01 = besselk(0,lambda*C.L);
I11 = besseli(1,lambda*C.L);
I01 = besseli(0,lambda*C.L);

for i = 1:length(C.Template)
A = C.Template(i)*C.RNAP/C.k_Dss/C.RNase; % /min-1

% analytic solution
c1 = A*(1 - lambda*C.L*K11*besseli(0,lambda*xl));
c2 = A*lambda*C.L*I11*besselk(0,lambda*xu);
c = [c1 c2];

% num = plot(sol.x,sol.y(1,:)*1e24,'ko');
hold on
ana = plot(xmesh,c*1e24,"LineWidth",1.5);
% title('Solution of 2D Diffusion Equation with bvp5c');

% legend([num ana],{'numerical solution','analytic solution'}) 
end

f = gcf;
exportgraphics(f,'fig4e_2D.jpg','Resolution',300)

clear all; close all

% define constants
C.D = 1; % [=] um2/min
C.RNAP = 1; %
C.k_p = 1; % [=] min-1 Txn rate constant (fitted)
C.RNase = 1; % [=] U/um3; 0.5 U/mL  1e6 uL = 1L; 1000L = 1m3; 1m3 = 1e18 um3; 1uL = 1e9 um3; 1mL = 1e12 um3
kc = logspace(-1,1,3); % [=] (U/um3)-1 min-1 ; 20 (U/uL)-1 s-1  1uL = 1e9 um3; 1e9 uL = 1m3; 2.4e-2 s-1
C.Template = 1; % [=] mol/um3; 400nM initial; rejection rate 40%;  1mol/L = 1e-15 mol/um3    1L = 1e15 um3
C.L = 1; % um

% variables 
C.Lstar = logspace(-2,2,101); 

fig = figure();
newcolors = {'k','r','b'};
colororder(newcolors)
   axe = gca;
        fig.Units = "inches";
        axe.LineWidth = 1.5;
        axe.FontSize  = 16;
        axe.NextPlot  = "add";
        axe.Box       = "on";
        axe.XLim = [C.Lstar(1) C.Lstar(end)]; %[-100, 400]; 
        axe.YLim = [C.Lstar(1) 1]; %[-1, 6];
        axe.XLabel.String = "R'";
        axe.XAxis.FontWeight = 'bold';
        axe.XLabel.FontSize = 24;
        axe.XLabel.FontWeight = 'normal';
        axe.YLabel.String = "C'_{surface}";
        axe.YAxis.FontWeight = 'bold';
        axe.YLabel.FontSize = 24;
        axe.YLabel.FontWeight = 'normal';
        axe.XMinorTick = "on";
        axe.YMinorTick = "on";
        axe.TickLength = [0.03 0.05];
        axe.XTick = logspace(-2, 2, 3);
        axe.YTick = logspace(-2, 0, 3);
        legend('FontSize',12, 'Location','southeast')
        legend('boxoff')
        axis square
        set(axe, 'XScale', 'log')
        set(axe, 'YScale', 'log')

% analytic solution
c_max_1d = (1 - exp(-C.Lstar).*cosh(C.Lstar));
p1 = plot(C.Lstar , c_max_1d,'LineWidth',1);

% analytic solution
K11 = besselk(1,C.Lstar);
I01 = besseli(0,C.Lstar);
c_max_2d = (1 - C.Lstar.*K11.*I01); % nM
p2 = plot(C.Lstar, c_max_2d,'LineWidth',1);

% analytic solution
c_max_3d = (1 - (C.Lstar + 1).*exp(-C.Lstar).*sinh(C.Lstar)./C.Lstar);
p3 = plot(C.Lstar, c_max_3d,'LineWidth',1);

legend([p1 p2 p3],{'1D','2D','3D'}) 
exportgraphics(fig,'fig3a.jpg','Resolution',300)




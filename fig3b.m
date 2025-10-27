clear all; close all

% define constants
C.D = 1; % [=] um2/min
C.RNAP = 1; %
C.k_p = 1; % [=] min-1 Txn rate constant (fitted)
C.RNase = 1; % [=] U/um3; 0.5 U/mL  1e6 uL = 1L; 1000L = 1m3; 1m3 = 1e18 um3; 1uL = 1e9 um3; 1mL = 1e12 um3
kc = logspace(-2,0,3); % [=] (U/um3)-1 min-1 ; 20 (U/uL)-1 s-1  1uL = 1e9 um3; 1e9 uL = 1m3; 2.4e-2 s-1
C.Template = 1; % [=] mol/um3; 400nM initial; rejection rate 40%;  1mol/L = 1e-15 mol/um3    1L = 1e15 um3
C.L = 1; % um

% variables 
C.Lstar = logspace(-2,2,101); 

fig = figure();
newcolors = {'k','k','r','r','b','b'};
colororder(newcolors)
   axe = gca;
        fig.Units = "inches";
        axe.LineWidth = 1.5;
        axe.FontSize  = 16;
        axe.NextPlot  = "add";
        axe.Box       = "on";
        axe.XLim = [C.Lstar(1) C.Lstar(end)]; %[-100, 400]; 
        % axe.YLim = [C.Lstar(1) C.Lstar(end)]; %[-1, 6];
        axe.XLabel.String = "R'";
        axe.XAxis.FontWeight = 'bold';
        axe.XLabel.FontSize = 24;
        axe.XLabel.FontWeight = 'normal';
        axe.YLabel.String = "HWHM / R";
        axe.YAxis.FontWeight = 'bold';
        axe.YLabel.FontSize = 24;
        axe.YLabel.FontWeight = 'normal';
        axe.XMinorTick = "on";
        axe.YMinorTick = "on";
        axe.TickLength = [0.03 0.05];
        axe.XTick = logspace(-2, 2, 3);
        % axe.YTick = logspace(-3, 1, 3);
        legend('FontSize',12, 'Location','northeast')
        legend('boxoff')
        axis square
        set(axe, 'XScale', 'log')
        set(axe, 'YScale', 'log')

% analytic solution
for i = 1:length(C.Lstar)
    HWHM(i) = fzero(@(HW) 2*sinh(C.Lstar(i)).*exp(-HW.*C.Lstar(i)) - (1 - exp(-C.Lstar(i)).*cosh(C.Lstar(i))), [0.0001 1000]);
    K11 = besselk(1,C.Lstar(i));
    K01 = besselk(0,C.Lstar(i));
    I11 = besseli(1,C.Lstar(i));
    I01 = besseli(0,C.Lstar(i));
    HWHM2(i) = fzero(@(HW) 2*besselk(0,HW*C.Lstar(i)).*(C.Lstar(i).*I11) - (1 - C.Lstar(i).*K11.*I01), [0.0001 1000]); 
    HWHM3(i) = fzero(@(HW) (1 - (C.Lstar(i) + 1).*exp(-C.Lstar(i)).*sinh(C.Lstar(i))./C.Lstar(i)) - 2*(C.Lstar(i).*cosh(C.Lstar(i)) - sinh(C.Lstar(i))).*exp(-C.Lstar(i)*HW)./HW./C.Lstar(i), [0.0001 10000]);
end

p1 = plot(axe,C.Lstar, HWHM,'LineWidth',1);
p2 = plot(axe,C.Lstar, HWHM2,'r','LineWidth',1);
p3 = plot(axe,C.Lstar, HWHM3,'b','LineWidth',1);

% p4 = plot(axe,C.Lstar,C.Lstar,':k','LineWidth',1);
legend([p1 p2 p3],{'1D','2D','3D'}) 

exportgraphics(fig,'fig3b.jpg','Resolution',300)




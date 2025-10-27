clear all; close all

% define constants
C.D = 1; % [=] um2/min
C.RNAP = 1; %
C.k_p = 1; % [=] min-1 Txn rate constant (fitted)
C.RNase = 1; % [=] U/um3; 0.5 U/mL  1e6 uL = 1L; 1000L = 1m3; 1m3 = 1e18 um3; 1uL = 1e9 um3; 1mL = 1e12 um3
kc = logspace(-1,1,3); % [=] (U/um3)-1 min-1 ; 20 (U/uL)-1 s-1  1uL = 1e9 um3; 1e9 uL = 1m3; 2.4e-2 s-1
C.Template = 1; % [=] mol/um3; 400nM initial; rejection rate 40%;  1mol/L = 1e-15 mol/um3    1L = 1e15 um3
C.L = 1; % um

% plotting
fig = figure();
newcolors = {'k','k','r','r','b','b'};
colororder(newcolors)
        axe = axes();
        fig.Units = "inches";
        axe.LineWidth = 1.5;
        axe.FontSize  = 16;
        axe.NextPlot  = "add";
        axe.Box       = "on";
        axe.XLim = [0 5]; %[-100, 400];
        axe.YLim = [1e-6 1]; %[-1, 6];
        axe.XLabel.String = "r / R";
        axe.XAxis.FontWeight = 'bold';
        axe.XLabel.FontSize = 24;
        axe.XLabel.FontWeight = 'normal';
        axe.YLabel.String = "C'";
        axe.YAxis.FontWeight = 'bold';
        axe.YLabel.FontSize = 24;
        axe.YLabel.FontWeight = 'normal';
        axe.XMinorTick = "on";
        axe.YMinorTick = "on";
        axe.TickLength = [0.03 0.05];
        legend('FontSize',12, 'Location','northeast')
        legend('boxoff')
        axis square
        % xref = xline(C.L,":",'Boundary');
        % xref.LineWidth = 1.5;
        % xref.FontSize = 12;
        % xref.LabelOrientation = 'aligned';
        % xref.LabelHorizontalAlignment = 'right';
        % xref.LabelVerticalAlignment = 'bottom';
hold on
    set(axe, 'YScale', 'log')

%% kc = 0.1
C.k_Dss = 0.1;
lambda = sqrt(C.k_Dss*C.RNase/C.D); % um
A = C.k_p*C.Template*C.RNAP/C.k_Dss/C.RNase; % /min-1
lambda
% Solve series ODE
xl = linspace(0.001,C.L,11);
xu = linspace(C.L, C.L+100 ,1000);
xmesh = [xl xu];
yinit = [C.Template; C.Template];
sol = bvpinit(xmesh,yinit);
sol = bvp4c(@(x,y,region) PolarDiffusion(x,y,region,C), @bc, sol);

K11 = besselk(1,lambda*C.L);
K01 = besselk(0,lambda*C.L);
I11 = besseli(1,lambda*C.L);
I01 = besseli(0,lambda*C.L);

% analytic solution
c1 = (1 - lambda*C.L*K11*besseli(0,lambda*xl));
c2 = lambda*C.L*I11*besselk(0,lambda*xu);
c = [c1 c2];

% num1 = plot(sol.x,sol.y(1,:)/A,'o');
ana1 = plot(xmesh,c,'k-','LineWidth',1);

%% kc = 1
C.k_Dss = 1;
lambda = sqrt(C.k_Dss*C.RNase/C.D); % um
A = C.k_p*C.Template*C.RNAP/C.k_Dss/C.RNase; % /min-1
lambda
% Solve series ODE
xl = linspace(0.001,C.L,11);
xu = linspace(C.L, C.L+100 ,1000);
xmesh = [xl xu];
yinit = [C.Template; C.Template]*1e-3;
sol = bvpinit(xmesh,yinit);
sol = bvp4c(@(x,y,region) PolarDiffusion(x,y,region,C), @bc, sol);

K11 = besselk(1,lambda*C.L);
K01 = besselk(0,lambda*C.L);
I11 = besseli(1,lambda*C.L);
I01 = besseli(0,lambda*C.L);

% analytic solution
c1 = (1 - lambda*C.L*K11*besseli(0,lambda*xl));
c2 = lambda*C.L*I11*besselk(0,lambda*xu);
c = [c1 c2];

% num2 = plot(sol.x,sol.y(1,:)/A,'o');
ana2 = plot(xmesh,c,'r-','LineWidth',1);

%% kc = 10
C.k_Dss = 10;
lambda = sqrt(C.k_Dss*C.RNase/C.D); % um
A = C.k_p*C.Template*C.RNAP/C.k_Dss/C.RNase; % /min-1
lambda
% Solve series ODE
xl = linspace(0.001,C.L,11);
xu = linspace(C.L, C.L+100 ,1000);
xmesh = [xl xu];
yinit = [C.Template; C.Template]*1e-3;
sol = bvpinit(xmesh,yinit);
sol = bvp4c(@(x,y,region) PolarDiffusion(x,y,region,C), @bc, sol);

K11 = besselk(1,lambda*C.L);
K01 = besselk(0,lambda*C.L);
I11 = besseli(1,lambda*C.L);
I01 = besseli(0,lambda*C.L);

% analytic solution
c1 = (1 - lambda*C.L*K11*besseli(0,lambda*xl));
c2 = lambda*C.L*I11*besselk(0,lambda*xu);
c = [c1 c2];

% num3 = plot(sol.x,sol.y(1,:)/A,'o');
ana3 = plot(xmesh,c,'b-','LineWidth',1);

patch = fill([0 1 1 0], [c2(end) c2(end) 1 1], 'k');
set(patch,'edgecolor', 'none');
set(patch,'FaceAlpha', 0.10);

legend([ana1 ana2 ana3],{"R' = 0.3","R' = 1","R' = 3"}) 

f = gcf;
exportgraphics(f,'fig2e.jpg','Resolution',300)

function dydx = PolarDiffusion(x,y,region,C)
% constants
lambda = sqrt(C.k_Dss*C.RNase/C.D);
A = C.Template*C.k_p/C.D;

dydx = zeros(2,1);
dydx(1) = y(2);

switch region
    case 1    % x in [0 1]
        dydx(2) = -1/x*y(2) + lambda^2*y(1) - A;
    case 2    % x in [1 inf]
        dydx(2) =-1/x*y(2) + lambda^2*y(1);
end
end

function res = bc(YL,YR)
res = [YL(2,1) - 0  % y'(0) = 0
       YR(1,1) - YL(1,2)  % Continuity of c(x) at x=1
       YR(2,1) - YL(2,2)  % Continuity of dc(x)/dx at x=1
       YR(2,end) - 0];    % c(end) = 0
end

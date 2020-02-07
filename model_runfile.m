clear; close all; clc;
%%% runfile to solve the ODEs for the VQ matching model
%%% Copyright University of Michigan, Andrew (Drew) Marquis, Feb 11th 2020
%%% If you find a bug, have a question/comment/suggestion, email me at admarqui@umich.edu

%%% Settings - BE VERY CAREFUL OF THE SAVE INITIAL CONDITIONS AND SAVE SIMULATIONS SETTINGS!!!
%%%          - It is extremely easy to accidentally overwrite files if you are not careful
%%%          - General Rules: 
%%%          - (1) Save one thing at a time
%%%          - (2) Avoid saving new initial conditions if at all possible. Similarly, I would advise always use the saved initial conditions unless you're trying to do something new/weird
%%%          - Loose Warning: The model is pretty stiff and can take a long time to equilibrate and reach steady state (~3 min of smulation time).
%%%            I would avoid running the model with all of the vessel flows, pressures, and O2 concentrations initialized at zero - relaxation to 
%%%            steady state typically requires these variabless to transiently become negative and may requite fidling in the VQ_RHS file 
%%%            and/or not enforcing the non-negative condition to make sure ode15s does not get stuck. I do recomend enforcing non-negative 
%%%            for O2 concentrations when running simulations to make predictions (it's numerically/mathematically possible to go below zero, non-neg enforces conservation of mass)

%what to simulate
FLAG  = 1; %set FLAG to select which network you want to simulate (1 == reuleaux, 2 == quarter circ, 3 == half circ)
fFodP = 1; %"fixed Flow or change in Pressure (delta P)" (0 = fixed flow, 1 = fixed primary inlet-outlet pressure drop)
iIC   = 1; %boolean to initialize model with saved initial conditions (0 == no, 1 == yes)
HPV   = 0; %"triple" boolean to use HPV control (0 == no, 1 == yes, 2 == uniform vasoconstriction)

%what to save
sIC   = 0; %boolean to save end points as initial conditions (0 == no, 1 == yes)
sSS   = 0; %boolean to save results from steady state simulation (0 == no, 1 == yes)
sHPV  = 0; %boolean to save results from HPV equilibration simulation (0 == no, 1 == yes)
sUVC  = 0; %boolean to save results from global uniform vasoconstriction simulation (0 == no, 1 == yes)

%%% loading info parameters and connectivity info
[CON, Par, PDATA, LOOK, X0] = load_network_par(FLAG);
Nt   = CON.Nt; %CON is a struct that stores network geometry and connectivity info
Nseg = CON.Nseg; segL = CON.segL; L    = CON.L; TI   = CON.TI; LtT  = CON.LtT;
Net  = CON.Net; COORD = CON.COORD; MAP = CON.MAP; Node = CON.Node;
Element = CON.Element; A = CON.A;

%%% loading optimized parameters from BELUGA
if FLAG == 1
    load('reuleaux_par_opt_logscaled_par.mat','optPAR')
elseif FLAG == 2
    load('quarter_circ_par_opt_logscaled_par.mat','optPAR')
elseif FLAG == 3
    load('half_circ__par_opt_logscaled_par.mat','optPAR')
end
optPAR = exp(optPAR);
Par.DO   = optPAR(1);
Par.Tmax = optPAR(2);
Par.Tmin = optPAR(3);

if iIC == 1 %yes, we want to use saved initial conditions
    if FLAG == 1
        if fFodP == 0
            load('IC_reuleaux_fF.mat','XIC') %loading initial conditions
        elseif fFodP == 1
            load('IC_reuleaux_fdP.mat','XIC') %loading initial conditions
        end
    elseif FLAG == 2
        if fFodP == 0
            load('IC_quarter_circ_fF.mat','XIC')
        elseif fFodP == 1
            load('IC_quarter_circ_fdP.mat','XIC') %loading initial conditions
        end
    elseif FLAG == 3
        if fFodP == 0
            load('IC_half_circ_fF.mat','XIC')
        elseif fFodP == 1
            load('IC_half_circ_fdP.mat','XIC') %loading initial conditions
        end
    end
    T0 = X0((end-Nseg+1):end);
    X0 = [XIC(1:end-Nseg) T0'];
end

%%% load data from hypoxia-hemodyanics experiments
addpath Cheng_etal_raw_data/
load('Cheng_DATA.mat') %loading the data we want to fit the model to
%%% apply DATA to boundary conditions
IND = [2 8]; 
% Par.Pair = 75; %if we want to make the inspired air hypoxic (10% O2)
i = 1; %i = 1 --> normoxia, i = 2 --> hypoxia

Par.pdrive = DATA.SPAP(IND(i));
Par.Qvent  = DATA.VENT(IND(i));
Par.Qco    = DATA.CO(IND(i))/2;

%%% solving the system of ODEs
TOL     = 1e-6; %error tolerance
NNIND = [(2*Nseg+Nt+1):(4*Nseg+2*Nt) (4*Nseg+2*Nt+1):(5*Nseg+3*Nt)]; %indices for ODE states (pressure/volume and oxygen mass) that need the non-negative condition enforced
% NNIND = [];
options = odeset('RelTol',TOL,'AbsTol',TOL,'NonNegative',NNIND); %ode solving options
tic;
[t,X] = ode15s(@VQ_RHS, [0 100*Par.R/(2*pi)], X0, options, Par, CON, PDATA, LOOK, HPV,fFodP); 
rt_mat = toc;
disp(rt_mat)

XIC = X(end,:); %variable end points to use as IC's
if sIC == 1     %yes we want to save initial conditions
    if FLAG == 1
        if fFodP == 0
            save('IC_reuleaux_fF.mat','XIC') %loading initial conditions
        elseif fFodP == 1
            save('IC_reuleaux_fdP.mat','XIC') %loading initial conditions
        end
    elseif FLAG == 2
        if fFodP == 0
            save('IC_quarter_circ_fF.mat','XIC')
        elseif fFodP == 1
            save('IC_quarter_circ_fdP.mat','XIC') %loading initial conditions
        end
    elseif FLAG == 3
        if fFodP == 0
            save('IC_half_circ_fF.mat','XIC')
        elseif fFodP == 1
            save('IC_half_circ_fdP.mat','XIC') %loading initial conditions
        end
    end
end

%%% post-ODE solving organization and computing other relevant quantities
qinA  = X(:,1:Nseg);                 %inlet flow (arterial)
qinV  = X(:,(Nseg+1):(2*Nseg));      %inlet flow (venous)
qinC  = X(:,(2*Nseg+1):(2*Nseg+Nt)); %inlet flow (capillary)
qin   = [qinA qinV qinC];

VvA = X(:,(2*Nseg+Nt+1):(3*Nseg+Nt));   %vascular (intraluminal) pressure (arterial)
VvV = X(:,(3*Nseg+Nt+1):(4*Nseg+Nt));   %vascular (intraluminal) pressure (venous)
VvC = X(:,(4*Nseg+Nt+1):(4*Nseg+2*Nt)); %vascular (intraluminal) pressure (capillary)
Vv  = [VvA VvV VvC];

T    = X(:,5*Nseg+4*Nt+3:end); %HPV tone
Cv_A = Par.Cv_A';
Cv_A = Cv_A./T;
Cv   = Par.Cv;  %compliance parameters
Cvc  = Par.Cvc';

%%% back to processing model solutions
PvA = VvA./Cv_A;%arterial pressures
PvV = VvV./Cv';  %venous pressures
PvC = VvC./Cvc; %capillary pressures
Pv  = [PvA PvV PvC];

VA = sum(VvA,2); %total arterial volume
VV = sum(VvV,2); %total venous volume
VC = sum(VvC,2); %total capillary volume
VTOT = VA+VV+VC; %total blood volume

Ppl   = Par.Pscale*sin(Par.R*t);      %alveolar pressure
Pout  = GET_POUT(t, Pv, qin, T, Ppl, Par, CON); %outlet pressures

PoutA = Pout(:,1:Nseg);
PoutV = Pout(:,(Nseg+1):(2*Nseg));
PoutC = Pout(:,(2*Nseg+1):end);

CO    = X(:,(4*Nseg+2*Nt+1):(4*Nseg+3*Nt)); %capillary oxygen concentrations
COven = X(:,(4*Nseg+3*Nt+1):(5*Nseg+3*Nt)); %venous oxygen concentrations
COTOT = COven(:,CON.xind);                  %total oxygen concentration leaving pulm circulation
POTOT = interp1(LOOK.Clookup,LOOK.Plookup,COTOT); %oxygen tension leaving pulm circulation
POven = interp1(LOOK.Clookup,LOOK.Plookup,COven);
n = Par.n; P50 = Par.P50; %parameters for Hb saturation
SOTOT = POTOT.^n./(POTOT.^n+P50.^n);           %total oxygen saturation leaving pulm circulation
PO    = interp1(LOOK.Clookup,LOOK.Plookup,CO); %capillary oxygen tension
SO2   = PO.^n./(PO.^n+P50.^n);                 %capillary oxygen saturation

VQ    = abs(Par.qVent./qinC); %VQ ratios

VENT     = X(:,(5*Nseg+3*Nt+1):5*Nseg+4*Nt+2); %ventilation states
Pairway1 = VENT(:,1);                %partial pressure of O2 in the airways (mmHg)
Pairway2 = VENT(:,2);                %partial pressure of O2 in the airways (mmHg)
Malv     = VENT(:,3:end);            %mass of O2 in the alveoli (mol)

qsi   = Par.qVent;                           %flow magnitude per each alveoli
V0i   = Par.V0*CON.Ft;                       %volume per each alveoli
Valv  = V0i'+qsi'*(1-cos(Par.R.*t'))./Par.R; %alveolar volumes - this could be formulated as another ODE state, but we're going with the analytic solution
Valv  = Valv';
Palv  = Par.gamma*Malv./Valv;      %alveolar oxygen partial pressure

%%%% checking mass conservation
MO    = CO.*VvC;
MOven = COven.*VvV;
MOtot = MOven(:,CON.xind);

%%% empirical HPV control stuff
%some parameters
Phpv    = Par.Phpv; %(mmHg of O2)
lambda  = Par.lambda; %(mm)
Tmin    = Par.Tmin;
Tmax    = Par.Tmax;

S = exp(-Palv/Phpv); %HPV signal from each terminal capillary
Tinf_u = zeros(length(t),Nseg); %unitless tone for each vessel segment
for i = 1:Nseg
   if i <= Nt
      Tinf_u(:,i) = S(:,i).*exp(-L(i)/2/lambda);%./(1 + S(:,i).*exp(-L(i)/2/lambda));
   else
       vi = i-Nt; %non-terminal vessel segment index
       sS = zeros(length(t),1); %preallocate signal sum
       for k = 1:length(TI{vi})
          sS = sS + S(:,k).*exp(-LtT{vi}(k)/lambda); 
       end
       Tinf_u(:,i) = sS;%sS./(1+sS);
   end
end
Tinf = Tinf_u*(Tmax-Tmin)+Tmin; %unitless tone, but scaled with Tmin and Tmax

%%% O2 flux
Aref = 6.714326703426993e7*1e-6; %total area of the perfused domain (mm^2)
JO2  = Par.alpha*Par.DO*(Palv-PO)./(Aref*CON.Ft);

%%% average value of quantities of interest
tL   = 2*Par.R/(2*pi); %"time length" - amount of time to integrate over - 2 respiratory cycles
tIND = find(t<(t(end)-tL), 1, 'last' ); %index of time vector to start integrating over
tINT = t(tIND:end); %relevant portion of the time vector

nIND = length(t(tIND:end));

PoutC_1 = PoutC(1,:); %capillary pressure
PoutC_2 = trapz(tINT,PoutC(tIND:end,:))/tL; %trapezoid rule approximation of integral, then divided by tL (whats up fundamental theorem of calculus)

PvC_1 = PvC(1,:);
PvC_2 = trapz(tINT,PvC(tIND:end,:))/tL;

qinC_1  = qinC(1,:); %capillary flow
qinC_2  = trapz(tINT,qinC(tIND:end,:))/tL;

JO2_1   = trapz(tINT, JO2(1:nIND,:))/tL; %O2 flux
JO2_2   = trapz(tINT, JO2(tIND:end,:))/tL;

VQ_1    = VQ(1,:); %VQ ratio
VQ_2    = trapz(tINT,VQ(tIND:end,:))/tL;

qinA_2  = trapz(tINT,qinA(tIND:end,:))/tL;
qinV_2  = trapz(tINT,qinV(tIND:end,:))/tL;
VvA_2   = trapz(tINT,VvA(tIND:end,:))/tL;
VvC_2   = trapz(tINT,VvC(tIND:end,:))/tL;
VvV_2   = trapz(tINT,VvV(tIND:end,:))/tL;

%%% RBC transit time distributions (before and after HPV activation)
[~, ~, hb, tb, mub] = TT_deconvolution(Nt, CON.NETsegA,qinA(1,:),qinC(1,:),qinV(1,:),VvA(1,:),VvC(1,:),VvV(1,:));
[~, ~, ha, ta, mua] = TT_deconvolution(Nt, CON.NETsegA,qinA_2,qinC_2,qinV_2,VvA_2,VvC_2,VvV_2);

%%% saving steady state simulation results
if sSS == 1
    if FLAG == 1
        save('reuleaux_SS_results.mat')
    elseif FLAG == 2
        save('quarter_circ_SS_results.mat')
    elseif FLAG == 3
        save('half_circ_SS_results.mat')
    end
end

%%% saving pulsatile HPV equilibration simulation results
if sHPV == 1
    if FLAG == 1
        save('reuleaux_HPV_results.mat')
    elseif FLAG == 2
        save('quarter_circ_HPV_results.mat')
    elseif FLAG == 3
        save('half_circ_HPV_results.mat')
    end
end

%%% saving pulsatile HPV equilibration simulation results
if sUVC == 1
    if FLAG == 1
        save('reuleaux_UVC_results.mat')
    elseif FLAG == 2
        save('quarter_circ_UVC_results.mat')
    elseif FLAG == 3
        save('half_circ_UVC_results.mat')
    end
end

%%% plotting
figure; %flows
subplot(3,1,1)
plot(t,qinA,'b','linewidth',2) %arterial flows 
hold on
plot(t,qinV,'r','linewidth',2) %venous flows 
hold on
plot(t,qinC,'k','linewidth',2) %capillary flows 
set(gca,'fontsize',18)
ylabel('q_{in} (ml/s)')
xlabel('Time (s)')
grid on

% figure; %luminal pressures
subplot(3,1,2)
plot(t,PoutA, 'b', 'linewidth',2)%arterial pressures
hold on
plot(t,PoutV, 'r', 'linewidth',2)%venous pressures
hold on
plot(t,PoutC, 'k', 'linewidth',2)%capilllary pressures 
set(gca,'fontsize',18)
ylabel('P_{out} (mmHg)')
xlabel('Time (s)')
grid on

% figure; %volumes
subplot(3,1,3)
plot(t,VA, 'b', 'linewidth',2)%arterial volumes
hold on
plot(t,VV, 'r', 'linewidth',2)%venous volumes
hold on
plot(t,VC, 'k', 'linewidth',2)%capilllary volumes 
set(gca,'fontsize',18)
ylabel('V_v (mmHg)')
xlabel('Time (s)')
grid on

% figure; %wall pressures
% plot(t,PvA, 'b', 'linewidth',2)%arterial pressures
% hold on
% plot(t,PvV, 'r', 'linewidth',2)%venous pressures
% hold on
% plot(t,PvC, 'k', 'linewidth',2)%capilllary pressures 
% set(gca,'fontsize',18)
% ylabel('P_v (mmHg)')
% xlabel('Time (s)')
% grid on

figure; %airway O2 tension
subplot(2,1,1)
plot(t, Palv,'color',[0,0,0]+0.75,'linewidth',2)
hold on
plot(t, Pairway1, 'k:', t, Pairway2, 'k-.', 'linewidth',2)
set(gca,'fontsize',18)
xlabel('Time (s)')
ylabel('Airway Oxygen(mmHg)')
% legend('Pairway1','Pairway2')
grid on

% figure;
subplot(2,1,2)
plot(t,PO,'color',[0,0,0]+0.75,'linewidth',2)
hold on
plot(t,POTOT,'k','linewidth',2)
set(gca,'fontsize',18)
ylabel('Capillary Oxygen (mmHg)')
xlabel('Time (s)')
% ylim([0 160])
grid on


% figure; %alveolar O2 mass
% plot(t,Malv,'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time (s)')
% ylabel('Alveolar O2 (mol)')
% grid on

% figure;
% plot(t,Valv,'linewidth',2)
% set(gca,'fontsize',18)
% xlabel('Time (s)')
% ylabel('Alveolar volume (ml)')
% grid on

figure;
subplot(1,2,1)
plot(qinC(end,:)./mean(qinC(end,:)),PO(end,:),'o','linewidth',2)
set(gca,'fontsize',18)
ylabel('Venous Oxygen Tension (mmHg)')
xlabel('Capillary Flow/Mean Capillary Flow')
grid on

subplot(1,2,2)
semilogx(VQ_2,PO(end,:),'o','linewidth',2)
set(gca,'fontsize',18)
xlabel('V/Q Ratio')
xlim([1e-2 1e3])
% ylabel('Oxygen Tension (mmHg)')
grid on

figure;
semilogx(VQ_2,qinC_2,'o','linewidth',2,'markersize',5)
hold on
semilogx(VQ_2,Par.qVent,'o','linewidth',2,'markersize',5)
set(gca,'fontsize',18)
xlabel('V/Q ratio')
ylabel('Flow (ml/s)')
legend('Blood','Air')
grid on

figure;
plot(t,S)
set(gca,'fontsize',18)
xlabel('Time (s)')
ylabel('HPV signal (S) (unitless)')
grid on

figure;
plot(t,Tinf)
set(gca,'fontsize',18)
xlabel('Time (s)')
ylabel('HPV Tone (Tinf) (unitless)')
grid on

figure;
plot(t,T)
set(gca,'fontsize',18)
xlabel('Time (s)')
ylabel('HPV Tone (T) (unitless)')
grid on

figure;
plot(tb,hb,':','linewidth',1.5)
hold on
plot(ta,ha,'linewidth',1.5)
set(gca,'fontsize',18)
xlabel('Time (s)')
ylabel('h(t)')
legend('w/o HPV','w/ HPV')
grid on
xlim([0 2])

figure;
plot(tb/mub,hb,':','linewidth',1.5)
hold on
plot(ta/mua,ha,'linewidth',1.5)
set(gca,'fontsize',18)
xlabel('Time/mean TT (unitless)')
ylabel('h(t)')
legend('w/o HPV','w/ HPV')
grid on
xlim([0 2])

fh =  findobj('type','figure'); %figure handles
fN = length(fh); %number of figures made

% PLOT_PERFUSION_ZONES(Nt, MAP, Element, Node, PoutC(1,:), 'Capillary Pressure (mmHg)', 0, fN+1)
% hold on
% gplot2(A, COORD,'k-','linewidth',2)
% set(gca,'xtick',[],'ytick',[],'linewidth',2)
% axis equal
% axis off
% 
% PLOT_PERFUSION_ZONES(Nt, MAP, Element, Node, PoutC(end,:), 'Capillary Pressure (mmHg)', 0, fN+2)
% hold on
% gplot2(A, COORD,'k-','linewidth',2)
% set(gca,'xtick',[],'ytick',[],'linewidth',2)
% axis equal
% axis off

space = 5000;
COORDA = COORD; COORDA(:,1) = COORD(:,1)+space;
COORDB = COORD; COORDB(:,1) = COORD(:,1)-space;
COMPARE_PERFUSION_ZONES(Nt, MAP, Element, Node,  PoutC_1, PoutC_2, 'Capillary Pressure (mmHg)', 5000, 0, fN+1)
hold on
gplot2(A, COORDA,'k-','linewidth',2)
hold on
gplot2(A, COORDB,'k-','linewidth',2)
set(gca,'xtick',[],'ytick',[],'linewidth',2)
axis equal
axis off


% COMPARE_PERFUSION_ZONES(Nt, MAP, Element, Node,  PvC_1, PvC_2, 'Capillary Wall Pressure (mmHg)', 5000, 0, fN+100)
% hold on
% gplot2(A, COORDA,'k-','linewidth',2)
% hold on
% gplot2(A, COORDB,'k-','linewidth',2)
% set(gca,'xtick',[],'ytick',[],'linewidth',2)
% axis equal
% axis off

% PLOT_PERFUSION_ZONES(Nt, MAP, Element, Node, log10(qinC(1,:)), 'Capillary Flow (ml/s)', 1, fN+3)
% hold on
% gplot2(A, COORD,'k-','linewidth',2)
% set(gca,'xtick',[],'ytick',[])
% axis equal
% axis off
% 
% PLOT_PERFUSION_ZONES(Nt, MAP, Element, Node, log10(qinC(end,:)), 'Capillary Flow (ml/s)', 1, fN+4)
% hold on
% gplot2(A, COORD,'k-','linewidth',2)
% set(gca,'xtick',[],'ytick',[])
% axis equal
% axis off

COMPARE_PERFUSION_ZONES(Nt, MAP, Element, Node,  log10(qinC_1), log10(qinC_2), 'Capillary Flow (ml/s)', 5000, 1, fN+2)
hold on
gplot2(A, COORDA,'k-','linewidth',2)
hold on
gplot2(A, COORDB,'k-','linewidth',2)
set(gca,'xtick',[],'ytick',[],'linewidth',2)
axis equal
axis off

COMPARE_PERFUSION_ZONES(Nt, MAP, Element, Node,  log10(VQ_1), log10(VQ_2), 'VQ ratio (unitless)', 5000, 1, fN+3)
hold on
gplot2(A, COORDA,'k-','linewidth',2)
hold on
gplot2(A, COORDB,'k-','linewidth',2)
set(gca,'xtick',[],'ytick',[],'linewidth',2)
axis equal
axis off

% PLOT_PERFUSION_ZONES(Nt, MAP, Element, Node, JO2b, 'Alveolar-Capillary O_2 flux (mol/(s*mm^2))', 0, fN+5)
% hold on
% gplot2(A, COORD,'k-','linewidth',2)
% set(gca,'xtick',[],'ytick',[],'linewidth',2)
% axis equal
% axis off
% 
% PLOT_PERFUSION_ZONES(Nt, MAP, Element, Node, JO2a, 'Alveolar-Capillary O_2 flux (mol/(s*mm^2))', 0, fN+6)
% hold on
% gplot2(A, COORD,'k-','linewidth',2)
% set(gca,'xtick',[],'ytick',[],'linewidth',2)
% axis equal
% axis off

COMPARE_PERFUSION_ZONES(Nt, MAP, Element, Node,  log10(JO2_1), log10(JO2_2), 'Alveolar-Capillary O_2 flux (mol/(s*mm^2))', 5000, 1, fN+4)
hold on
gplot2(A, COORDA,'k-','linewidth',2)
hold on
gplot2(A, COORDB,'k-','linewidth',2)
set(gca,'xtick',[],'ytick',[],'linewidth',2)
axis equal
axis off

% PLOT_PERFUSION_ZONES(Nt, MAP, Element, Node, PO(1,:), 'Capillary O_2 (mmHg))', 0, fN+7)
% hold on
% gplot2(A, COORD,'k-','linewidth',2)
% set(gca,'xtick',[],'ytick',[],'linewidth',2)
% axis equal
% axis off
% 
% PLOT_PERFUSION_ZONES(Nt, MAP, Element, Node, PO(end,:), 'Capillary O_2 (mmHg))', 0, fN+8)
% hold on
% gplot2(A, COORD,'k-','linewidth',2)
% set(gca,'xtick',[],'ytick',[],'linewidth',2)
% axis equal
% axis off

% PLOT_PERFUSION_ZONES(Nt, MAP, Element, Node, MO(end,:), 'Mass of O_2)', 0, fN+40)
% hold on
% gplot2(A, COORD,'k-','linewidth',2)
% set(gca,'xtick',[],'ytick',[],'linewidth',2)
% axis equal
% axis off

% COMPARE_PERFUSION_ZONES(Nt, MAP, Element, Node,  PO(1,:), PO(end,:), 'Capillary O_2 (mmHg))', 5000, 0, fN+4)
% hold on
% gplot2(A, COORDA,'k-','linewidth',2)
% hold on
% gplot2(A, COORDB,'k-','linewidth',2)
% set(gca,'xtick',[],'ytick',[],'linewidth',2)
% ylim([min(Node(:,2)) max(Node(:,2))])
% axis equal
% axis off

NETWORK_PLOT(segL, COORD, Net, T(end,:), 'Tone from HPV (unitless)', 0, fN+5)
axis equal
axis off

% NETWORK_PLOT(segL, COORD, Net, PoutA(1,:), 'P_{out} (mmHg)', fN+10)
% NETWORK_PLOT(segL, COORD, Net, PoutA(end,:), 'P_{out} (mmHg)', fN+11)

COMPARE_NETWORK_PLOT(segL, COORD, Net, PoutA(1,:), PoutA(end,:), space, 'P_{out} (mmHg)', 0, fN+6)
set(gca,'xtick',[],'ytick',[],'linewidth',2)
axis equal
axis off

% NETWORK_PLOT(segL, COORD, Net, PoutV(1,:), 'P_{out} (mmHg)', fN+12)
% NETWORK_PLOT(segL, COORD, Net, PoutV(end,:), 'P_{out} (mmHg)', fN+13)

COMPARE_NETWORK_PLOT(segL, COORD, Net, PoutV(1,:), PoutV(end,:), space, 'P_{out} (mmHg)', 0, fN+7)
set(gca,'xtick',[],'ytick',[],'linewidth',2)
axis equal
axis off

[fp1, xp1] = ksdensity(PoutC_1);
[fp2, xp2] = ksdensity(PoutC_2);
figure;
plot(xp1,fp1,xp2,fp2,'linewidth',2)
set(gca,'fontsize',20)
legend('w/o HPV', 'w/ HPV')
xlabel('Capillary Pressure (mmHg)')
grid on

[fq1, xq1] = ksdensity(log10(qinC_1));
[fq2, xq2] = ksdensity(log10(qinC_2));
figure;
plot(xq1,fq1,xq2,fq2,'linewidth',2)
set(gca,'fontsize',20)
legend('w/o HPV', 'w/ HPV')
xlabel('log_{10}[Capillary Flow]')
grid on

[fVQ1, xVQ1] = ksdensity(log10(VQ_1));
[fVQ2, xVQ2] = ksdensity(log10(VQ_2));
figure;
plot(xVQ1,fVQ1,xVQ2,fVQ2,'linewidth',2)
set(gca,'fontsize',20)
legend('w/o HPV', 'w/ HPV')
xlabel('log_{10}[V/Q ratio]')
grid on

[fJO1, xJO1] = ksdensity(log10(JO2_1));
[fJO2, xJO2] = ksdensity(log10(JO2_2));
figure;
plot(xJO1,fJO1,xJO2,fJO2,'linewidth',2)
set(gca,'fontsize',20)
legend('w/o HPV', 'w/ HPV')
xlabel('log_{10}[Alveolar-Capillary O_2 flux]')
grid on

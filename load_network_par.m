function [CON, Par, PDATA, LOOK, X0] = load_network_par(FLAG)
%%% This script computes and returns the parameters for the model in VQ_RHS

%%% load anatomical geometry infomration made in "network_run"
addpath VascNetwork
if FLAG == 1
    load('NETWORK_reuleaux.mat','i3','SO','segL','COORD','Net', 'Ft', 'NETseg', 'NETsegV', 'TI','LtT','X', 'MAP', 'Element', 'Node', 'A') %i3 is the list of terminal nodes - all we need for this file at least the number of outlets ( length(i3) )
elseif FLAG == 2
    load('NETWORK_quarter_circ.mat','i3','SO','segL','COORD','Net', 'Ft', 'NETseg', 'NETsegV', 'TI','LtT','X', 'MAP', 'Element', 'Node', 'A') %i3 is the list of terminal nodes - all we need for this file at least the number of outlets ( length(i3) )
elseif FLAG == 3
    load('NETWORK_half_circ.mat','i3','SO','segL','COORD','Net', 'Ft', 'NETseg', 'NETsegV', 'TI','LtT','X', 'MAP', 'Element', 'Node', 'A') %i3 is the list of terminal nodes - all we need for this file at least the number of outlets ( length(i3) )
end

%NETseg = [segmentIndex radius Length parentIndex daughter1Index daughter2Index]
%the primary vessel has parentIndex == 0. Terminal vessels have daughterIndex == 0.5
%X is the index of the primary vessel segment (feed artery)
%Ft is the idealized terminal flow for each termianl segment

load('RAT031518_HD.mat','DATA') %hemodynamic data from a rat experiment - using the pulmonary artery waveform to drive the model
i1 = find(DATA.t == 0.25);
i2 = find(DATA.t == 0.454);
PDATA = [DATA.t(i1:i2)'-DATA.t(i1) DATA.Ppa(i1:i2)]; %pulmonary artery pressure data - the first column is corresponding time points - going to use interp1 to find the appropriate pressure value in the RHS of the ODE

%%% some initial info
Nt    = length(i3); %number of terminal segments
xind  = X;          %index of the primary artery/vein segment 
Pla   = -5;          %fixed left atrial pressure (mmHg)
pdrive = 20;        %fixed inlet pressure (mmHg)

%%% anatomical info from vascular network
r = NETseg(:,2)*1e-3;  %radii for each segment (mm)
L = NETseg(:,3)*1e-3;  %length for each segment (mm)
NETseg(:,2:3) = [];    %deleting radius and length from this struct to save space/time when called in the ODE solver
Nseg = size(NETseg,1); %total number of segments (for just the arterial or venous network)

%%% invariant physical information
hct = 0.45;      %hematocrit
E   = 1e4;       %youngs modulus for a blood vessel (mmHg)
rho = 1060;      %blood density (kg/m^3)
f   = 0.01;      %time constant of vessel relaxation (sec)

%%% parameters for each vessel segment (for topologically equivalent arterial and venous segments)
%ntertance
Iv  = 1e3*rho.*L./(pi.*r.^2); %Intertance for each vessel segment (Kg/m^4 or Pa*s^2/m^3)
Iv  = Iv./(133.3e6);          %(mmHg*s^2/ml) converting from Pa to mmHg

%hydraulic resistance
eta = 1.2e-3*VISCOSITY(2.*r*1e3,hct)/133.322; %apparent viscosity in each vessel segment: 1.2e-3 - viscosity of plasma (Pa*s), 133.322 - conversion factor for Pa to mmHg
Rv  = 1e3*8.*eta.*L./(pi*r.^4);               %resistance for each vessel segment (mmHg*s/ml) - multiplying by 1e-3 to convert mm^3 to ml

%compliance
a  = 0.2802; b = -0.5053; c = 0.1324; d = -0.01114; %parameters for empirical vessel wall thickness as function of luminal radius
h  = r.*(a.*exp(b.*r)+c.*exp(d.*r)); %thickness of vessel wall (mm) as an empirical function of luminal radius 
Cv = (2*pi/E).*(r.^3).*L./h;         %compliance for each vessel segment (ml/mmHg)

%%% saving arterial parameters
Rv_A = Rv;
Cv_A = Cv;

%global venous vessel tone
Tone = 1/2;
Rv   = Rv*Tone;
Cv   = Cv/Tone;
Rd   = f./Cv; %viscous dampening of the vessel wall (mmHg*s/ml)

%left atrial resistance
Qco    = 80/60/2;          %average cardiac ouput of a rat (80 ml/min is normal - converted to ml/s)- used to scale normalized flow to physical units
dP     = 2;                %pressure drop across capillary - a made up number
Rla    = dP./(Qco);        %resistance leading to left atrium

%%% capillary parameters
VBC = 2*1e-3*pi*sum(r.^2.*L); %volume of blood in pulmonary capilaries (ml)

Ivc = mean(Iv(1:Nt))*ones(Nt,1); %assuming capillaries have the same intertance - mean of terminal arteriole intertances
Rvc = ones(Nt,1)*dP./(Qco*Ft');  %assumed pressure drop over a priori terminal flows
Cvc = (ones(Nt,1).*Ft').*(0.7*VBC)./dP; %stressed volume of capillaries over assumed pressure drop. assuming stressed volume is 70% of total capilary volume
Rdc = f./Cvc;            %estimate dampening resistance in capillaries the same as art/vein segments.

%%% gas exchange parameters
alpha   = 1.3e-6;  % O2 solubility by Henry's law(M/mmHg)
beta    = 3.08e-5; % CO2 solubility by Hency's law (M/mmHg)
CHb     = 0.021;   % Hb binding site conc (mol/L of RBC's)
Hct     = 0.40;    % hematocrit (unitless)
C0      = CHb*Hct; % blood oxygen binding capacity (mol/L)
n       = 2.7;     % Hill exponent
P50     = 30;      % half-max saturation of Hb
DO      = 1;       % (apparent) O2 diffusion coefficient (ml/s)

% O2 look up table (concetration to partial pressure conversion)
Plookup = 0:1:10000; %look up table
Clookup = alpha*Plookup + C0*(Plookup.^n)./(Plookup.^n + P50.^n);

LOOK.Plookup = Plookup;
LOOK.Clookup = Clookup;

PAO2 = 100; %resting systemic artial O2 partial pressure (mmHg)
PVO2 = 45;  %resting systemic venous O2 partial pressure (mmHg)
CAO2 = alpha*PAO2 + C0*(PAO2.^n)./(PAO2.^n + P50.^n);   %converted to mol/L
CVO2 = alpha*PVO2 + C0*(PVO2.^n)./(PVO2.^n + P50.^n);   %coverted to mol/L

MO = Qco*(CAO2 - CVO2); % O2 metabolism coefficient (mol/s) - from Fick's Law

%%% Ventilation parameters (also gas exchange)
Mass    = .350;% (kg)
RR      = 53.46*Mass^(-0.2599); %Respiration Rate (bpm)- numbers from rovent ventilator
R       = 2*pi*RR/60;           %breath per sec*(2pi) - parameter for model
Qvent   = Qco;                  %setting minute ventilation 
Pscale  = 0.75;                 %alveolar pressure amplitude (mmHg)
Vairway = 0.6;   %airway volume of lung
gamma   = 16800; %mmHg * L / mol (convert M to mmHg, for gas)

%%% empirical HPV control parameters
if FLAG == 1
    %network 1
    Phpv   = 80/log(2);  % HPV O2 sensitivity (mmHg of O2)
    lambda = 1e2;  % vessel length conduction constant (um)
    tau    = 1.5;  % time averaging constant for tone (s)
    Tmin   = 0.35; % minimum Tone
    Tmax   = 1.2; % maximum Tone
elseif FLAG == 2
    %network 2
    Phpv   = 80/log(2);  % HPV O2 sensitivity (mmHg of O2)
    lambda = 1e2;  % vessel length conduction constant (mm)
    tau    = 1.5;   % time averaging constant for tone (s)
    Tmin   = 0.35; % minimum Tone
    Tmax   = 1; % maximum Tone
elseif FLAG == 3
    %network 3
    Phpv   = 80/log(2);  % HPV O2 sensitivity (mmHg of O2)
    lambda = 1e2;  % vessel length conduction constant (mm)
    tau    = 1.5;   % time averaging constant for tone (s)
    Tmin   = 0.45; % minimum Tone
    Tmax   = 1.3; % maximum Tone
end

%%% set up initial conditions
%hemodynamics
Q0 = zeros(2*Nseg+Nt,1);        %initial flows (per segment)
P0 = pdrive.*ones(2*Nseg+Nt,1); %initial pressures (per segment)
Vv0 = pdrive.*[Cv_A/Tone; Cv; Cvc];

%oxygen transport in blood vessels
CO     = CVO2.*ones(Nt,1);   %initial oxygen concentration for capillaries (per capillary segment)
COven0 = CVO2.*ones(Nseg,1); %initial oxygen concentration for venous segments
cOin   = CVO2;

%oxygen transport in airways
RFt = Ft'; %if we want to use the perfusion zone size distribution

Pair     = 150;            %oxygen partial pressure in normal air
V0       = 3;              %total alveolar volume initial condition
qVent    = Qvent.*RFt';    %air flow through each perfusion zone
Pairway1 = 142;
Pairway2 = 140;
Malvi    = 0.001*ones(Nt,1).*RFt;

%%%
T0 = ones(Nseg,1)/2; %HPV tone IC

%initial conditions stacked into a column vector
X0 = [Q0; Vv0; CO; COven0; Pairway1; Pairway2; Malvi; T0]; 

%%% assemble structs
%parameter struct
Par.Pla = Pla;
Par.pdrive = pdrive;
Par.Qco = Qco;

Par.Iv  = Iv;
Par.Rv  = Rv;
Par.VA0 = 1e-3*pi*r.^2.*L*12; %basic definition of volume for each artery (ml)
Par.VV0 = 1e-3*pi*r.^2.*L; %basic definition of volume for each vein (ml)
Par.VC0 = VBC*Ft'; %total volume of capillary blood dispered through the organ by Ft
Par.Cv  = Cv;
Par.Rd  = Rd;
Par.Rla = Rla;

Par.Rv_A = Rv_A;
Par.Cv_A = Cv_A;

Par.Ivc = Ivc;
Par.Rvc = Rvc;
Par.Cvc = Cvc;
Par.Rdc = Rdc;

Par.alpha = alpha;
Par.beta  = beta;
Par.CHb   = CHb;
Par.Hct   = Hct;
Par.C0    = C0;
Par.n     = n;
Par.P50   = P50;
Par.DO    = DO;
Par.MO    = MO;
Par.cOin  = cOin;

Par.gamma   = gamma;
Par.Vairway = Vairway;
Par.Qvent   = Qvent;
Par.Pscale  = Pscale;
Par.R       = R;
Par.Pair    = Pair;
Par.V0      = V0;
Par.qVent   = qVent;

Par.Phpv   = Phpv;
Par.lambda = lambda;
Par.tau    = tau;
Par.Tmin   = Tmin;
Par.Tmax   = Tmax;

%network geometry and connectivity struct
CON.Nt      = Nt;
CON.Nseg    = Nseg;
CON.NETsegA = NETseg;
CON.NETsegV = NETsegV;
CON.segL    = segL;
CON.xind    = xind;
CON.MAP     = MAP;
CON.Element = Element;
CON.Node    = Node;
CON.Net     = Net;
CON.A       = A;
CON.Ft      = Ft;
CON.COORD   = COORD;
CON.i3      = i3;
CON.SO      = SO;
CON.r       = r;
CON.L       = L;
CON.TI      = TI;
CON.LtT     = LtT;

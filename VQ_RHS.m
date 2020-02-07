function dX = VQ_RHS(t, X, Par, CON, PDATA, LOOK, HPV, fFodP)

%%% connectivity info
Nt      = CON.Nt;
Nseg    = CON.Nseg;
NETsegA = CON.NETsegA;
NETsegV = CON.NETsegV;
xind    = CON.xind;
L       = CON.L;
LtT     = CON.LtT;
Ft      = CON.Ft;
RFt     = Ft';

%%%look up table
Clookup = LOOK.Clookup;
Plookup = LOOK.Plookup;

%%%initialize parameters
Pla     = Par.Pla;    %venous/left atrium outlet BC parameters
pdrive  = Par.pdrive;

Iv      = Par.Iv;     %aterial and venous segment parameters
Rv      = Par.Rv; 
Cv      = Par.Cv;
Rd      = Par.Rd;
Rla     = Par.Rla;
Ivc     = Par.Ivc;    %capillary segment parameters
Rvc     = Par.Rvc;
Cvc     = Par.Cvc;
Rdc     = Par.Rdc;
alpha   = Par.alpha;  %gas exchange parameters
gamma   = Par.gamma;
DO      = Par.DO;
MO      = Par.MO;
cOin    = Par.cOin;
R       = Par.R;      %airway parameters
Pscale  = Par.Pscale; 
Qvent   = Par.Qvent;
Qco     = Par.Qco;
Vairway = Par.Vairway;
Pair    = Par.Pair;
V0      = Par.V0;
Rv_A    = Par.Rv_A;   %venous resitance
Cv_A    = Par.Cv_A;   %venous compliance
Phpv    = Par.Phpv;   %HPV control parameters
lambda  = Par.lambda;
tau     = Par.tau;
Tmin    = Par.Tmin;
Tmax    = Par.Tmax;

%%% State Variables (for each segment)
qinA  = X(1:Nseg);                 %inlet flow (arterial)
qinV  = X((Nseg+1):(2*Nseg));      %inlet flow (venous)
qinC  = X((2*Nseg+1):(2*Nseg+Nt)); %inlet flow (capillary)
qin   = [qinA; qinV; qinC];        %stacking inflows into one vector

VvA = X((2*Nseg+Nt+1):(3*Nseg+Nt));   %vascular (intraluminal) pressure (arterial)
VvV = X((3*Nseg+Nt+1):(4*Nseg+Nt));   %vascular (intraluminal) pressure (venous)
VvC = X((4*Nseg+Nt+1):(4*Nseg+2*Nt)); %vascular (intraluminal) pressure (capillary)

CO    = X((4*Nseg+2*Nt+1):(4*Nseg+3*Nt)); %capillary oxygen concentration
COven = X((4*Nseg+3*Nt+1):(5*Nseg+3*Nt)); %venous segment oxygen concentration

%%% ventilation driving function
Qair  = Qvent*sin(R*t); %assumed driving function
qi    = Qair*RFt;        %driving airflow for each alveolar compartment
V0i   = V0*RFt;          %volumes for each alveolar compartment

Valv = V0i+Qvent.*RFt*(1-cos(R*t))/R; %alveolar volumes - this could be formulated as another ODE state, but we're going with the analytic solution
Ppl  = Pscale*sin(R*t);                %alveolar pressure

%%% ventilation mechanics/O2 transport states
VENT     = X((5*Nseg+3*Nt+1):(5*Nseg+4*Nt+2)); %ventilation states 
Pairway1 = VENT(1);                %airway 1 - interacts with atmosphere
Pairway2 = VENT(2);                %airway 2 - big ol' mixing chamber
Malv     = VENT(3:end);            %mass of O2 in the alveoli (mol)
PalvO    = gamma*Malv./Valv;       %alveolar oxygen tension (mmHg)

f = 0.7; %fraction of airway volume to be in the mixing compartment - this is made up and don't know if it's reasonable. Probably worth doing some local sensitivity analysis on it
Vairway1 = Vairway*(1-f);
Vairway2 = Vairway*f;

%%% HPV tone states
if HPV == 2 % if we want to globally vasoconstrict the system
    T = 0.8*ones(Nseg,1);
else
    T = X(5*Nseg+4*Nt+3:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% math to set up RHS equations
%%% myogenic constriction (volume-dependent compliance)
VA0 = Par.VA0; %reference baseline arterial volumes
VC0 = Par.VC0; %reference baseline capillary volumes
VV0 = Par.VV0; %reference baseline venous volumes

%%% apply HPV tone to resistances and compliances
Rv_A = Rv_A.*T;
Cv_A = Cv_A./T;
Rd_A = 0.01./Cv_A;

%%% stuff for non-linear resistors
PvA = VvA./Cv_A; %instantaneous arterial volumes
PvC = VvC./Cvc;  %instantaneous capillary volumes
PvV = VvV./Cv;   %instantaneous venous volumes
Pv = [PvA; PvV; PvC];

Rva = Rv_A.*(VA0./VvA).^2;
Rvc = Rvc.*(VC0./VvC).^2; %dynamic resistances
Rvv = Rv.*(VV0./VvV).^2;

%%%%%% vascular 'per segment' quantities
%%% outflows
%arterial qout
qoutA             = zeros(Nseg,1);                        %preallocate outlet flow for each arterial segment
qoutA(1:Nt)       = qinC(1:Nt);                           %outflow of terminal arterioles is the inflow to capillaries - these are the flows that will determine RBC transit times and gas exchange
qoutA((Nt+1):end) = sum(qinA(NETsegA((Nt+1):end,3:4)),2); %sum flow through daughters for outlet flow - using NETsegA data struct to index daughter vessels

%capillary qout
qoutC = qinV(1:Nt); %outflow of capillaries is the inflow to terminal venules

%venous qout
qoutV = zeros(Nseg,1); %preallocate outlet flow for each venous segment
if fFodP == 0
    qoutV(xind) = Qco; %fixed venous Flow
elseif fFodP == 1
    qoutV(xind) = (PvV(xind) - Ppl + Pla + Rd(xind)*qinV(xind))./(Rla+Rd(xind)); %fixed outlet pressure drop
end
qoutV(1:(xind-1)) = (PvV(NETsegV(1:(xind-1),1))-PvV(NETsegV(1:(xind-1),2))+...
    Rd(NETsegV(1:(xind-1),2)).*(qinV(NETsegV(1:(xind-1),3))-qinV(NETsegV(1:(xind-1),2)))+...
    Rd(NETsegV(1:(xind-1),1)).*qinV(NETsegV(1:(xind-1),1)))./(Rd(NETsegV(1:(xind-1),1))+Rd(NETsegV(1:(xind-1),2)));
qoutV((xind+1):end) = (PvV(NETsegV((xind+1):end,1))-PvV(NETsegV((xind+1):end,2))+...
    Rd(NETsegV((xind+1):end,2)).*(qinV(NETsegV((xind+1):end,3))-qinV(NETsegV((xind+1):end,2)))+...
    Rd(NETsegV((xind+1):end,1)).*qinV(NETsegV((xind+1):end,1)))./(Rd(NETsegV((xind+1):end,1))+Rd(NETsegV((xind+1):end,2)));

qout = [qoutA; qoutV; qoutC]; %stack outflows into one vector

%%% outlet pressures
Pout  = Pv + Ppl + [Rd_A; Rd; Rdc].*(qin - qout); %outlet pressures - from conservation of energy/force balance
PoutA = Pout(1:Nseg);       %arterial outlet pressures - need this to determine arterial and input pressures
PoutC = Pout(2*Nseg+1:end); %capillary outlet pressures

%%% inlet pressures
%arterial Pin
PinA               = zeros(Nseg,1); %preallocate input pressures
PinA(xind)         = pdrive; %fixed input pressure
% PinA(xind)         = interp1(PDATA(:,1),PDATA(:,2),mod(t,PDATA(end,1))); %assign feed pressure by interpolating from data
PinA(1:(xind-1))   = PoutA(NETsegA(1:(xind-1),2)); %Pin is the Pout in the immediately distal vessel segment - using NETseg to index the corresponding parent
PinA((xind+1):end) = PoutA(NETsegA((xind+1):end,2));

%capillary Pin
PinC = PoutA(1:Nt);

%venous Pin
PinV           = zeros(Nseg,1);
PinV(1:Nt)     = PoutC(1:Nt); %venules have outlet pressure of arterioles
% PinV1 = Ppl+PvV(NETsegV(Nt+1:end,4)) + ... %venous inlet pressure is based on a presure balance of partner vessels - using NETsegV to index upstream vessel(s)
%                  Rd(NETsegV(Nt+1:end,4)).*(qinV(NETsegV(Nt+1:end,4)) - qoutV(NETsegV(Nt+1:end,4))); %this is one of the equations used to solve for venous outflows
PinV2 = Ppl+PvV(NETsegV(Nt+1:end,5)) + ...
                 Rd(NETsegV(Nt+1:end,5)).*(qinV(NETsegV(Nt+1:end,5)) - qoutV(NETsegV(Nt+1:end,5))); 
PinV(Nt+1:end) = PinV2;%min([PinV1 PinV2],[],2);     

Pin = [PinA; PinV; PinC]; %stacking inlet pressures

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% evalutating RHS
%%% hemodynamics RHS
dQ  = (Pin - Pout - [Rva; Rvv; Rvc].*qin)./[Iv; Iv; Ivc]; %conservation of momentum (inflow state)                   
dVv = qin - qout;                                         %conservation of mass (vascular volume state)

%%% capillary oxygen RHS
PO  = interp1(Clookup,Plookup,CO); %capillary oxygen tension
COven_t = COven(1:Nt); %venous terminal O2 concentrations

% cOin = COven(xind) - MO/(2*qinV(xind));% setting inlet O2 oxygen concentration from Fick's Law
% we multiply the venous flow rate by 2, as this model only simulates one lung which receives half of the cardiac output

capi = qinC  >= 0; %logical indices for capillary inflows >= 0
capo = qoutC >= 0; capo_n = logical(1-capo);  %logical indices for capillary outflows >= 0 and < 0

JO2cap         = zeros(Nt,1); %preallocate capillary O2 flux (from sources other than exchange with alveoli)
JO2cap(capi)   = JO2cap(capi) + qinC(capi).*(cOin - CO(capi)); %O2 from arteriers
JO2cap(capo_n) = JO2cap(capo_n) - qoutC(capo_n).*(COven_t(capo_n) - CO(capo_n)); %O2 from venous backflow
dCO            = (JO2cap + alpha*DO*(PalvO - PO))./VvC; %all together with alv-cap O2 transport
% dCO            = JO2cap./VvC; %all together with alv-cap O2 transport

% dCO = (qinC.*(cOin-CO)+alpha*DO*(PalvO-PO))./VvC; %if no backflow, this works for the RHS

%%% venous oxygen RHS
veno    = qoutV >= 0;   % logical for positive venous outflow
veno_n  = qoutV < 0;%logical(1 - veno);     % logical for negative venous outflow
veno_nt = logical(veno_n(1:Nt)); % extracting logicals for just terminal venules
veno_s  = logical([zeros(Nt+1,1); veno(Nt+2:end)]); %extracting logicals for non-terminal venous segments

c = NETsegA(:,2); %indices up downsream segments c
a = NETsegA(:,3); %indices up upstream segments a
b = NETsegA(:,4); %indices of upstream segments b

JO2ven = zeros(Nseg,1);
disp(find(c(veno_n)<=0))
JO2ven(veno_n) = JO2ven(veno_n) - qoutV(veno_n).*(COven(c(veno_n)) - COven(veno_n)); %O2 from backflow
JO2ven(c(veno_s)) = JO2ven(c(veno_s)) + qoutV(a(c(veno_s))).*(COven(a(c(veno_s))) - COven(c(veno_s))); %O2 inflow from vein a
JO2ven(c(veno_s)) = JO2ven(c(veno_s)) + qoutV(b(c(veno_s))).*(COven(b(c(veno_s))) - COven(c(veno_s))); %o2 inflwo from vein b

JO2ven_t          = zeros(Nt,1); %O2 flux for terminal venules
JO2ven_t(capo)    = JO2ven_t(capo) + qoutC(capo).*(CO(capo)-COven_t(capo)); %O2 inflow from capillaries
JO2ven_t(veno_nt) = JO2ven_t(veno_nt) - qoutV(veno_nt).*(CO(veno_nt) - COven_t(veno_nt)); %O2 from backflow
JO2ven(1:Nt)      = JO2ven_t; %storing venous O2 fluxes into one vector

dCOven = JO2ven./VvV; %venous O2 concentration RHS

% dCOven = zeros(Nseg,1); %if no backflow, this works for the RHS
% dCOven(1:Nt) = qinV(1:Nt).*(CO - COven(1:Nt))./VvV(1:Nt);
% dCOven(Nt+1:end) = (qout(NETsegV(Nt+1:end,4)).*(COven(NETsegV(Nt+1:end,4))-COven(Nt+1:end))+qout(NETsegV(Nt+1:end,5)).*(COven(NETsegV(Nt+1:end,5))-COven(Nt+1:end)))./VvV(Nt+1:end);

%%% airway oxygen (ventilation) RHS - mass transport through airways and exchange at alv-cap boundary
dV = zeros(2+Nt,1);
if Qair > 0
  dV(1)     = (Qair/Vairway1)*(Pair - Pairway1);          %airway1 RHS
  dV(2)     = (Qair/Vairway2)*(Pairway1 - Pairway2);      %airway2 RHS
  dV(3:end) = qi.*Pairway2/gamma - alpha*DO*(PalvO - PO); %Malv RHS
else
  dV(1)     = -(Qair/Vairway1)*(Pairway2 - Pairway1);  %airway1 RHS
  dV(2)     = -sum(qi.*(PalvO - Pairway2))./Vairway2;  %airway2 RHS
  dV(3:end) = qi.*PalvO/gamma - alpha*DO*(PalvO - PO); %Malv RHS
end

%%% HPV RHS
S       = exp(-PalvO/Phpv);  %HPV signal from each terminal capillary
Tinf_u  = zeros(Nseg,1);  %unitless tone for each vessel segment [0 1]
for a = 1:Nseg            %iterate by segment
   if a <= Nt
      Tinf_u(a) = S(a)*exp(-L(a)/2/lambda);%/(1 + S(i)*exp(-L(i)/2/lambda)); %signals for terminal segments
   else
       vi = a-Nt; %non-terminal vessel segment index
       sS = 0;    %preallocate signal sum
       for k = 1:length(LtT{vi})
          sS = sS + S(k)*exp(-LtT{vi}(k)/lambda); %compute cumulative signal
       end
       Tinf_u(a) = sS;%sS/(1+sS); %convert signal into unitless tone
   end
end
Tinf = Tinf_u*(Tmax-Tmin)+Tmin; %unitless tone, but scaled with Tmin and Tmax [0.4 0.5]
dT   = (Tinf-T)./tau;            %HPV Tone RHS - this equation is conceptually similar to mechanical creep
if HPV == 0 || 2
    dT = 0*dT;
end

% disp(t)
% disp([t sum(VvA <0) sum(VvC <0) sum(VvV <0)])
%%% final RHS stacked into one monster vector
dX  = [dQ; dVv; dCO; dCOven; dV; dT];
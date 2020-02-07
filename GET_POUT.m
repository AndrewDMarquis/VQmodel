function Pout = GET_POUT(t, Pv, qin, T, Ppl, Par, CON)

N = length(t);
%%% connectivity info
Nt      = CON.Nt;
Nseg    = CON.Nseg;
NETsegA = CON.NETsegA;
NETsegV = CON.NETsegV;
xind    = CON.xind;

%%%initialize parameters
Pla  = Par.Pla;
Qco  = Par.Qco;
Cv_A = Par.Cv_A;
Rv_A = Par.Rv_A;

Rd      = Par.Rd';
Rla     = Par.Rla;
Rdc     = Par.Rdc';

%%% time-varying resistance and compliance
Rv_A = Rv_A'.*T;
Cv_A = Cv_A'./T;
Rd_A = 0.01./Cv_A;


%%%
qinA = qin(:,1:Nseg);
qinV = qin(:,(Nseg+1):(2*Nseg));
qinC = qin(:,(2*Nseg+1):(2*Nseg+Nt));
PvV  = Pv(:,(Nseg+1):(2*Nseg));


%arterial qout
qoutA             = zeros(N,Nseg);                        %preallocate outlet flow for each arterial segment
qoutA(:,1:Nt)       = qinC(:,1:Nt);                           %outflow of terminal arterioles is the inflow to capillaries - these are the flows that will determine RBC transit times and gas exchange
for i = 1:N
    for k = (Nt+1):Nseg
        qoutA(i,k) =  sum(qinA(i,NETsegA(k,3:4)));
    end
end


%capillary qout
qoutC = qinV(:,1:Nt); %outflow of capillaries is the inflow to terminal venules

%venous qout
qoutV       = zeros(N,Nseg); %preallocate outlet flow for each venous segment
% qoutV(:,xind) = (PvV(:,xind) - Pla + Rd(xind)*qinV(:,xind))./(Rla+Rd(xind)); %flow into left aftrium
qoutV(:,xind) = Qco;
for j = 1:N
    qoutV(j,1:(xind-1)) = (Pv(j,NETsegV(1:xind-1,1))-Pv(j,NETsegV(1:xind-1,2)) + ...
        Rd(NETsegV(1:xind-1,2)).*qinV(j,NETsegV(1:xind-1,3)) + ...
        Rd(NETsegV(1:xind-1,1)).*qinV(j,NETsegV(1:xind-1,1)) - ...
        Rd(NETsegV(1:xind-1,2)).*qinV(j,NETsegV(1:xind-1,2)))...
        ./(Rd(NETsegV(1:xind-1,1))+Rd(NETsegV(1:xind-1,2)));
    qoutV(j,(xind+1):end) = (Pv(j,NETsegV((xind+1):end,1))-Pv(j,NETsegV((xind+1):end,2)) + ... %venous outflows - using NETsegV struct to index partner and confluent vessels
        Rd(NETsegV((xind+1):end,2)).*qinV(j,NETsegV((xind+1):end,3)) + ... %this equation comes from resolving the boundary condition for vessel segments joining together - mathematically it breaks down to solving a system of 2 algebraic equations (pressure and flow conservation)
        Rd(NETsegV((xind+1):end,1)).*qinV(j,NETsegV((xind+1):end,1)) - ...
        Rd(NETsegV((xind+1):end,2)).*qinV(j,NETsegV((xind+1):end,2)))...
        ./(Rd(NETsegV((xind+1):end,1))+Rd(NETsegV((xind+1):end,2)));
end
                  
qout = [qoutA qoutV qoutC]; %stack outflows into one vector

%%%
Pout  = Pv + Ppl + [Rd Rd Rdc].*(qin - qout); %outlet pressures - from conservation of energy/force balance

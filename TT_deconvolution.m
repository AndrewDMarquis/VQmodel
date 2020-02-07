function [TAUf, INDf, hf, t, mu] = TT_deconvolution(Nt, NETsegA,qinA,qinC,qinV,VvA,VvC,VvV)
%%% this function takes segment flows and volumes and calculates the
%%% transit time distribution probability density function - based on the
%%% convolution operation used to interpret indicator-dilution experiments

%local transit times
tauA = VvA(end,:)./qinA(end,:);
tauC = VvC(end,:)./qinC(end,:);
tauV = VvV(end,:)./qinV(end,:);

%%% calculate tau for each unique path
TAU = zeros(Nt,1); %taus for each unique path
IND{Nt} = [];
for i = 1:Nt
   IND{i} = i;
   TAU(i) = tauC(i) + tauA(i) + tauV(i); %initalize the local capillary arterial and venous taus
   par    = NETsegA(i,2); %parent index
   while par ~= 0 %loop to add the taus from the other paths - exit when par == 0
       IND{i} = [IND{i} par];
       TAU(i) = TAU(i) + tauA(par) + tauV(par);
       par = NETsegA(par,2);
   end
end

%%% deconvolution
dt = 0.00005;
t  = 0:dt:200;
h  = zeros(size(t));
f  = qinC(end,:); %capillary flows
F  = sum(f); %total flow
RD = 0.3;
for i = 1:length(TAU)
   h = h + (f(i)/F).*normpdf(t,TAU(i),RD*TAU(i));
end

mu = dt*sum(t.*h); %mean TT

%%% variables to return
TAUf = TAU; %transit times for the unique paths
INDf = IND; %indices of segments for each of the unique paths
hf   = h; %transit time distribution function
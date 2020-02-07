%clear; close all; clc;
%%% This scipt loads the results from quasi-steady state simulations of the
%%% VQ model and plots them. Namely we want to see how HPV affects the
%%% spatial distribution of blood flow and oxygen transport
addpath VascNetwork

FLAG = 1; % select the network we want to plot (1, 2, or 3)

%%% load everything
if FLAG == 1
    load('reuleaux_SS_results.mat','qinC_2','JO2_2','VQ_2','ta','ha','mua',...
        'COORD','Nt','MAP','Element','Node','A','CON','Par','t','qinA_2','qinC_2','qinV_2','VvA_2','VvC_2','VvV_2','POTOT','SOTOT')
    
    [~, ~, ha, ta, mua] = TT_deconvolution(Nt, CON.NETsegA,qinA_2,qinC_2,qinV_2,VvA_2,VvC_2,VvV_2);
    
    qinC_no = qinC_2; JO2_no = JO2_2; VQ_no = VQ_2; ta_no = ta; ha_no = ha; mua_no = mua; POTOT_no = POTOT; SOTOT_no = SOTOT; %no HPV results
    COORD_1 = COORD; Nt_1 = Nt; MAP_1 = MAP; Element_1 = Element; Node_1 = Node; A_1 = A; tno = t; % network strucutre info
    load('reuleaux_HPV_results.mat','qinC_2','JO2_2','VQ_2','ta','ha','mua','POTOT','SOTOT','t')
    qinC_hpv = qinC_2; JO2_hpv = JO2_2; VQ_hpv = VQ_2; ta_hpv = ta; ha_hpv = ha; mua_hpv = mua; POTOT_HPV = POTOT; SOTOT_HPV = SOTOT;  thpv = t;%with HPV results
    load('reuleaux_UVC_results.mat','qinC_2','JO2_2','VQ_2','ta','ha','mua','POTOT','SOTOT','t')
    qinC_uvc = qinC_2; JO2_uvc = JO2_2; VQ_uvc = VQ_2; ta_uvc = ta; ha_uvc = ha; mua_uvc = mua; POTOT_UVC = POTOT; SOTOT_UVC = SOTOT;  tuvc = t;%with UVC results
elseif FLAG == 2
    load('quarter_circ_SS_results.mat','qinC_2','JO2_1','VQ_2','ta','ha','mua',...
        'COORD','Nt','MAP','Element','Node','A','POTOT','SOTOT','t','Par')
    qinC_no = qinC_2; JO2_no = JO2_1; VQ_no = VQ_2; ta_no = ta; ha_no = ha; mua_no = mua;  POTOT_no = POTOT;  SOTOT_no = SOTOT;  %no HPV results
    COORD_1 = COORD; Nt_1 = Nt; MAP_1 = MAP; Element_1 = Element; Node_1 = Node; A_1 = A; tno = t;% network strucutre info
    load('quarter_circ_HPV_results.mat','qinC_2','JO2_2','VQ_2','ta','ha','mua','POTOT','SOTOT','t','CON')
    qinC_hpv = qinC_2; JO2_hpv = JO2_2; VQ_hpv = VQ_2; ta_hpv = ta; ha_hpv = ha; mua_hpv = mua; POTOT_HPV = POTOT;  SOTOT_HPV = SOTOT; thpv = t; %with HPV results
    load('quarter_circ_UVC_results.mat','qinC_2','JO2_2','VQ_2','ta','ha','mua','POTOT','SOTOT','t')
    qinC_uvc = qinC_2; JO2_uvc = JO2_2; VQ_uvc = VQ_2; ta_uvc = ta; ha_uvc = ha; mua_uvc = mua; POTOT_UVC = POTOT; SOTOT_UVC = SOTOT;  tuvc = t;%with UVC results
elseif FLAG == 3
    load('half_circ_SS_results.mat','qinC_2','JO2_1','VQ_2','ta','ha','mua',...
        'COORD','Nt','MAP','Element','Node','A','POTOT','SOTOT','t','Par')
    qinC_no = qinC_2; JO2_no = JO2_1; VQ_no = VQ_2; ta_no = ta; ha_no = ha; mua_no = mua;  POTOT_no = POTOT;  SOTOT_no = SOTOT; %no HPV results
    COORD_1 = COORD; Nt_1 = Nt; MAP_1 = MAP; Element_1 = Element; Node_1 = Node; A_1 = A; tno = t;% network strucutre info
    load('half_circ_HPV_results.mat','qinC_2','JO2_2','VQ_2','ta','ha','mua','POTOT','SOTOT','t','CON')
    qinC_hpv = qinC_2; JO2_hpv = JO2_2; VQ_hpv = VQ_2; ta_hpv = ta; ha_hpv = ha; mua_hpv = mua; POTOT_HPV = POTOT;  SOTOT_HPV = SOTOT; thpv = t; %with HPV results
    load('half_circ_UVC_results.mat','qinC_2','JO2_2','VQ_2','ta','ha','mua','POTOT','SOTOT','t')
    qinC_uvc = qinC_2; JO2_uvc = JO2_2; VQ_uvc = VQ_2; ta_uvc = ta; ha_uvc = ha; mua_uvc = mua; POTOT_UVC = POTOT; SOTOT_UVC = SOTOT;  tuvc = t;%with UVC results
end

%%% average value of oxygen tension and Hb saturation
tL   = 2*Par.R/(2*pi); %"time length" - amount of time to integrate over - 2 respiratory cycles
tINDno = find(tno<(tno(end)-tL), 1, 'last' ); %index of time vector to start integrating over
tINTno = tno(tINDno:end); %relevant portion of the time vector
tINDhpv = find(thpv<(thpv(end)-tL), 1, 'last' ); %index of time vector to start integrating over
tINThpv = thpv(tINDhpv:end); %relevant portion of the time vector
tINDuvc = find(tuvc<(tuvc(end)-tL), 1, 'last' ); %index of time vector to start integrating over
tINTuvc = tuvc(tINDuvc:end); %relevant portion of the time vector


POTOT_no  = trapz(tINTno, POTOT_no(tINDno:end,:))/tL;
POTOT_HPV = trapz(tINThpv, POTOT_HPV(tINDhpv:end,:))/tL;
POTOT_UVC = trapz(tINTuvc, POTOT_UVC(tINDuvc:end,:))/tL;

SOTOT_no  = trapz(tINTno, SOTOT_no(tINDno:end,:))/tL;
SOTOT_HPV = trapz(tINThpv, SOTOT_HPV(tINDhpv:end,:))/tL;
SOTOT_UVC = trapz(tINTuvc, SOTOT_UVC(tINDuvc:end,:))/tL;

%%% CV's
qCV_no  = std(qinC_no)/mean(qinC_no);
qCV_hpv = std(qinC_hpv)/mean(qinC_hpv);
qCV_uvc = std(qinC_uvc)/mean(qinC_uvc);

VQCV_no  = std((sum(qinC_no)/Par.Qvent)*VQ_no)/mean((sum(qinC_no)/Par.Qvent)*VQ_no);
VQCV_hpv = std((sum(qinC_hpv)/Par.Qvent)*VQ_hpv)/mean((sum(qinC_hpv)/Par.Qvent)*VQ_hpv);
VQCV_uvc = std((sum(qinC_uvc)/Par.Qvent)*VQ_uvc)/mean((sum(qinC_uvc)/Par.Qvent)*VQ_uvc);

JCV_no  = std(JO2_no)/mean(JO2_no);
JCV_hpv = std(JO2_hpv)/mean(JO2_hpv);
JCV_uvc = std(JO2_uvc)/mean(JO2_uvc);

%%%
Aref = 6.714326703426993e7*1e-6; %total area of the perfused domain (mm^2) - needed to calculate total o2 mass tranfer across alv-cap boundary
JO2TOT_no  = sum(JO2_no.*(Aref.*CON.Ft));
JO2TOT_hpv = sum(JO2_hpv.*(Aref.*CON.Ft));
JO2TOT_uvc = sum(JO2_uvc.*(Aref.*CON.Ft));

%%% plotting
fh =  findobj('type','figure'); %figure handles
fN = length(fh); %number of figures made

SPACE = 10000;
COORDA = COORD_1; COORDA(:,1) = COORD_1(:,1)-SPACE;
COORDB = COORD_1;
COORDC = COORD_1; COORDC(:,1) = COORD_1(:,1)+SPACE;

COMPARE_3PERFUSION_ZONES(Nt_1, MAP_1, Element_1, Node_1,  log10(qinC_no./(Aref.*CON.Ft)), log10(qinC_hpv./(Aref.*CON.Ft)), log10(qinC_uvc./(Aref.*CON.Ft)),...
    'Capillary Perfusion (ml/s/mm^2)', SPACE, 1, 1, fN+1)
hold on
gplot2(A, COORDA,'k-','linewidth',2)
hold on
gplot2(A, COORDB,'k-','linewidth',2)
hold on
gplot2(A, COORDC,'k-','linewidth',2)
set(gca,'xtick',[],'ytick',[],'linewidth',2)
axis equal
axis off

COMPARE_3PERFUSION_ZONES(Nt, MAP, Element, Node,  log10(VQ_no), log10(VQ_hpv), log10(VQ_uvc),...
    'log_{10}(VQ ratio) (unitless)', SPACE, 0, 3, fN+2)
hold on
gplot2(A, COORDA,'k-','linewidth',2)
hold on
gplot2(A, COORDB,'k-','linewidth',2)
hold on
gplot2(A, COORDC,'k-','linewidth',2)
set(gca,'xtick',[],'ytick',[],'linewidth',2)
axis equal
axis off

COMPARE_3PERFUSION_ZONES(Nt_1, MAP_1, Element_1, Node_1,  1e3*JO2_no, 1e3*JO2_hpv, 1e3*JO2_uvc,...
    'Alveolar-Capillary O_2 flux (mmol/s/mm^2)', SPACE, 0, 2, fN+3)
hold on
gplot2(A, COORDA,'k-','linewidth',2)
hold on
gplot2(A, COORDB,'k-','linewidth',2)
hold on
gplot2(A, COORDC,'k-','linewidth',2)
set(gca,'xtick',[],'ytick',[],'linewidth',2)
axis equal
axis off

fs = 14; %fontsize

[fq1, xq1] = ksdensity(log10(qinC_no));
[fq2, xq2] = ksdensity(log10(qinC_hpv));
[fq3, xq3] = ksdensity(log10(qinC_uvc));
figure(fN+4);
subplot(2,2,1)
plot(xq1,fq1,'k',xq2,fq2,'k--',xq3,fq3,'k:','linewidth',2)
set(gca,'fontsize',fs)
legend('w/o HPV', 'w/ HPV', 'UVC')
xlabel('log_{10}(Capillary Flow)')
ylabel('Probability Density')
% grid on

%figure;
subplot(2,2,2)
plot(ta_no/mua_no,ha_no,'k','linewidth',1.5)
hold on
plot(ta_hpv/mua_hpv,ha_hpv,'k--','linewidth',1.5)
hold on
plot(ta_uvc/mua_uvc,ha_uvc,'k:','linewidth',1.5)
set(gca,'fontsize',fs)
xlabel('$$t/\overline{t}$$','Interpreter','latex')
ylabel('$$h(t/ \overline{t})$$','Interpreter','latex')
legend('w/o HPV','w/ HPV','UVC')
% grid on
% ylim([0 1.1])
xlim([0 3])

[fVQ1, xVQ1] = ksdensity(log10(VQ_no));
[fVQ2, xVQ2] = ksdensity(log10(VQ_hpv));
[fVQ3, xVQ3] = ksdensity(log10(VQ_uvc));
% figure(fN+5);
subplot(2,2,3)
plot(xVQ1,fVQ1,'k',xVQ2,fVQ2,'k--',xVQ3,fVQ3,'k:','linewidth',2)
set(gca,'fontsize',fs)
legend('w/o HPV', 'w/ HPV', 'UVC')
xlabel('log_{10}(V/Q ratio)')
ylabel('Probability Density')
% grid on

[fJO1, xJO1] = ksdensity((1e3*JO2_no));
[fJO2, xJO2] = ksdensity((1e3*JO2_hpv));
[fJO3, xJO3] = ksdensity((1e3*JO2_uvc));
% figure(fN+6);
subplot(2,2,4)
plot(xJO1,fJO1,'k',xJO2,fJO2,'k--',xJO3,fJO3,'k:','linewidth',2)
set(gca,'fontsize',fs)
legend('w/o HPV', 'w/ HPV', 'UVC')
xlabel('Alveolar-Capillary O_2 flux (mmol/s/mm^2)')
ylabel('Probability Density')
% grid on

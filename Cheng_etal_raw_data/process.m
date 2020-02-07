clear; close all; clc;
%%% this script take raw data from Cheng, H. Y., Croft, Q. P., Frise, M. C., Talbot, N. P., Petousi, N., Robbins, P. A., & Dorrington, K. L. (2017). Human hypoxic pulmonary vasoconstriction is unaltered by 8 h of preceding isocapnic hyperoxia. Physiological reports, 5(17), e13396.
%%% and processes it into a useable format

%%% from digitized images of their data
O2 = [0.13361638361638573, 100.64935064935065
1.1390172327672357, 98.3116883116883
2.150739885114888, 99.48051948051948
3.123595154845158, 99.0909090909091
4.133912962037965, 99.48051948051948
5.143528346653349, 99.48051948051948
6.117554320679324, 99.74025974025975
7.091346153846156, 99.87012987012987
8.100727397602398, 99.74025974025975
9.218750000000004, 99.87012987012987
10.077578671328677, 96.23376623376623
11.062609265734267, 62.5974025974026
11.976226898101899, 49.35064935064935
13.021665834165834, 49.22077922077922
14.067104770229768, 49.09090909090909
15.040896603396604, 49.22077922077922
16.05121441058941, 49.6103896103896
17.024772102897103, 49.6103896103896
18.034387487512486, 49.6103896103896
19.00794517982018, 49.6103896103896
20.09459290709291, 52.33766233766234
21.12574925074925, 84.28571428571428
22.131384240759242, 102.07792207792207
23.13889235764236, 100.9090909090909
24.11315247252747, 101.2987012987013
25.156718281718284, 100.12987012987013];

SPAP = [0.011965614222812171, 26.070824868233643
1.9611316362164075, 25.966930078879066
3.9814406454924676, 26.252499765025703
5.966178161127296, 26.343337213421734
7.950264978707713, 26.540368584297934
9.955824832083755, 29.233000513328466
11.926138540845761, 31.677803243368302
13.914888694482803, 31.113778169803275
15.933787857974307, 31.629434687990283
17.91895917231208, 31.6494761880662
19.995662012970577, 28.73157259259795
21.957950445728166, 26.486100366559906
23.943338659417414, 26.470743892475763];

CO = [-0.0045888894887458065, 4.88104518706473
2.0186254926307816, 4.912476380715866
4.006532419154563, 4.983544782162716
5.999622091453866, 5.195583868703773
7.981374507369215, 5.099249581601252
10.027425363062136, 5.751832856448739
12.07218053231118, 6.369173460022674
14.054256869837499, 6.281649840738541
16.02499595097986, 5.88575284781083
18.05274523565297, 6.040533390919397
20.043729417480968, 5.195303136640932
22.018031636344006, 4.896323489715488
24.03930248879771, 4.874890676456297];

VENT = [0.4850746268656714, 8.60651696472592
2.014925373134327, 11.138227780018827
3.9925373134328357, 11.330285509389988
6.0074626865671625, 11.072116893012417
7.98507462686567, 11.309219667428625
9.999999999999998, 11.501501501501505
12.014925373134327, 24.80189144368249
13.992537313432836, 22.201156380260862
16.007462686567166, 20.861906682802207
17.985074626865668, 19.702613060822017
20, 16.87687687687688
22.014925373134325, 11.393483035274084
24.029850746268657, 12.396575680157774];

%%% manipulating things
O2t   = O2(:,1);   O2   = O2(:,2);
SPAPt = SPAP(:,1); SPAP = SPAP(:,2);
COt   = CO(:,1);   CO   = CO(:,2);
VENTt = VENT(:,1); VENT = VENT(:,2);

COr = 80/1000;      %(assumed) resting cardiac output of an adult rat (L/min)
COh = 5;            %(assumed) resting cardiac output of an adult human (L/min) - (Cheng data is around this range)
K   = COr/COh;      %proportionality factor to convert human CO to rat OC per --> COr = K*COh
CO  = K*CO*1000/60; %converting the human data from Cheng et al to rat sized quantities (ml/s)

VENT = K*VENT*1000/60; %converting ventilation to rat valus with K, also dividing ventilation data by 2 to account for airflow being didirectional (ml/s)

T = mean([SPAPt COt VENTt],2); %averaging time so we can have on consistent time vector for SPAP, CO, and VENT

%%% averaging data over normoxic and hypoxic conditions
norm   = find(T<9);                 %indices for the normoxic period
hypo   = find((T>11) & (T<20));     %indices for the normoxic period
normO2 = find(O2t<9);               %indices for the hypoxic period
hypoO2 = find((O2t>11) & (O2t<20)); %indices for the hypoxic period

%normoxic
O2n   = mean(O2(normO2));
SPAPn = mean(SPAP(norm));
COn   = mean(CO(norm));
VENTn = mean(VENT(norm));

%hypoxic
O2h   = mean(O2(hypoO2));
SPAPh = mean(SPAP(hypo));
COh   = mean(CO(hypo));
VENTh = mean(VENT(hypo));

%%% making a struct to put everything in
DATA.T    = T;
DATA.O2t  = O2t;
DATA.O2   = O2;
DATA.SPAP = SPAP;
DATA.CO   = CO;
DATA.VENT = VENT;

DATA.O2n   = O2n;   DATA.O2h = O2h;
DATA.SPAPn = SPAPn; DATA.SPAPh = SPAPh;
DATA.COn   = COn;   DATA.COh = COh;
DATA.VENTn = VENTn; DATA.VENTh = VENTh;

save('Cheng_DATA.mat','DATA')

%%% plots to confirm we got the right info
figure;
subplot(4,1,1)
plot(O2t,O2,'ko')
hold on
plot(O2t,O2,'k')
set(gca,'fontsize',15)
ylabel('O_2 Tension (mmHg)')
xlim([0 25])
grid on

subplot(4,1,2)
plot(T,SPAP,'ko')
hold on
plot(T,SPAP,'k')
set(gca,'fontsize',15)
ylabel('SPAP (mmHg)')
xlim([0 25])
grid on

subplot(4,1,3)
plot(T,CO,'ko')
hold on
plot(T,CO,'k')
set(gca,'fontsize',15)
ylabel('Cardiac Output (L/min)')
xlim([0 25])
grid on

subplot(4,1,4)
plot(T,VENT,'ko')
hold on
plot(T,VENT,'k')
set(gca,'fontsize',15)
xlabel('Time (min)')
ylabel('Ventilation (L/s)')
xlim([0 25])
grid on

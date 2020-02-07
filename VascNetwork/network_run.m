clear; close all; clc;
%%% This script makes a vascular network to cover some defined domain. Made with lungs in mind, but is theoretically generalizable to any tissue
%%% Copyright University of Michigan, Andrew (Drew) Marquis, Feb 11th 2020
%%% If you find a bug, have a question/comment/suggestion, email me at admarqui@umich.edu

addpath PolyMesher

FLAG    = 1;   %set FLAG to select which network you want to construct (1 == reuleaux, 2 == quarter circ, 3 == half circ)
display = 100; %display paremeter - how often do we want to check that NET_GEN algorithm is progressing. For practial use, display << K but really only matters when you're making a huge network
savB    = 0;   %boolean to save network results or not (0 == no, 1 == yes)


if FLAG == 1
    load('SEED_reuleaux_triangle.mat', 'P', 'Element', 'Node') %load info for finite element mesh over domain of interest (made with PolyMesher)
    sc        = 1;              %scaling - can be used to make the lungs bigger/smaller
    root      = [4000 7000]*sc; %root coordinates - where to start the network
    grow_size = 75*sc;          %how far a new node should be placed from exisiting network (with units of um)
    K         = 800;            %number of nodes to make
    
elseif FLAG == 2
    load('SEED_quarter_circ.mat', 'P', 'Element', 'Node') %load info for finite element mesh over domain of interest (made with PolyMesher)
    sc        = 1;             %scaling - can be used to make the lungs bigger/smaller
    root      = [500 -500]*sc; %root coordinates - where to start the network
    grow_size = 75*sc;         %how far a new node should be placed from exisiting network (with units of um)
    K         = 1000;          %number of nodes to make
elseif FLAG == 3
    load('SEED_half_circ.mat', 'P', 'Element', 'Node') %load info for finite element mesh over domain of interest (made with PolyMesher)
    sc        = 1;           %scaling - can be used to make the lungs bigger/smaller
    root      = [0 7000]*sc; %root coordinates - where to start the network
    grow_size = 75*sc;       %how far a new node should be placed from exisiting network (with units of um)
    K         = 1000;        %number of nodes to make
end

P    = P*sc;
Node = Node*sc; %note that "Node" refers to nodes in the finite element graph from polymesher - not the nodes made by the branching algorithm
Xc   = P(:,1);  %x coordinates
Yc   = P(:,2);  %y coordinates

%%%%%% generate the network topology and geometry
tic;
[Net, i3, COORD, A] = NET_GEN(grow_size, K, Xc, Yc, root, display);
toc;

%%%%%% compute anatomical statistics for the network
tic;
[Ft, FLOW, NpO, segL, AS, L, R, MAP, FS, NpOS] = ANATOMY_STATS(Xc, Yc, Node, Element, i3, Net, COORD, K+2, grow_size);
toc;

%%%%%%%%% diamter definied Strahler ordering
tic;
[SO, DO, iter] = DdStrahler_order(Net, segL, R);
toc;

%%%%%%%%% compute mean length and diamter for each strahler order - used to compare to data from Jiang, et al.
MO  = max(SO);      %max Strahler order
LSO = zeros(1,MO);  %mean length per order
DSO = zeros(1,MO);  %mean diamter per order
CSA = zeros(1,MO);  %total cross sectional area for each order
Lse = zeros(1,MO);  %standard deviaion for length per order
Dse = zeros(1,MO);  %standard deviaion for diamter per order
Cstd = zeros(1,MO); %standard deviaion for cross sectional area per order
for i = 1:MO
   clear ind temp
   ind     = find(SO == i);
   temp    = L(ind);
   LSO(i)  = mean(temp);
   Lse(i)  = std(temp)/sqrt(length(temp));
   DSO(i)  = mean(DO{i});
   Dse(i)  = std(DO{i})/sqrt(length(temp));
   CSA(i)  = length(DO{i})*pi*(1e-3*DSO(i)/2)^2;
end

%%%% assemble data struct for kinetic system of ODEs
%NETseg = [segmentIndex radius Length parentIndex daughter1Index daughter2Index]
%the primary vessel has parentIndex == 0. Terminal vessels have daughterIndex == 0.5
Nseg   = length(segL);  %number of segments
NETseg = zeros(Nseg,6); %preallocae look-up table
clear i
for i = 1:Nseg
    clear IND par pari
    IND = find(AS(:,i)==1);
    if length(IND) == 1 %if terminal segment
        NETseg(i,:) = [i R(i) L(i) IND 0.5 0.5];
    elseif length(IND) == 2 %if the primary feed vessel
        NETseg(i,:) = [i R(i) L(i) 0 IND'];
        X = i;
    elseif length(IND) == 3 %if an arbitrary segment
        par_node = Net(segL{i},1);                           %list of parent nodes
        TF       = logical(1 - ismember(par_node,segL{i}));  %logical index of parent node outside of segment
        c        = par_node(TF);                             %extract parent node index
        for j = 1:Nseg
            if ismember(c,segL{j}) %loop to determine which segment that parent is in
                par = j;
            end
        end
        IND(IND == par) = []; %remove parent index from IND
        NETseg(i,:) = [i R(i) L(i) par IND'];
    end
end
%%% Make venous look-up tables of vessel indices
NETsegV = zeros(size(NETseg,1),5);
for i = 1:size(NETsegV,1)
    clear confluent partner upstream
    upstream = NETseg(i,5:6); %one of the upsteam vessels 
    if NETseg(i,4) ~= 0
        confluent = NETseg(i,4);           %find confluent vessel index (topologically equivalent to arterial parent)
        partner   = NETseg(confluent,5:6); %partner vessel(s) - redundant in that this is a vector with i index and its partner
        partner(partner == i) = [];        %delete i index so partner is no longer redundant
        NETsegV(i,:) = [i partner confluent upstream]; %store into venous look-up table
    else
       NETsegV(i,:) = [i 0 0 upstream]; 
    end
end
%%% Make cell arrays that stores indices of terminal vessel perfused by non-terminal segments, and pathlengths to terminals
Nt      = length(i3); %number of terminals
clear IND
IND{Nt} = []; %cell array to store indices for all of the unique paths
for i = 1:Nt
   IND{i} = i;
   par    = NETseg(i,4); %parent index
   while par ~= 0 %loop to add the taus from the other paths - exit when par == 0
       IND{i} = [IND{i} par];
       par    = NETseg(par,4);
   end
end
TI{Nseg-Nt}  = []; %preallocate cell array 'Terminal Indices'
LtT{Nseg-Nt} = []; %perallocate cell array 'Length to Terminals'
for i = 1:Nseg-Nt  %index of non-terminal vessel for array
    vi = Nt+i;     %index of non-terminal vessel for indexing
    for j = 1:Nt   %iterate through all of the terminal segments
        if ismember(vi,IND{j})
            TI{i}  = [TI{i} j]; %store index if the terminal is perfused by the i'th non-terminal segment
        end
    end
    
    LtT{i} = zeros(1,length(TI{i}));   %preallocate zeros for pathlengths
    for k = 1:length(TI{i})            %iterate by distinct path
       LtT{i}(k) = L(TI{i}(k));        %length of terminal segment
       par       = NETseg(TI{i}(k),4); %parent segment
       while par ~= vi
           LtT{i}(k) = LtT{i}(k) + L(par); %add length
           par2      = NETseg(par,4);      %find next parent
           if par2 == vi
              LtT{i}(k) = LtT{i}(k) - L(par)/2; %subtract half so we're only going to the mid point of the relevant segment
           end
           par = par2; %parent is now new parent
       end
    end
end

if savB == 1
    if FLAG == 1
        save('NETWORK_reuleaux.mat','Net','i3','COORD','A','Ft','FLOW','segL', 'AS', 'L', 'R','SO', 'DO', 'MAP', 'FS', 'NETseg', 'NETsegV','TI','LtT','X', 'Node', 'Element','DSO','LSO');
    elseif FLAG == 2
        save('NETWORK_quarter_circ.mat','Net','i3','COORD','A','Ft','FLOW','segL', 'AS', 'L', 'R','SO', 'DO', 'MAP', 'FS', 'NETseg', 'NETsegV','TI','LtT','X', 'Node', 'Element','DSO','LSO');
    elseif FLAG == 3
        save('NETWORK_half_circ.mat','Net','i3','COORD','A','Ft','FLOW','segL', 'AS', 'L', 'R','SO', 'DO', 'MAP', 'FS', 'NETseg', 'NETsegV','TI','LtT','X', 'Node', 'Element','DSO','LSO');
    end
end
%%% NETWORK_[].mat file are accessed by scripts in the VQmodel directory

%%%%% Data from Jiang, et al. - Journal of Applied Physiology 76.2 (1994): 882-892.
%%%%% using their data to validate the anatomical characteristics of the network
D_DATA    = [13.3 31.1 43.7 60.5 89.9 156 262 404 613 938 1638]; % um - diamter measurement for each Strahler order by segment
D_DATA_SE = [3.2 5.3 4.0 5.1 13 24 36 44 73 136 193];            %standard error for diamter measurements
N_DE      = [63961 12789 4098 1540 777 346 146 49 28 13 8];      %number of segments per order

L_DATA    = [0.04 0.11 0.18 0.23 0.29 0.39 0.47 0.63 0.69 1.11 2.26]; %(mm) - length measurements for each Strahler order
L_DATA_SE = [0.02 0.07 0.08 0.09 0.12 0.18 0.21 0.30 0.36 0.65 1.73]; % standard error for length measurements
L_DATA    = L_DATA*1e3;    %convert to um
L_DATA_SE = L_DATA_SE*1e3;

%%%%% Plotting
figure; %visualizing the network itself
plot(Xc,Yc,'.','markersize',0.75) %plot seed points
hold on
ms = 30;
gplot2(A,COORD,'k.-','linewidth',2,'markersize',ms); %overlay vascular network
hold on
gh = plot(COORD(i3,1), COORD(i3,2),'r.','markersize',ms); %overlay terminal nodes
axis equal
set(gh,'linewidth',2)
axis off
scalebar('location','northeast','Unit','\mum')

figure; %hisogram of flow in terminal branches
histogram(Ft/mean(Ft),100)
set(gca,'fontsize',18)
xlabel('Flow through termianl branches (Flow/(mean(Flow)))')

fp = 1:MO;
fp = fp+0.25+(11-MO);
figure; %compare network anatomy to Fung data (diameter by Strahler orders)
h1 = semilogy(D_DATA,'k.','linewidth',1.5);
hold on
semilogy(D_DATA+D_DATA_SE,'k^','linewidth',1.5)
hold on
semilogy(D_DATA-D_DATA_SE,'kv','linewidth',1.5)
hold on
h2 = semilogy(fp,DSO,'bo', 'linewidth',1.5);
% hold on
% semilogy(fp,DSO+Dse,'b^', 'linewidth',1.5)
% hold on
% semilogy(fp,DSO-Dse,'bv', 'linewidth',1.5)
set(gca,'fontsize',18)
xlabel('Strahler Order')
ylabel('Diameter (\mum)')
H = [h1(1) h2(1)];
legend(H, 'Jiang et al.','Model','location','southeast')
grid on

figure; %compare network anatomy to Fung data (length by Strahler orders)
hh1 = semilogy(L_DATA,'k.','linewidth',1.5);
hold on
semilogy(L_DATA+L_DATA_SE,'k^','linewidth',1.5)
hold on
semilogy(L_DATA-L_DATA_SE,'kv','linewidth',1.5)
hold on
hh2 = semilogy(fp,LSO,'bo', 'linewidth',1.5);
% hold on
% semilogy(fp,LSO+Lse,'b^', 'linewidth',1.5)
% hold on
% semilogy(fp,LSO-Lse,'bv', 'linewidth',1.5)
set(gca,'fontsize',18)
xlabel('Strahler Order')
ylabel('Length (\mum)')
HH = [hh1(1) hh2(1)];
legend(HH, 'Jiang et al.','Model','location','southeast')
grid on

figure; %plot to confirm that we collected segments properly
gplot(A,COORD,'k.-')
grid on
for ii = 1:size(segL,2)
    hold on
    gplot(A(segL{ii},segL{ii}),COORD(segL{ii},:),'o--')
end
axis equal

fh =  findobj('type','figure'); %figure handles
fN = length(fh); %number of figures made
Nt = length(i3); %number of terminal nodes

PLOT_PERFUSION_ZONES(Nt, MAP, Element, Node, Ft, 'Perfusion Zone Area', 0, fN+1)
hold on
gplot2(A, COORD,'k-','linewidth',2);
set(gca,'xtick',[],'ytick',[])
axis equal
axis off
scalebar('location','northwest','Unit','\mum')

NETWORK_PLOT(segL, COORD, Net, SO, 'Strahler Order', 0, fN+2)
set(gca,'xtick',[],'ytick',[])
axis equal
axis off

NETWORK_PLOT(segL, COORD, Net, R, 'radius', 0, fN+3)
set(gca,'xtick',[],'ytick',[])
axis equal
axis off


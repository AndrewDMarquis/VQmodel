function [FtF, FLOWF, NpOF, segLF, ASF, L, R, MAP, FS, NpOS] = ANATOMY_STATS(Xc, Yc, Node, Element, i3, Net, COORD, K, grow_size)
%%% This function computates anatomical information for the vascular
%%% network generated with NET_GEN.
%%% last updated: 2/25/19 - Drew Marquis

%(1) assign a normalized flow to every node of the network based on the assumption that flow in the terminal branches is proportional the area of the perfusion zone
%(2) determine vessel segments (discrete portions of network in between bifurcations)
%(3) compute length and radius for each vessel segment

% Xc  and Yc and (x,y) coordiantes of the domain covered by the network
% Node and Element are made from Polymesher - they denote the vertices of the voronoi tesselation used to make the seed points (Xc, Yc) - need these to calculate perfusion zones
% Net and i3 are outputs from NET_GEN, and K is the number of nodes - also used as an input to NET_GEN. Note that this function assumes K = K +2 wrt to the K used as an input to NET_GEN
% grow_size is also an input to NET_GEN that is needed here to calculate
% length of vessel segments

% FtF is the normalized flow in the terminal nodes - assigned based on the
% area of the perfusion zone. This corresponds to the indices in i3
% FLOWF is the flow through each node - computed by summing flows going upstream through the network
% Note: FLOWF(i3) = FtF
% NpOF is the number of perfused outlets each vessel segment gives rise to. Terminal outlets have NpOF(i) = 1
% segLF is a cell array containing the indices of nodes corresponding to the vessel segment they are a part of
% ASF is the adjaceny matrix for segment connectivity
% L and R are the length and radius for each vessel segment in micrometers (um)
% MAP contains the indices for each seed point closest to each terminal branch - used in PLOT_PERFUSION_ZONES
% FS is analagous to FLOWF but for each segment rather than each node

Np = length(Xc); %number of seed points
Nt = length(i3); %number of terminal branches

%%% (1) Determine perfusion zones
SAP = zeros(Np,1); %preallocate vector to store area of the polygon from the voronoi diagram for seed point in P (Xc,Yc) covers
VIt = zeros(Np,1); %Vector of Indices for terminal branches nearest each sample point
for i = 1:Np       %iterate for each seed point
    clear temp D
    temp   = Node(Element{i},:);            %get the verticies for i'th point in P - "Node" and "Element" correspond to the points in P - part of the output from Polymesher
    SAP(i) = polyarea(temp(:,1),temp(:,2)); %compute and save area
    D      = sqrt((COORD(i3,1)-Xc(i)).^2 + (COORD(i3,2)-Yc(i)).^2); %distance between i'th point of P and all terminal branches
    [~, VIt(i)] = min(D); %save the index for the terminal branch closest to each seed point
end

ApTB = zeros(1,Nt); %preallocate area covered by each terminal  - "Area per Terminal Branch"
for h = 1:Nt        %iterate by each terminal node
    clear IND
    IND = VIt == h;          %find all the points covered by terminal branch
    ApTB(h) = sum(SAP(IND)); %sum the area of the points 
end
Ft = ApTB./(sum(ApTB)); % normilizing area perfused by total area - this allows us to assign flow to each terminal branch -
                        % based on the assumption that the flow within a terminal branch is proportional to the area it perfuses --> sum(Ft) == 1
                        
MAP{Nt} = []; %indices for each seed point closest to each terminal branch - used in PLOT_PERFUSION_ZONES
for i = 1:Nt
    for j = 1:Np
       if i == VIt(j)
          MAP{i} = [MAP{i} j]; 
       end
    end
end

%%%% more (1) - calculating flow through the whole network - summing flows by going upstream through the network
FLOW     = zeros(1,K); %preallocate - flow within each node (note that flow is unitless/normalized here)
FLOW(i3) = Ft;         %store terminal node flows
NpO      = zeros(1,K); %preallocae - number of perfused outlets each node could leads to
NpO(i3)  = 1;          %terminal nodes have only one outlet

i3t   = i3;     %temporary terminal indices - should eventually collapse to a scalar as the loop below progresses
count = Nt;     %number of already stored flows
while count < K %end loop once we've computed a flow for every node
    clear i
    for i = 1:length(i3t)      %iterate for each terminal node
        par = Net(i3t(i),1);   %find parent node index
        if par == 1            %if par is the root node
            FLOW(1) = FLOW(2); %assign final flow
            NpO(1)  = NpO(2);  %assign final number of perfused outlets
            count   = count+1; %update count - this should be very last count update to exit the large while loop
            break;             %per the above comment, we don't want to continue through this for-loop, or else the nested while condition throws an error
        end
        while Net(par,2) == 1         %if parent normal, keep going
            FLOW(par) = FLOW(i3t(i)); %store same flow into parent
            NpO(par)  = NpO(i3t(i));  %store same number of terminals into parent        
            i3t(i)    = par;          %update temporary terminal index to keep track of what node we're on             
            count     = count+1;      %successfully stored a flow - increase count
            par       = Net(par,1);   %find next parent
            if par == 0
               break; 
            end
        end %the point of this while loop is to propogate flow upstream until we reach a bifurcation
    end
    %%% handling bifurcations
    bif_par = Net(i3t,1); %indices of parents (bifurcation points)
    ubif    = []; %preallocating - bifurcation nodes we can deal with this iteration - 'useful' bifucation
    dau     = []; %preallocating - indices of daughter nodes that come from the bifurcation point
    for j = 1:length(bif_par)         %iterate by bifuracatin points
        TF = bif_par == bif_par(j);   %'True-False' - logical comparision of every index in bif_par with itself
        if sum(TF) == 2               %if TF == 2, then bif_par(j) is a bifurcation node we can deal with
            ubif = [ubif bif_par(j)]; %add index to list
            dau  = [dau i3t(TF)];     %add indices of daughter nodes to list
        end
    end
    [ubif,ind] = unique(ubif); %by design, we double count the useful bifurcation nodes - so we need to use 'unique' to remove redundancy from the list(s)
    dau        = dau(:,ind);   %remove redundant entries by indexing with output from previous line
    for k = 1:length(ubif)     %iterate by useful bifurcation nodes
        FLOW(ubif(k)) = FLOW(dau(1,k)) + FLOW(dau(2,k)); %summing flows - hurray for flow conservation
        count = count + 1;     %successfully stored a flow - increase count
        ti1 = i3t == dau(1,k); %index of i3t to update with new parent
        ti2 = i3t == dau(2,k); %index of i3t to delete - the choice between ti1 and and ti2 is arbitrary and can be used interchangeably
        i3t(ti1) = ubif(k);    %updating a node -  which node is updated and which is deleted is arbitrary, but it is important to update a node BEFORE you delete the other, otherwise indexing will be offset and that is no good
        i3t(ti2) = [];         %deleting the other node - the above comment applies to this line as well
        
        NpO(ubif(k)) = NpO(dau(1,k))+NpO(dau(2,k)); %adding number of perfused outlets together at bifurcations
    end
end

%%% - "final" version of outputs to ensure that they are not prematurely returned
FLOWF = FLOW;
FtF   = Ft;
NpOF  = NpO;

%%%% (2) determine segments
bInd       = find(Net(:,2) == 2); %find indices of bifucations
segIND     = [i3; bInd];          %stack terminal indices and bifucations indicies
Nseg       = length(segIND);      %total number of segments
segL{Nseg} = [];                  %preallocate cell array of segment indices

clear j
for j = 1:Nseg
    segL{j} = segIND(j);        %store first point (terminal or bifurcation)
    par     = Net(segIND(j),1); %find parent index
    if par == 1
        segL{j} = [segL{j} par]; %the parent index of node 1 is 0, which messes up the while condition. So we need this if statement 
    else
        while Net(par,2) == 1 %see if parent is normal
            segL{j} = [segL{j} par]; %add parent to segment list
            par     = Net(par,1);    %find next parent.
            if par == 0
                break;
            end
        end
    end
end
segLF = segL; %final version of segL to return

%%% flow and NpO for each segment
FS   = zeros(1,Nseg);
NpOS = zeros(1,Nseg);
for i = 1:Nseg
    FS(i) = mean(FLOW(segL{i}));
    NpOS(i) = mean(NpO(segL{i}));
end

%%%%%% (3) calculate length and radius for each segment - also make segment adjaceny matrix
L  = zeros(Nseg,1); 
R  = zeros(Nseg,1);
AS = zeros(Nseg);

Fscale = 75*1e12/60; %average cardiac ouput of an adult rat (75 ml/min - converted to um^3/s)- used to scale unitless/normalized flow to physical units
clear i j
for i = 1:Nseg
    clear f dp
    L(i) = grow_size*length(segL{i});           %store length
    f    = FS(i);                               %normalized flow through a segment
    par  = Net(segL{i},1);                      %parents of nodes in the i'th segment
    TF   =  logical(1 - ismember(par,segL{i})); %find the parent that is outside of this segment
    c = par(TF);                                %extract the parent in segment distinct from i'th segment
    if c == 0
        ref = max(NpO)+1;
    else
        for j = 1:Nseg
            if ismember(c,segL{j})
                AS(i,j) = 1;
                ref = NpO(segL{j}(1));
                break;
            end
        end
    end
    dP   = (L(i)*Fscale*f/NpOS(i))/1e12; %pressure gradient
    R(i) = fzero(@RAD_MIN,[0.01 100000],[],L(i),Fscale*f,dP); %calculate radius - using a root finder as we're using Preis et al.'s "in vivo" nonlinear viscosity law which precludes an analytic derivation
end
ASF = AS'+AS;

